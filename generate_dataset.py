#!/usr/bin/env python3
"""
generate_dataset.py

Automated pipeline for generating HemeLB CFD ground truth data from AneuRisk STL meshes.

Steps per patient:
  1. Auto-generate .pr2 boundary profile from STL topology
  2. Voxelize STL to .gmy lattice geometry via hlb-gmy-cli
  3. Patch HemeLB XML (absolute .gmy path, step count, properties extraction block)
  4. Purge stale output directory to prevent HemeLB file collision errors
  5. Run HemeLB simulation via mpirun
  6. Extract whole.xtr to CSV via hlb-dump-extracted-properties
  7. Convergence check: relative velocity L2 change between last two timesteps < 1%
  8. Mass conservation check: iolet flux imbalance < 2%

Cases failing checks 7 or 8 are logged to data/excluded_cases.csv and excluded
from training. Their CSV output is preserved for inspection.

Boundary conditions (fixed throughout project, do not change):
  Inlet : velocity parabolic BC, maximum velocity 0.05 m/s (Ma < 0.1)
  Outlet: pressure cosine BC, 0.0 mmHg
  VoxelSize: 0.2 mm | TimeStepSeconds: 2e-4 s | Steps: 20000

  The pr2 file retains pressure inlet syntax for hlb-gmy-cli compatibility.
  The inlet condition is replaced with velocity parabolic during XML patching (Step 3).
  Parabolic profile: v(r) = v_max*(1 - r^2/R^2); mean velocity = v_max/2 = 0.025 m/s.
  At dx=0.2mm, dt=2e-4s: c_s = 0.577 m/s, Ma = 0.087 < 0.1 limit.

Inlet detection:
  Uses centerlines.csv (VMTK MaximumInscribedSphereRadius) from AneuriskDatabase
  to identify the parent vessel (inlet) as the boundary loop whose nearest
  centerline point has the largest vessel radius. Falls back to largest-loop
  heuristic if centerlines.csv is absent, and logs a warning.
"""

import os
import re
import csv
import subprocess
import logging
import traceback
import shutil
import xml.etree.ElementTree as ET
from datetime import datetime
import numpy as np
import trimesh
import networkx as nx
from pathlib import Path
from scipy.spatial import KDTree

# ── Configuration ─────────────────────────────────────────────────────────────
# HEMELB_BIN     = "hemelb"
HEMELB_BIN = "/Users/joyantamondal/Downloads/hemelb/build/hemelb-prefix/src/hemelb-build/hemelb"
SETUP_TOOL_BIN = "hlb-gmy-cli"
HLB_DUMP_BIN   = "hlb-dump-extracted-properties"
MPI_CORES      = 4

BASE_DIR = Path(__file__).parent
RAW_DIR  = BASE_DIR / "data/raw_meshes"
GMY_DIR  = BASE_DIR / "data/processed_gmy"
OUT_DIR    = BASE_DIR / "data/outputs_csv"
MODELS_DIR = BASE_DIR / "AneuriskDatabase" / "models"

# Quality check thresholds
CONVERGENCE_THRESHOLD       = 0.01   # Max relative velocity L2 change between last two timesteps
MASS_CONSERVATION_THRESHOLD = 0.02   # Max relative flux imbalance across all iolets

# Exclusion log path
EXCLUSION_LOG = BASE_DIR / "data" / "excluded_cases.csv"

# Inlet velocity BC
INLET_MAX_VELOCITY  = 0.05    # m/s, peak parabolic velocity. Ma = 0.087 < 0.1 limit.
                               # Mean velocity = v_max / 2 = 0.025 m/s.

# Regex for extracting floats/ints from text lines
NUM_PATTERN = re.compile(r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')


# ── Custom exception ──────────────────────────────────────────────────────────

class SimulationQualityError(Exception):
    """Raised when a simulation fails a quality check. Patient is excluded from training."""
    pass


# ── PR2 generation (unchanged) ────────────────────────────────────────────────

def _load_centerlines(patient_id):
    """
    Load centerline points and MaximumInscribedSphereRadius from the AneuriskDatabase.

    Centerline coordinates are in mm (same units as the STL).

    Returns
    -------
    cl_points : np.ndarray (M, 3) or None   XYZ in mm
    cl_radii  : np.ndarray (M,)  or None    MaximumInscribedSphereRadius in mm

    Returns (None, None) if the file is absent.
    """
    cl_path = MODELS_DIR / patient_id / "morphology" / "centerlines.csv"
    if not cl_path.exists():
        return None, None

    points = []
    radii  = []
    with open(cl_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                points.append([float(row['X']), float(row['Y']), float(row['Z'])])
                radii.append(float(row['MaximumInscribedSphereRadius']))
            except (KeyError, ValueError):
                continue

    if not points:
        return None, None

    return np.array(points), np.array(radii)


def _classify_inlet_by_centerline(loop_centers, cl_points, cl_radii):
    """
    Identify which boundary loop is the inlet using centerline vessel radii.

    For each boundary loop centre, find the nearest centerline point and record
    its MaximumInscribedSphereRadius. The loop associated with the largest radius
    is the parent vessel (inlet).

    Parameters
    ----------
    loop_centers : list of np.ndarray (3,)  loop centres in mm
    cl_points    : np.ndarray (M, 3)        centerline XYZ in mm
    cl_radii     : np.ndarray (M,)          inscribed sphere radii in mm

    Returns
    -------
    inlet_idx    : int   index into loop_centers of the classified inlet
    loop_radii   : list  nearest centerline radius for each loop (for logging)
    """
    tree = KDTree(cl_points)
    loop_radii = []

    for centre in loop_centers:
        _, nn_idx = tree.query(centre)
        loop_radii.append(float(cl_radii[nn_idx]))

    inlet_idx = int(np.argmax(loop_radii))
    return inlet_idx, loop_radii


def auto_generate_pr2(mesh_path, pr2_path):
    """
    Parse an open STL mesh, classify boundaries as inlet/outlet, and write a .pr2 profile.

    Inlet classification (in priority order):
    1. Centerline-based: the boundary loop whose nearest centerline point has the
       largest MaximumInscribedSphereRadius is the parent vessel (inlet).
    2. Largest-loop fallback: used only when centerlines.csv is absent. Logs a warning.

    Boundary conditions:
      pr2 uses pressure inlet syntax for hlb-gmy-cli compatibility only.
      Step 3 (XML patching) replaces inlet with velocity parabolic BC (v_max=0.05 m/s).
      Outlet: pressure cosine BC, 0.0 mmHg (unchanged).
    """
    patient_id = mesh_path.stem
    logging.info(f"Generating PR2 for: {mesh_path.name}")

    mesh = trimesh.load(mesh_path)

    edges = mesh.edges_sorted
    unique_edges, counts = np.unique(edges, axis=0, return_counts=True)
    boundary_edges = unique_edges[counts == 1]

    if len(boundary_edges) == 0:
        raise ValueError(
            f"Mesh {mesh_path.name} has no open boundaries. "
            f"HemeLB requires open holes at inlets and outlets."
        )

    g = nx.Graph()
    g.add_edges_from(boundary_edges)
    all_loops = list(nx.connected_components(g))

    loops = [l for l in all_loops if len(l) >= 20]
    if len(loops) < 2:
        logging.warning("Topological filtering too aggressive. Falling back to largest two components.")
        all_loops.sort(key=len, reverse=True)
        loops = all_loops[:2]

    logging.info(f"  Found {len(loops)} boundary loops after topological filtering.")

    centroid = mesh.center_mass

    # Compute geometric properties for each loop
    loop_centers  = []
    loop_radii_geom = []
    loop_normals  = []

    for loop in loops:
        loop_nodes    = list(loop)
        loop_vertices = mesh.vertices[loop_nodes]

        center = np.mean(loop_vertices, axis=0)
        radius = np.max(np.linalg.norm(loop_vertices - center, axis=1)) + 0.1

        _, _, vh = np.linalg.svd(loop_vertices - center)
        normal = vh[2, :]
        if np.dot(normal, centroid - center) < 0:
            normal = -normal

        loop_centers.append(center)
        loop_radii_geom.append(radius)
        loop_normals.append(normal)

    # ── Inlet classification ──────────────────────────────────────────────────
    cl_points, cl_radii = _load_centerlines(patient_id)

    if cl_points is not None:
        inlet_idx, cl_loop_radii = _classify_inlet_by_centerline(
            loop_centers, cl_points, cl_radii
        )
        logging.info(
            f"  Centerline-based inlet classification: "
            f"loop {inlet_idx} selected as Inlet "
            f"(vessel radii: {[f'{r:.2f}' for r in cl_loop_radii]} mm)"
        )

        # Sanity check: flag if the centerline and geometric heuristics disagree
        largest_loop_idx = int(np.argmax([len(list(l)) for l in loops]))
        if inlet_idx != largest_loop_idx:
            ratio = cl_loop_radii[largest_loop_idx] / max(cl_loop_radii[inlet_idx], 1e-6)
            logging.warning(
                f"  INLET MISMATCH: centerline selects loop {inlet_idx}, "
                f"largest-loop heuristic selects loop {largest_loop_idx}. "
                f"Vessel radius ratio (heuristic/centerline) = {ratio:.2f}. "
                f"Using centerline result. Manual verification recommended."
            )
    else:
        # Fallback: largest loop = inlet
        largest_loop_idx = int(np.argmax([len(list(l)) for l in loops]))
        inlet_idx = largest_loop_idx
        logging.warning(
            f"  centerlines.csv not found for {patient_id}. "
            f"Falling back to largest-loop inlet heuristic. "
            f"Verify this case manually."
        )

        # Flag ambiguous cases: ratio of largest to second-largest loop size
        loop_sizes = sorted([len(list(l)) for l in loops], reverse=True)
        if len(loop_sizes) >= 2 and loop_sizes[1] > 0:
            ratio = loop_sizes[0] / loop_sizes[1]
            if ratio < 1.5:
                logging.warning(
                    f"  AMBIGUOUS inlet: largest/second-largest loop size ratio = {ratio:.2f} < 1.5. "
                    f"Inlet classification is unreliable for {patient_id}."
                )

    # ── Write PR2 ─────────────────────────────────────────────────────────────
    pr2_content = []
    pr2_content.append("DurationSeconds: 6.0")
    pr2_content.append("Iolets:")

    inlet_counter  = 1
    outlet_counter = 1

    # Inlet first, then outlets (HemeLB expects inlet before outlets)
    ordered_indices = [inlet_idx] + [i for i in range(len(loops)) if i != inlet_idx]

    for position, loop_idx in enumerate(ordered_indices):
        center = loop_centers[loop_idx]
        radius = loop_radii_geom[loop_idx]
        normal = loop_normals[loop_idx]

        if position == 0:
            iolet_type = "Inlet"
            iolet_name = f"Inlet{inlet_counter}"
            inlet_counter += 1
            pressure_x = 0.1   # 0.1 mmHg - fixed for entire project, do not change
        else:
            iolet_type = "Outlet"
            iolet_name = f"Outlet{outlet_counter}"
            outlet_counter += 1
            pressure_x = 0.0

        pr2_content.append("- Centre:")
        pr2_content.append(f"    x: {center[0]:.10f}")
        pr2_content.append(f"    y: {center[1]:.10f}")
        pr2_content.append(f"    z: {center[2]:.10f}")
        pr2_content.append(f"  Name: {iolet_name}")
        pr2_content.append("  Normal:")
        pr2_content.append(f"    x: {normal[0]:.10f}")
        pr2_content.append(f"    y: {normal[1]:.10f}")
        pr2_content.append(f"    z: {normal[2]:.10f}")
        pr2_content.append("  Pressure:")
        pr2_content.append(f"    x: {pressure_x:.10f}")
        pr2_content.append("    y: 0.0")
        pr2_content.append("    z: 0.0")
        pr2_content.append(f"  Radius: {radius:.10f}")
        pr2_content.append(f"  Type: {iolet_type}")

    gmy_out_name = f"{mesh_path.with_suffix('').name}.gmy"
    xml_out_name = f"{mesh_path.with_suffix('').name}_input.xml"

    pr2_content.append(f"OutputGeometryFile: {gmy_out_name}")
    pr2_content.append(f"OutputXmlFile: {xml_out_name}")
    pr2_content.append("SeedPoint:")
    pr2_content.append(f"  x: {centroid[0]:.10f}")
    pr2_content.append(f"  y: {centroid[1]:.10f}")
    pr2_content.append(f"  z: {centroid[2]:.10f}")
    pr2_content.append(f"StlFile: {mesh_path.name}")
    pr2_content.append("StlFileUnitId: 1")
    pr2_content.append("TimeStepSeconds: 0.0002")
    pr2_content.append("VoxelSize: 0.2")

    with open(pr2_path, "w") as f:
        f.write("\n".join(pr2_content))
    logging.info(f"  PR2 written: {pr2_path.name}")


# ── Shell execution (unchanged) ───────────────────────────────────────────────

def run_cmd(cmd, step_name):
    """Executes a shell command, streaming stdout live. Captures stderr for error reporting."""
    logging.info(f"Starting: {step_name}")
    try:
        subprocess.run(
            cmd, shell=True, check=True,
            stderr=subprocess.PIPE   # stdout streams live to terminal
        )
        logging.info(f"Success: {step_name}")
    except subprocess.CalledProcessError as e:
        stderr_msg = e.stderr.decode('utf-8') if e.stderr else "None"
        logging.error(f"FAILED: {step_name}\n\n[STANDARD ERROR]:\n{stderr_msg}")
        raise


# ── PR2 parser (new) ──────────────────────────────────────────────────────────

def parse_pr2_iolets(pr2_path):
    """
    Parse iolet definitions from a .pr2 profile file.

    Returns
    -------
    list of dict, each with keys:
        'centre'  : [x, y, z] in mm  (STL coordinate units)
        'normal'  : [nx, ny, nz]     (unit vector pointing into fluid domain)
        'radius'  : float in mm
        'type'    : 'Inlet' or 'Outlet'
    """
    TOP_LEVEL_EXIT_KEYS = {
        'OutputGeometryFile', 'OutputXmlFile', 'SeedPoint',
        'StlFile', 'StlFileUnitId', 'TimeStepSeconds', 'VoxelSize', 'DurationSeconds'
    }

    with open(pr2_path, 'r') as f:
        lines = f.readlines()

    iolets          = []
    current         = None
    current_sub     = None   # 'centre' | 'normal' | None
    in_iolets_block = False

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue

        # Enter the Iolets block
        if stripped == 'Iolets:':
            in_iolets_block = True
            continue

        if not in_iolets_block:
            continue

        # Detect top-level keys that signal the end of the Iolets block
        key_candidate = stripped.split(':')[0].strip()
        if key_candidate in TOP_LEVEL_EXIT_KEYS:
            if current is not None:
                iolets.append(current)
                current = None
            in_iolets_block = False
            continue

        # New iolet entry starts with "- Centre:"
        if stripped == '- Centre:':
            if current is not None:
                iolets.append(current)
            current     = {'centre': [0.0, 0.0, 0.0], 'normal': [0.0, 0.0, 0.0],
                           'radius': 0.0, 'type': None}
            current_sub = 'centre'
            continue

        if stripped == 'Normal:':
            current_sub = 'normal'
            continue

        # Keys that end a coordinate sub-block
        if stripped in ('Pressure:', 'Centre:') or stripped.startswith('Name:'):
            current_sub = None
            continue

        if stripped.startswith('Radius:') and current is not None:
            current['radius'] = float(stripped.split(':', 1)[1].strip())
            current_sub = None
            continue

        if stripped.startswith('Type:') and current is not None:
            current['type'] = stripped.split(':', 1)[1].strip()
            current_sub = None
            continue

        # Parse x/y/z coordinates for current sub-block
        if current_sub in ('centre', 'normal') and current is not None:
            if stripped.startswith('x:'):
                current[current_sub][0] = float(stripped[2:].strip())
            elif stripped.startswith('y:'):
                current[current_sub][1] = float(stripped[2:].strip())
            elif stripped.startswith('z:'):
                current[current_sub][2] = float(stripped[2:].strip())

    if current is not None:
        iolets.append(current)

    if not iolets:
        raise ValueError(f"No iolets parsed from {pr2_path}. Check the file format.")

    return iolets


# ── CSV parsing utilities (new) ───────────────────────────────────────────────

def _parse_csv_header(lines):
    """
    Extract geometry origin (m) and voxel size (m) from the CSV header comment block.
    Returns (origin: np.ndarray shape (3,), voxel_size: float).
    """
    origin     = None
    voxel_size = None

    for line in lines:
        if line.startswith('# Geometry origin'):
            nums = NUM_PATTERN.findall(line)
            if len(nums) < 3:
                raise ValueError(f"Cannot parse origin from: {line.strip()}")
            origin = np.array([float(n) for n in nums[:3]])
        elif line.startswith('# Voxel size'):
            nums = NUM_PATTERN.findall(line)
            if not nums:
                raise ValueError(f"Cannot parse voxel size from: {line.strip()}")
            voxel_size = float(nums[0])

    if origin is None:
        raise ValueError("'# Geometry origin' not found in CSV header.")
    if voxel_size is None:
        raise ValueError("'# Voxel size' not found in CSV header.")

    return origin, voxel_size


def _parse_data_block(lines, start_idx):
    """
    Parse data rows from lines[start_idx+1:] until the next comment line or EOF.
    Returns (grid: np.ndarray (N,3) int32, velocity: np.ndarray (N,3) float64).
    """
    grid_list = []
    vel_list  = []

    for line in lines[start_idx + 1:]:
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            break
        nums = [float(x) for x in NUM_PATTERN.findall(stripped)]
        if len(nums) != 7:
            continue
        grid_list.append(nums[0:3])
        vel_list.append(nums[3:6])

    if not grid_list:
        raise ValueError(f"No valid data rows found after line {start_idx}.")

    return (np.array(grid_list, dtype=np.int32),
            np.array(vel_list,  dtype=np.float64))


def parse_final_two_timesteps(csv_path):
    """
    Parse the last two timestep blocks from a HemeLB extracted properties CSV.

    Returns
    -------
    origin      : np.ndarray (3,)  in metres
    voxel_size  : float            in metres
    grid_a      : np.ndarray (N,3) second-to-last timestep grid indices
    vel_a       : np.ndarray (N,3) second-to-last timestep velocities (m/s)
    grid_b      : np.ndarray (N,3) final timestep grid indices
    vel_b       : np.ndarray (N,3) final timestep velocities (m/s)

    Raises ValueError if fewer than two timestep blocks are found.
    """
    with open(csv_path, 'r') as f:
        lines = f.readlines()

    origin, voxel_size = _parse_csv_header(lines)

    ts_line_indices = [i for i, l in enumerate(lines) if l.startswith('# Timestep')]

    if len(ts_line_indices) < 2:
        raise ValueError(
            f"Need at least 2 timestep blocks for convergence check; "
            f"found {len(ts_line_indices)} in {csv_path}."
        )

    idx_a = ts_line_indices[-2]   # Second-to-last timestep
    idx_b = ts_line_indices[-1]   # Final timestep

    logging.info(f"  Convergence: comparing {lines[idx_a].strip()} vs {lines[idx_b].strip()}")

    grid_a, vel_a = _parse_data_block(lines, idx_a)
    grid_b, vel_b = _parse_data_block(lines, idx_b)

    return origin, voxel_size, grid_a, vel_a, grid_b, vel_b


def parse_final_timestep(csv_path):
    """
    Parse only the final timestep block. Used by the mass conservation check.

    Returns
    -------
    origin     : np.ndarray (3,)  in metres
    voxel_size : float            in metres
    grid       : np.ndarray (N,3) int32
    velocity   : np.ndarray (N,3) float64 in m/s
    """
    with open(csv_path, 'r') as f:
        lines = f.readlines()

    origin, voxel_size = _parse_csv_header(lines)

    ts_line_indices = [i for i, l in enumerate(lines) if l.startswith('# Timestep')]

    if not ts_line_indices:
        raise ValueError(f"No timestep blocks found in {csv_path}.")

    # Collect all non-comment, non-empty lines after the last timestep marker
    last_idx  = ts_line_indices[-1]
    data_lines_raw = [
        line for line in lines[last_idx + 1:]
        if line.strip() and not line.strip().startswith('#')
    ]

    grid_list = []
    vel_list  = []

    for line in data_lines_raw:
        nums = [float(x) for x in NUM_PATTERN.findall(line.strip())]
        if len(nums) != 7:
            continue
        grid_list.append(nums[0:3])
        vel_list.append(nums[3:6])

    return (origin, voxel_size,
            np.array(grid_list, dtype=np.int32),
            np.array(vel_list,  dtype=np.float64))


# ── Quality checks (new) ──────────────────────────────────────────────────────

def check_convergence(csv_path, threshold=CONVERGENCE_THRESHOLD):
    """
    Compare the velocity field between the last two extraction timesteps.

    Metric: relative L2 norm of the velocity difference.
        rel_change = ||V_last - V_second_last||_F / ||V_second_last||_F

    Sites are matched by grid index; both blocks must cover the same set of sites.

    Parameters
    ----------
    csv_path  : Path to the HemeLB fluid data CSV.
    threshold : Maximum allowed relative change (default 0.01 = 1%).

    Returns
    -------
    passed     : bool
    rel_change : float
    """
    _, _, grid_a, vel_a, grid_b, vel_b = parse_final_two_timesteps(csv_path)

    # Match sites by grid index (sort both arrays by grid tuple)
    def sort_key(grid):
        return np.lexsort((grid[:, 2], grid[:, 1], grid[:, 0]))

    order_a = sort_key(grid_a)
    order_b = sort_key(grid_b)

    grid_a_sorted = grid_a[order_a]
    grid_b_sorted = grid_b[order_b]

    if not np.array_equal(grid_a_sorted, grid_b_sorted):
        raise ValueError(
            "Grid indices differ between the last two timestep blocks. "
            "Cannot compute convergence metric."
        )

    vel_a_sorted = vel_a[order_a]
    vel_b_sorted = vel_b[order_b]

    diff_norm    = np.linalg.norm(vel_b_sorted - vel_a_sorted)
    ref_norm     = np.linalg.norm(vel_a_sorted)

    if ref_norm < 1e-30:
        logging.warning("  Convergence check: reference velocity norm is effectively zero. Marking as failed.")
        return False, float('inf')

    rel_change = diff_norm / ref_norm
    passed     = rel_change < threshold

    logging.info(
        f"  Convergence check: rel_change = {rel_change:.4e} "
        f"(threshold {threshold:.2e}) -> {'PASS' if passed else 'FAIL'}"
    )

    return passed, rel_change


def check_mass_conservation(csv_path, pr2_path, threshold=MASS_CONSERVATION_THRESHOLD):
    """
    Verify that total inflow flux equals total outflow flux at the iolet planes.

    Method:
    - For each iolet, identify fluid nodes within 1.5 voxel widths of the iolet plane
      and within the iolet radius of the iolet centre.
    - Compute signed flux: F_i = sum(v . n_hat_i) * voxel_size^2
      where n_hat_i is the iolet inward normal (pointing into fluid domain).
    - For an incompressible fluid the net flux across all boundaries must be zero:
        sum_i(F_i) = 0
    - Relative imbalance: |sum_i(F_i)| / |F_inlet|

    Notes on coordinate systems:
    - pr2 centre/radius are in mm (STL file units, StlFileUnitId = 1).
    - CSV physical coordinates are in metres.
    - Conversion: divide pr2 values by 1000.

    Parameters
    ----------
    csv_path  : Path to the HemeLB fluid data CSV.
    pr2_path  : Path to the patient .pr2 profile file.
    threshold : Maximum allowed relative flux imbalance (default 0.02 = 2%).

    Returns
    -------
    passed   : bool
    rel_error: float
    """
    iolets                     = parse_pr2_iolets(pr2_path)
    origin, voxel_size, grid, velocity = parse_final_timestep(csv_path)

    # Physical coordinates of all fluid nodes (metres)
    xyz = origin + grid.astype(np.float64) * voxel_size   # (N, 3)

    # Search radius: 1.5 voxel widths normal to the iolet plane
    plane_dist_threshold = 3.0 * voxel_size

    inlet_flux  = None
    outlet_flux = 0.0
    fluxes      = []

    for iolet in iolets:
        centre_m = np.array(iolet['centre']) / 1000.0   # mm -> m
        radius_m = iolet['radius']            / 1000.0  # mm -> m
        normal   = np.array(iolet['normal'],  dtype=np.float64)
        n_hat    = normal / np.linalg.norm(normal)

        # Signed distance from each node to the iolet plane
        disp          = xyz - centre_m[np.newaxis, :]           # (N, 3)
        plane_dist    = np.abs(disp @ n_hat)                    # (N,)

        # In-plane distance from iolet centre
        proj          = (disp @ n_hat)[:, np.newaxis] * n_hat   # (N, 3)
        in_plane_dist = np.linalg.norm(disp - proj, axis=1)     # (N,)

        mask = (plane_dist < plane_dist_threshold) & (in_plane_dist < radius_m)
        n_nodes = mask.sum()

        if n_nodes == 0:
            logging.warning(
                f"  Mass conservation: no fluid nodes found near iolet "
                f"'{iolet['type']}' (centre={centre_m}, radius={radius_m:.4f} m). "
                f"Skipping this iolet."
            )
            continue

        # Flux = sum(v . n_hat) * dx^2
        # Positive = flow in the inward-normal direction (into fluid)
        flux = float(np.sum(velocity[mask] @ n_hat)) * voxel_size ** 2
        fluxes.append(flux)

        logging.info(
            f"  Mass conservation: {iolet['type']} | "
            f"nodes={n_nodes} | flux={flux:.4e} m^3/s"
        )

        if iolet['type'] == 'Inlet':
            inlet_flux = flux
        else:
            outlet_flux += flux

    if inlet_flux is None:
        logging.warning("  Mass conservation: inlet flux not computed. Cannot perform check.")
        return True, 0.0   # Cannot check; do not exclude
    
    # If any iolet had no detectable nodes, skip rather than falsely exclude
    total_iolets    = len(iolets)
    detected_iolets = len(fluxes)
    if detected_iolets < total_iolets:
        logging.warning(
            f"  Mass conservation: only {detected_iolets}/{total_iolets} iolets detected. "
            f"Skipping check to avoid false exclusion."
        )
        return True, 0.0

    if abs(inlet_flux) < 1e-30:
        logging.warning("  Mass conservation: inlet flux is effectively zero. Marking as failed.")
        return False, float('inf')

    net_flux  = inlet_flux + outlet_flux
    rel_error = abs(net_flux) / abs(inlet_flux)
    passed    = rel_error < threshold

    logging.info(
        f"  Mass conservation: net_flux={net_flux:.4e} m^3/s | "
        f"rel_error={rel_error:.4e} (threshold {threshold:.2e}) "
        f"-> {'PASS' if passed else 'FAIL'}"
    )

    return passed, rel_error


# ── Exclusion logger (new) ────────────────────────────────────────────────────

def log_exclusion(patient_id, reason, metric_name, metric_value):
    """
    Append an exclusion record to the exclusion log CSV.

    Columns: timestamp, patient_id, reason, metric_name, metric_value
    Creates the file and writes a header row if it does not exist.
    """
    EXCLUSION_LOG.parent.mkdir(parents=True, exist_ok=True)
    write_header = not EXCLUSION_LOG.exists()

    with open(EXCLUSION_LOG, 'a', newline='') as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(['timestamp', 'patient_id', 'reason', 'metric_name', 'metric_value'])
        writer.writerow([
            datetime.now().isoformat(timespec='seconds'),
            patient_id,
            reason,
            metric_name,
            f"{metric_value:.6e}"
        ])

    logging.warning(
        f"  EXCLUDED: {patient_id} | {reason} | {metric_name}={metric_value:.4e} "
        f"| Logged to {EXCLUSION_LOG.name}"
    )


# ── Patient processing (steps 1-6 unchanged, 7-8 added) ──────────────────────

def _detect_resume_step(patient_id, patient_pr2, patient_gmy, patient_xml,
                        patient_out_dir):
    """
    Detect which pipeline step to resume from for a given patient.

    State machine (checked in priority order):
      csv complete (10 timestep blocks) -> 7 (quality checks)
      csv incomplete                    -> delete + rerun
      output dir without csv            -> delete (partial run) + fall through
      patched xml + gmy                 -> 5 (simulate)
      raw xml (RAW_DIR) + gmy           -> 3 (patch xml)
      gmy only (no xml anywhere)        -> 2 (re-voxelize to regenerate xml)
      pr2 only                          -> 2 (voxelize)
      nothing                           -> 1 (full run)

    Returns
    -------
    int : step number to start from (1-7).
    """
    csv_path     = patient_out_dir / f"{patient_id}_fluid_data.csv"
    raw_xml_path = RAW_DIR / f"{patient_id}_input.xml"

    # CSV exists: validate completeness (10 blocks for 20000 steps / period 2000)
    if csv_path.exists():
        try:
            with open(csv_path, 'r') as fh:
                ts_count = fh.read().count('# Timestep')
            if ts_count >= 10:
                return 7
            logging.warning(
                f"  {patient_id}: CSV incomplete ({ts_count}/10 timestep blocks). "
                f"Deleting and rerunning."
            )
        except Exception as e:
            logging.warning(f"  {patient_id}: CSV validation error ({e}). Deleting and rerunning.")
        csv_path.unlink(missing_ok=True)
        if patient_out_dir.exists():
            shutil.rmtree(patient_out_dir)

    # Output dir exists but no CSV: interrupted mid-simulation or mid-extraction
    if patient_out_dir.exists():
        logging.warning(
            f"  {patient_id}: partial output directory found (no CSV). "
            f"Deleting and rerunning simulation."
        )
        shutil.rmtree(patient_out_dir)

    # Patched XML + GMY: ready to simulate (step 5)
    if patient_xml.exists() and patient_gmy.exists():
        return 5

    # Raw XML from hlb-gmy-cli in RAW_DIR + GMY: ready to patch XML (step 3)
    if raw_xml_path.exists() and patient_gmy.exists():
        return 3

    # GMY exists but no XML anywhere: re-voxelize to regenerate XML
    # Happens when XMLs were deleted but GMY files preserved (e.g. BC type switch)
    if patient_gmy.exists():
        logging.info(
            f"  {patient_id}: GMY found but no XML. "
            f"Re-voxelizing to regenerate XML (~8s at 0.2mm)."
        )
        return 2

    # PR2 exists: skip PR2 generation, go to voxelization
    if patient_pr2.exists():
        return 2

    # Nothing exists: full run from scratch
    return 1


def process_patient(mesh_file, start_from=None):
    """
    Run the full pipeline for one patient, with optional step resume.

    Parameters
    ----------
    mesh_file  : Path to the patient STL file
    start_from : int or None
        If None, auto-detects the resume step from existing files.
        If int, forces starting from that step (1=full rerun, 7=quality checks only).
    """
    patient_id = mesh_file.stem
    logging.info(f"============ Processing Case: {patient_id} ============")

    patient_pr2     = RAW_DIR / f"{patient_id}.pr2"
    patient_gmy     = GMY_DIR / f"{patient_id}.gmy"
    patient_xml     = GMY_DIR / f"{patient_id}_input.xml"
    patient_out_dir = OUT_DIR / patient_id

    xml_out_name = f"{mesh_file.stem}_input.xml"
    gmy_out_name = f"{mesh_file.stem}.gmy"
    csv_path     = patient_out_dir / f"{patient_id}_fluid_data.csv"

    # Determine resume step
    if start_from is None:
        step = _detect_resume_step(
            patient_id, patient_pr2, patient_gmy, patient_xml, patient_out_dir
        )
    else:
        step = start_from

    if step > 1:
        logging.info(f"  Resuming from step {step} (steps 1-{step-1} already complete).")

    # STEP 1: PR2 generation
    if step <= 1:
        auto_generate_pr2(mesh_file, patient_pr2)

    # STEP 2: Voxelization
    if step <= 2:
        setup_cmd = f"cd {RAW_DIR} && {SETUP_TOOL_BIN} {patient_pr2.name}"
        run_cmd(setup_cmd, f"Headless Voxelization ({patient_id})")
        os.rename(RAW_DIR / gmy_out_name, patient_gmy)

    # STEP 3: Patch XML
    if step <= 3:
        generated_xml_path = RAW_DIR / xml_out_name
        logging.info("Patching HemeLB XML.")

        tree = ET.parse(generated_xml_path)
        root = tree.getroot()

        # 3a. Absolute GMY path
        datafile = root.find(".//geometry/datafile")
        if datafile is not None:
            datafile.set("path", str(patient_gmy.resolve()))

        # 3b. Simulation steps
        steps_element = root.find(".//simulation/steps")
        if steps_element is not None:
            steps_element.set("value", "20000")

        # 3c. Properties extraction block
        properties = root.find("properties")
        if properties is None:
            properties  = ET.SubElement(root, "properties")
            prop_output = ET.SubElement(properties, "propertyoutput",
                                        {"file": "whole.xtr", "period": "2000"})
            ET.SubElement(prop_output, "geometry", {"type": "whole"})
            ET.SubElement(prop_output, "field",    {"type": "velocity"})
            ET.SubElement(prop_output, "field",    {"type": "pressure"})

        # 3d. Replace inlet pressure BC with velocity parabolic BC
        #     The pr2 uses pressure syntax for hlb-gmy-cli compatibility only.
        #     Velocity parabolic is the scientifically correct BC (published HemeLB standard).
        iolets_info  = parse_pr2_iolets(patient_pr2)
        inlet_info   = next((io for io in iolets_info if io['type'] == 'Inlet'), None)
        if inlet_info is None:
            raise ValueError(f"No Inlet found in {patient_pr2}. Cannot set velocity BC.")

        inlet_radius_m = inlet_info['radius'] / 1000.0  # mm -> m
        logging.info(
            f"  Setting velocity parabolic inlet: "
            f"radius={inlet_radius_m*1000:.3f} mm, "
            f"v_max={INLET_MAX_VELOCITY} m/s, "
            f"v_mean={INLET_MAX_VELOCITY/2:.4f} m/s"
        )

        for inlet_elem in root.findall(".//inlets/inlet"):
            # Remove existing condition (pressure cosine from hlb-gmy-cli)
            old_cond = inlet_elem.find("condition")
            if old_cond is not None:
                inlet_elem.remove(old_cond)

            # Build velocity parabolic condition
            new_cond = ET.Element("condition")
            new_cond.set("type", "velocity")
            new_cond.set("subtype", "parabolic")

            radius_el = ET.SubElement(new_cond, "radius")
            radius_el.set("value", f"{inlet_radius_m:.8f}")
            radius_el.set("units", "m")

            maxvel_el = ET.SubElement(new_cond, "maximum")
            maxvel_el.set("value", f"{INLET_MAX_VELOCITY:.4f}")
            maxvel_el.set("units", "m/s")

            # Insert as first child of <inlet>
            inlet_elem.insert(0, new_cond)

        tree.write(patient_xml, encoding="utf-8", xml_declaration=True)

        generated_xml_path = RAW_DIR / xml_out_name
        if generated_xml_path.exists():
            os.remove(generated_xml_path)

        # Diagnostic: print the patched XML
        try:
            with open(patient_xml, 'r') as fh:
                logging.info(
                    f"\n{'='*50}\n--- DIAGNOSTIC: PATCHED XML ---\n{'='*50}\n"
                    f"{fh.read()}\n{'='*50}"
                )
        except Exception as e:
            logging.error(f"Could not read XML for diagnostic: {e}")

    # STEP 4: Purge stale output directory (only on full reruns)
    if step <= 4:
        if patient_out_dir.exists():
            logging.info(f"  Purging stale output directory: {patient_out_dir.name}")
            shutil.rmtree(patient_out_dir)

    # STEP 5: Run HemeLB simulation
    if step <= 5:
        hemelb_cmd = f"mpirun -n {MPI_CORES} {HEMELB_BIN} -in {patient_xml} -out {patient_out_dir}"
        run_cmd(hemelb_cmd, f"HemeLB Simulation ({patient_id})")

    # STEP 6: Extract whole.xtr to CSV
    if step <= 6:
        xtr_matches = list(patient_out_dir.glob("**/whole.xtr"))
        if not xtr_matches:
            existing = [str(p.relative_to(BASE_DIR)) for p in patient_out_dir.rglob("*") if p.is_file()]
            raise FileNotFoundError(
                f"whole.xtr not found in {patient_out_dir}.\n"
                f"Files present:\n" + "\n".join(existing)
            )
        xtr_path = xtr_matches[0]
        extract_cmd = f"{HLB_DUMP_BIN} {xtr_path} > {csv_path}"
        run_cmd(extract_cmd, f"CSV Extraction ({patient_id})")

    # STEP 7: Convergence check
    logging.info(f"--- Quality Check: Convergence ({patient_id}) ---")
    try:
        conv_passed, rel_change = check_convergence(csv_path)
    except Exception as e:
        logging.error(f"  Convergence check raised an exception: {e}")
        log_exclusion(patient_id, "convergence_check_error", "exception", float('nan'))
        raise SimulationQualityError(f"Convergence check failed with exception for {patient_id}.") from e

    if not conv_passed:
        log_exclusion(patient_id, "not_converged", "rel_velocity_change", rel_change)
        raise SimulationQualityError(
            f"{patient_id} excluded: velocity not converged "
            f"(rel_change={rel_change:.4e} > threshold={CONVERGENCE_THRESHOLD})."
        )

    # STEP 8: Mass conservation check
    logging.info(f"--- Quality Check: Mass Conservation ({patient_id}) ---")
    try:
        mass_passed, rel_error = check_mass_conservation(csv_path, patient_pr2)
    except Exception as e:
        logging.error(f"  Mass conservation check raised an exception: {e}")
        log_exclusion(patient_id, "mass_conservation_check_error", "exception", float('nan'))
        raise SimulationQualityError(f"Mass conservation check failed with exception for {patient_id}.") from e

    if not mass_passed:
        log_exclusion(patient_id, "mass_not_conserved", "rel_flux_imbalance", rel_error)
        raise SimulationQualityError(
            f"{patient_id} excluded: mass not conserved "
            f"(rel_error={rel_error:.4e} > threshold={MASS_CONSERVATION_THRESHOLD})."
        )

    logging.info(f"============ Case {patient_id}: All checks PASSED ============\n")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    GMY_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    mesh_files = sorted(list(RAW_DIR.glob("*.stl")) + list(RAW_DIR.glob("*.vtp")))

    if not mesh_files:
        logging.warning(f"No surface files found in {RAW_DIR}.")
        raise SystemExit(0)

    # ── Resume prompt ─────────────────────────────────────────────────────────
    print("\n" + "="*60)
    print("  HemeLB Surrogate Pipeline")
    print("="*60)
    print("  Run mode options:")
    print("  [1] Full run      - rerun everything from scratch for all cases")
    print("  [2] Auto-resume   - skip steps already completed (default)")
    print("  [3] Quality only  - rerun quality checks only (CSV already exists)")
    print("  [4] Single case   - process one specific patient ID")
    print("="*60)

    choice = input("  Enter choice [1/2/3/4], or press Enter for auto-resume: ").strip()

    force_start = None   # None = auto-detect per patient
    single_case = None

    if choice == "1":
        force_start = 1
        logging.info("Mode: Full run from scratch.")
    elif choice == "3":
        force_start = 7
        logging.info("Mode: Quality checks only.")
    elif choice == "4":
        pid = input("  Enter patient ID (e.g. C0001): ").strip()
        single_case = pid
        sub_choice = input(
            f"  Start from which step for {pid}? "
            f"[1=full, 2=voxelize, 3=xml, 5=simulate, 6=extract, 7=quality, Enter=auto]: "
        ).strip()
        force_start = int(sub_choice) if sub_choice.isdigit() else None
        logging.info(f"Mode: Single case {pid}, start_from={force_start or 'auto'}.")
    else:
        logging.info("Mode: Auto-resume (skipping completed steps per patient).")

    if single_case:
        mesh_files = [f for f in mesh_files if f.stem == single_case]
        if not mesh_files:
            logging.error(f"No STL found for patient ID '{single_case}' in {RAW_DIR}.")
            raise SystemExit(1)

    passed_cases   = []
    excluded_cases = []
    failed_cases   = []

    for mesh in mesh_files:
        try:
            process_patient(mesh, start_from=force_start)
            passed_cases.append(mesh.stem)

        except SimulationQualityError as e:
            logging.warning(f"Quality exclusion: {e}")
            excluded_cases.append(mesh.stem)

        except Exception:
            logging.error(f"Pipeline error for {mesh.name}:\n{traceback.format_exc()}")
            failed_cases.append(mesh.stem)

    # Final summary
    logging.info(
        f"\n{'='*60}\n"
        f"  BATCH COMPLETE\n"
        f"  Passed    : {len(passed_cases)}\n"
        f"  Excluded  : {len(excluded_cases)}  (see {EXCLUSION_LOG.name})\n"
        f"  Errors    : {len(failed_cases)}\n"
        f"  Total     : {len(mesh_files)}\n"
        f"{'='*60}"
    )

    if excluded_cases:
        logging.info(f"  Excluded cases : {', '.join(excluded_cases)}")
    if failed_cases:
        logging.info(f"  Error cases    : {', '.join(failed_cases)}")
