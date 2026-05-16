#!/usr/bin/env python3
"""
compute_wss.py

Computes Wall Shear Stress (WSS) magnitude at wall-adjacent fluid nodes
from HemeLB extracted properties CSV output (produced by hlb-dump-extracted-properties).

Method
------
Wall-adjacent nodes:
    Fluid lattice nodes that have at least one absent 6-connected face neighbor.
    Absent neighbor = voxel not present in the fluid site list = solid wall voxel.

Wall normal estimation:
    For each wall-adjacent node, accumulate unit vectors pointing toward each absent
    neighbor. The mean of these vectors points into the solid. Negate to get the
    outward wall normal pointing into the fluid domain.

WSS computation:
    Under the standard mid-link bounce-back boundary condition, the physical wall
    surface sits at the midpoint between the fluid node and its solid neighbor,
    i.e., at distance 0.5 * voxel_size from the fluid node.

    The no-slip condition gives v_wall = 0. The tangential velocity gradient is:

        dv_tangential / dn = v_tangential_fluid / (0.5 * voxel_size)

    WSS magnitude:
        WSS = mu * |v_tangential| / (0.5 * voxel_size)

    where v_tangential = v - (v . n_hat) * n_hat

Timestep selection:
    The script uses ONLY the final timestep in the CSV. For fixed-pressure-BC
    steady-state LBM simulations this is the most converged state available.
    TAWSS averaging is not applicable here (no pulsatile waveform).

Limitations:
    - No sub-voxel wall distance correction (would require .gmy link cut distances).
      At 0.1 mm voxel size the error is bounded to half a voxel spacing.
    - Newtonian blood assumption: mu = 0.0035 Pa.s.
    - WSS accuracy is limited by simulation convergence (1000-step prototype runs
      may not be fully converged; final training data requires convergence checks).

Output
------
Per-wall-adjacent-node CSV: node_id, x, y, z, wss_magnitude
x, y, z are in metres. wss_magnitude is in Pascals.
"""

import re
import logging
import traceback
import numpy as np
import pandas as pd
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - [%(levelname)s] - %(message)s'
)

# ── Physical constants ──────────────────────────────────────────────────────
DYNAMIC_VISCOSITY = 0.0035          # Pa.s, Newtonian blood approximation
LOW_WSS_THRESHOLD = 0.4             # Pa, clinical rupture-risk threshold

# ── 6-connected face neighbors (voxel lattice) ──────────────────────────────
FACE_NEIGHBORS = np.array([
    [ 1,  0,  0],
    [-1,  0,  0],
    [ 0,  1,  0],
    [ 0, -1,  0],
    [ 0,  0,  1],
    [ 0,  0, -1],
], dtype=np.int32)

# ── Number extraction regex ─────────────────────────────────────────────────
NUM_PATTERN = re.compile(r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?')


# ── CSV parser ───────────────────────────────────────────────────────────────

def parse_csv(csv_path):
    """
    Parse a HemeLB extracted properties CSV produced by hlb-dump-extracted-properties.

    The file contains a comment header followed by one or more timestep blocks:

        # Timestep 100
        [ i  j  k], [ vx  vy  vz], pressure
        ...
        # Timestep 200
        ...

    Only the FINAL timestep block is returned.

    Parameters
    ----------
    csv_path : str or Path

    Returns
    -------
    origin      : np.ndarray (3,)  geometry origin in metres
    voxel_size  : float            voxel edge length in metres
    grid        : np.ndarray (N,3) integer voxel indices [i, j, k]
    velocity    : np.ndarray (N,3) velocity components in m/s
    pressure    : np.ndarray (N,)  pressure in Pa
    """
    csv_path = Path(csv_path)

    with open(csv_path, 'r') as fh:
        lines = fh.readlines()

    origin = None
    voxel_size = None
    last_ts_line = None  # Index of the last "# Timestep N" line

    for i, line in enumerate(lines):
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

        elif line.startswith('# Timestep'):
            last_ts_line = i  # Overwrite each time; ends on the last marker

    if origin is None:
        raise ValueError(f"'# Geometry origin' not found in {csv_path}")
    if voxel_size is None:
        raise ValueError(f"'# Voxel size' not found in {csv_path}")
    if last_ts_line is None:
        raise ValueError(f"No '# Timestep' markers found in {csv_path}")

    ts_label = lines[last_ts_line].strip()
    logging.info(f"  Selected timestep: {ts_label}")

    # Collect all data lines after the final timestep marker.
    # Skip blank lines and any remaining comment lines (e.g., a trailing newline).
    data_lines = [
        line.strip()
        for line in lines[last_ts_line + 1:]
        if line.strip() and not line.strip().startswith('#')
    ]

    logging.info(f"  Parsing {len(data_lines)} site records...")

    grid_list = []
    vel_list  = []
    pres_list = []

    for line in data_lines:
        nums = [float(x) for x in NUM_PATTERN.findall(line)]
        if len(nums) != 7:
            logging.warning(f"  Skipping malformed line (expected 7 values, got {len(nums)}): {line}")
            continue
        grid_list.append(nums[0:3])
        vel_list.append(nums[3:6])
        pres_list.append(nums[6])

    if not grid_list:
        raise ValueError(f"No valid data rows parsed from final timestep in {csv_path}")

    grid     = np.array(grid_list, dtype=np.int32)
    velocity = np.array(vel_list,  dtype=np.float64)
    pressure = np.array(pres_list, dtype=np.float64)

    logging.info(f"  Parsed {len(grid)} fluid sites.")
    return origin, voxel_size, grid, velocity, pressure


# ── WSS computation ──────────────────────────────────────────────────────────

def compute_wss(csv_path, output_path):
    """
    Compute WSS at every wall-adjacent fluid node and write output CSV.

    Parameters
    ----------
    csv_path    : str or Path  HemeLB fluid data CSV
    output_path : str or Path  Destination for the WSS output CSV

    Returns
    -------
    df : pd.DataFrame  columns: node_id, x, y, z, wss_magnitude
    """
    origin, voxel_size, grid, velocity, _ = parse_csv(csv_path)
    N = len(grid)

    # ── Build fluid-site lookup (O(1) membership test) ───────────────────────
    # Encode each (i,j,k) triple as a single int64 for fast set lookup.
    # Strides are computed from the actual domain extent so the encoding is
    # injective (no collisions) without wasting memory.
    logging.info(f"  Building neighbor lookup for {N} fluid sites...")

    i_stride = int(grid[:, 1].max()) + 2   # j_max
    j_stride = int(grid[:, 2].max()) + 2   # k_max

    def encode(ijk):
        """Encode integer array of shape (M, 3) to int64 array of shape (M,)."""
        return (ijk[:, 0].astype(np.int64) * i_stride * j_stride
                + ijk[:, 1].astype(np.int64) * j_stride
                + ijk[:, 2].astype(np.int64))

    fluid_encoded = set(encode(grid).tolist())

    # ── Identify wall-adjacent nodes and accumulate wall normal ───────────────
    # For each of the 6 face directions, shift the entire grid and check which
    # shifted sites are absent from the fluid set.
    is_wall_adjacent = np.zeros(N, dtype=bool)
    normal_sum       = np.zeros((N, 3), dtype=np.float64)

    for n_dir in FACE_NEIGHBORS:
        shifted = grid + n_dir                         # (N, 3)

        # Any shifted index < 0 is guaranteed outside the domain (solid).
        in_bounds = np.all(shifted >= 0, axis=1)       # (N,)

        # For out-of-bounds sites, use sentinel -1 which is never in fluid_encoded.
        enc_shifted = np.where(
            in_bounds,
            encode(shifted),
            np.int64(-1)
        )

        # Vectorised set membership check (one Python-level call per neighbor direction)
        is_missing = np.array(
            [e not in fluid_encoded for e in enc_shifted.tolist()],
            dtype=bool
        )

        is_wall_adjacent |= is_missing
        # Accumulate direction-into-solid vectors for wall normal estimation
        normal_sum[is_missing] += n_dir.astype(np.float64)

    wall_idx = np.nonzero(is_wall_adjacent)[0]
    logging.info(
        f"  Wall-adjacent nodes: {len(wall_idx)} "
        f"({100.0 * len(wall_idx) / N:.1f}% of fluid domain)"
    )

    # ── Compute WSS at each wall-adjacent node ────────────────────────────────
    dist_to_wall = 0.5 * voxel_size   # Mid-link bounce-back wall position

    node_ids   = []
    coords_out = []
    wss_out    = []

    for idx in wall_idx:
        n_solid = normal_sum[idx]                      # Points into solid
        n_norm  = np.linalg.norm(n_solid)

        if n_norm < 1e-12:
            # Degenerate: no net wall direction (should not occur for wall-adjacent nodes)
            logging.warning(f"  Degenerate wall normal at site index {idx}; skipping.")
            continue

        # Outward wall normal: from wall surface into fluid
        n_hat = -n_solid / n_norm

        v = velocity[idx]

        # Decompose velocity into wall-normal and tangential components
        v_normal_mag   = np.dot(v, n_hat)
        v_tangential   = v - v_normal_mag * n_hat
        v_tang_mag     = np.linalg.norm(v_tangential)

        # WSS = mu * |dv_tangential / dn| = mu * v_tang_mag / dist_to_wall
        wss = DYNAMIC_VISCOSITY * v_tang_mag / dist_to_wall

        # Physical coordinates (metres)
        xyz = origin + grid[idx].astype(np.float64) * voxel_size

        node_ids.append(int(idx))
        coords_out.append(xyz)
        wss_out.append(wss)

    coords_arr = np.array(coords_out)   # (M, 3)
    wss_arr    = np.array(wss_out)      # (M,)

    df = pd.DataFrame({
        'node_id'      : node_ids,
        'x'            : coords_arr[:, 0],
        'y'            : coords_arr[:, 1],
        'z'            : coords_arr[:, 2],
        'wss_magnitude': wss_arr,
    })

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)

    # ── Summary statistics ────────────────────────────────────────────────────
    n_low = int((wss_arr < LOW_WSS_THRESHOLD).sum())
    logging.info(
        f"  WSS (Pa): min={wss_arr.min():.4e}  "
        f"max={wss_arr.max():.4e}  "
        f"mean={wss_arr.mean():.4e}  "
        f"median={np.median(wss_arr):.4e}"
    )
    logging.info(
        f"  Nodes below {LOW_WSS_THRESHOLD} Pa (low-WSS zone): "
        f"{n_low} ({100.0 * n_low / len(wss_arr):.1f}%)"
    )
    logging.info(f"  Saved WSS output: {output_path}")

    return df


# ── Entry point: batch over all patients ─────────────────────────────────────

if __name__ == "__main__":
    BASE_DIR = Path(__file__).parent
    OUT_DIR  = BASE_DIR / "data/outputs_csv"
    WSS_DIR  = BASE_DIR / "data/wss"
    WSS_DIR.mkdir(parents=True, exist_ok=True)

    patient_dirs = sorted([d for d in OUT_DIR.iterdir() if d.is_dir()])

    if not patient_dirs:
        logging.warning(f"No patient directories found in {OUT_DIR}.")

    for patient_dir in patient_dirs:
        patient_id  = patient_dir.name
        csv_matches = list(patient_dir.glob(f"{patient_id}_fluid_data.csv"))

        if not csv_matches:
            logging.warning(
                f"No fluid data CSV for patient '{patient_id}' in {patient_dir}. Skipping."
            )
            continue

        csv_path    = csv_matches[0]
        output_path = WSS_DIR / f"{patient_id}_wss.csv"

        logging.info(f"======== WSS: {patient_id} ========")
        try:
            compute_wss(csv_path, output_path)
        except Exception:
            logging.error(f"Failed for {patient_id}:\n{traceback.format_exc()}")
            continue
