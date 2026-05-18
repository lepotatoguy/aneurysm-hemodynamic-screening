#!/usr/bin/env python3
"""
verify_simulation.py

Prints diagnostic statistics from a HemeLB fluid data CSV to verify
simulation physics. Run this on any patient CSV to check:

  1. Header metadata (origin, voxel size, timesteps present)
  2. Velocity statistics at the final timestep
  3. Inlet vs outlet vs interior velocity comparison
  4. Pressure gradient direction check
  5. Reynolds number estimate
  6. WSS stats from the corresponding wss CSV (if available)

Usage
-----
    python verify_simulation.py              # runs on C0001 by default
    python verify_simulation.py C0002        # runs on specified patient
"""

import re
import sys
import numpy as np
from pathlib import Path

# ── Config ────────────────────────────────────────────────────────────────────
BASE_DIR   = Path(__file__).parent
PATIENT_ID = sys.argv[1] if len(sys.argv) > 1 else "C0001"
CSV_PATH   = BASE_DIR / f"data/outputs_csv/{PATIENT_ID}/{PATIENT_ID}_fluid_data.csv"
WSS_PATH   = BASE_DIR / f"data/wss/{PATIENT_ID}_wss.csv"
PR2_PATH   = BASE_DIR / f"data/raw_meshes/{PATIENT_ID}.pr2"

# Physical constants
MU          = 0.0035      # Pa.s, dynamic viscosity
RHO         = 1060.0      # kg/m3, blood density
NUM_PATTERN = re.compile(r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?')

DIVIDER = "=" * 60

# ── Helpers ───────────────────────────────────────────────────────────────────

def parse_csv(csv_path):
    with open(csv_path, 'r') as f:
        lines = f.readlines()

    origin     = None
    voxel_size = None
    ts_indices = []

    for i, line in enumerate(lines):
        if line.startswith('# Geometry origin'):
            nums   = NUM_PATTERN.findall(line)
            origin = np.array([float(n) for n in nums[:3]])
        elif line.startswith('# Voxel size'):
            nums       = NUM_PATTERN.findall(line)
            voxel_size = float(nums[0])
        elif line.startswith('# Timestep'):
            ts_indices.append((i, int(NUM_PATTERN.findall(line)[0])))

    # Parse final timestep block
    last_idx = ts_indices[-1][0]
    data_lines = [
        l.strip() for l in lines[last_idx + 1:]
        if l.strip() and not l.strip().startswith('#')
    ]

    grid_list = []
    vel_list  = []
    pres_list = []

    for line in data_lines:
        nums = [float(x) for x in NUM_PATTERN.findall(line)]
        if len(nums) != 7:
            continue
        grid_list.append(nums[0:3])
        vel_list.append(nums[3:6])
        pres_list.append(nums[6])

    return (origin, voxel_size,
            [t for _, t in ts_indices],
            np.array(grid_list, dtype=np.int32),
            np.array(vel_list,  dtype=np.float64),
            np.array(pres_list, dtype=np.float64))


def parse_pr2_iolets(pr2_path):
    """Parse iolet centres, normals, radii and types from .pr2 file."""
    TOP_LEVEL = {'OutputGeometryFile','OutputXmlFile','SeedPoint',
                 'StlFile','StlFileUnitId','TimeStepSeconds','VoxelSize','DurationSeconds'}
    iolets = []
    current = None
    current_sub = None
    in_iolets = False

    with open(pr2_path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s == 'Iolets:':
                in_iolets = True
                continue
            if not in_iolets:
                continue
            if s.split(':')[0].strip() in TOP_LEVEL:
                if current:
                    iolets.append(current)
                    current = None
                in_iolets = False
                continue
            if s == '- Centre:':
                if current:
                    iolets.append(current)
                current = {'centre':[0.,0.,0.],'normal':[0.,0.,0.],'radius':0.,'type':None}
                current_sub = 'centre'
                continue
            if s == 'Normal:':
                current_sub = 'normal'
                continue
            if s in ('Pressure:', 'Centre:') or s.startswith('Name:'):
                current_sub = None
                continue
            if s.startswith('Radius:') and current:
                current['radius'] = float(s.split(':',1)[1].strip())
                current_sub = None
                continue
            if s.startswith('Type:') and current:
                current['type'] = s.split(':',1)[1].strip()
                current_sub = None
                continue
            if current_sub in ('centre','normal') and current:
                if s.startswith('x:'):
                    current[current_sub][0] = float(s[2:].strip())
                elif s.startswith('y:'):
                    current[current_sub][1] = float(s[2:].strip())
                elif s.startswith('z:'):
                    current[current_sub][2] = float(s[2:].strip())

    if current:
        iolets.append(current)
    return iolets


def get_iolet_nodes(xyz, voxel_size, iolets):
    """Return boolean masks for nodes near each iolet plane."""
    masks = []
    for iolet in iolets:
        centre_m = np.array(iolet['centre']) / 1000.0
        radius_m = iolet['radius']            / 1000.0
        n_hat    = np.array(iolet['normal'], dtype=np.float64)
        n_hat   /= np.linalg.norm(n_hat)

        disp          = xyz - centre_m
        plane_dist    = np.abs(disp @ n_hat)
        in_plane_dist = np.linalg.norm(disp - (disp @ n_hat)[:,None] * n_hat, axis=1)

        mask = (plane_dist < 3.0 * voxel_size) & (in_plane_dist < radius_m)
        masks.append(mask)
    return masks


# ── Main verification ─────────────────────────────────────────────────────────

def main():
    print(f"\n{DIVIDER}")
    print(f"  SIMULATION VERIFICATION: {PATIENT_ID}")
    print(DIVIDER)

    if not CSV_PATH.exists():
        print(f"  ERROR: CSV not found at {CSV_PATH}")
        return

    # ── 1. Parse CSV ──────────────────────────────────────────────────────────
    origin, voxel_size, timesteps, grid, velocity, pressure = parse_csv(CSV_PATH)
    xyz      = origin + grid.astype(np.float64) * voxel_size
    vel_mag  = np.linalg.norm(velocity, axis=1)
    N        = len(grid)

    print(f"\n  [1] HEADER METADATA")
    print(f"      Origin (m)       : {origin}")
    print(f"      Voxel size (m)   : {voxel_size}")
    print(f"      Timesteps in CSV : {timesteps}")
    print(f"      Final timestep   : {timesteps[-1]}")
    print(f"      Fluid sites      : {N:,}")
    print(f"      Physical time    : {timesteps[-1] * voxel_size / (voxel_size / (voxel_size * 5000)):.4f} s"
          if False else f"      Physical time    : {timesteps[-1] * 2e-4:.4f} s  (at dt=2e-4 s)")

    # ── 2. Velocity statistics ────────────────────────────────────────────────
    print(f"\n  [2] VELOCITY STATISTICS (final timestep, all fluid nodes)")
    print(f"      Min magnitude    : {vel_mag.min():.4e} m/s")
    print(f"      Max magnitude    : {vel_mag.max():.4e} m/s")
    print(f"      Mean magnitude   : {vel_mag.mean():.4e} m/s")
    print(f"      Median magnitude : {np.median(vel_mag):.4e} m/s")
    print(f"      Std deviation    : {vel_mag.std():.4e} m/s")

    # Top 10 highest velocity nodes
    top10_idx = np.argsort(vel_mag)[-10:][::-1]
    print(f"\n      Top 10 highest velocity nodes:")
    print(f"      {'Index':>8}  {'Vel mag (m/s)':>15}  {'x (m)':>10}  {'y (m)':>10}  {'z (m)':>10}")
    for i in top10_idx:
        print(f"      {i:>8}  {vel_mag[i]:>15.4e}  {xyz[i,0]:>10.5f}  {xyz[i,1]:>10.5f}  {xyz[i,2]:>10.5f}")

    # ── 3. Pressure statistics ────────────────────────────────────────────────
    print(f"\n  [3] PRESSURE STATISTICS (Pa)")
    print(f"      Min pressure     : {pressure.min():.4e} Pa")
    print(f"      Max pressure     : {pressure.max():.4e} Pa")
    print(f"      Mean pressure    : {pressure.mean():.4e} Pa")
    print(f"      Pressure range   : {pressure.max() - pressure.min():.4e} Pa")
    print(f"      0.1 mmHg in Pa   : {0.1 * 133.322:.4e} Pa  (expected range)")

    # ── 4. Iolet analysis ─────────────────────────────────────────────────────
    print(f"\n  [4] IOLET VELOCITY ANALYSIS")
    if PR2_PATH.exists():
        iolets = parse_pr2_iolets(PR2_PATH)
        masks  = get_iolet_nodes(xyz, voxel_size, iolets)

        for iolet, mask in zip(iolets, masks):
            n_nodes = mask.sum()
            if n_nodes == 0:
                print(f"      {iolet['type']:6s}: NO NODES DETECTED (radius={iolet['radius']:.2f} mm)")
                continue

            n_hat    = np.array(iolet['normal']) / np.linalg.norm(iolet['normal'])
            v_nodes  = velocity[mask]
            vm_nodes = vel_mag[mask]

            # Normal flux
            flux = float(np.sum(v_nodes @ n_hat)) * voxel_size**2

            print(f"\n      {iolet['type']:6s} | nodes={n_nodes} | "
                  f"centre radius={iolet['radius']:.2f} mm")
            print(f"        Vel mean    : {vm_nodes.mean():.4e} m/s")
            print(f"        Vel max     : {vm_nodes.max():.4e} m/s")
            print(f"        Flux        : {flux:.4e} m^3/s  "
                  f"({'INTO' if flux > 0 else 'OUT OF'} domain)")
    else:
        print(f"      PR2 file not found at {PR2_PATH}. Skipping iolet analysis.")

    # ── 5. Reynolds number estimate ───────────────────────────────────────────
    print(f"\n  [5] REYNOLDS NUMBER ESTIMATE")
    # Use max velocity and median vessel radius from iolets if available
    v_char = vel_mag.max()
    # Estimate characteristic length from domain extent
    domain_extent = (xyz.max(axis=0) - xyz.min(axis=0))
    L_char = domain_extent.min() / 2.0   # rough vessel radius estimate
    Re = RHO * v_char * L_char / MU
    print(f"      Max velocity     : {v_char:.4e} m/s")
    print(f"      Domain extent    : {domain_extent * 1000} mm")
    print(f"      Char. length est : {L_char*1000:.2f} mm  (half of smallest domain dim)")
    print(f"      Reynolds number  : {Re:.2f}")
    print(f"      Expected range   : 100-400 for cerebral vessels at physiol. pressure")
    print(f"      At 0.1 mmHg      : Re ~{Re:.1f}  (scales linearly with pressure)")

    # ── 6. Convergence spot check ─────────────────────────────────────────────
    print(f"\n  [6] SPATIAL VELOCITY DISTRIBUTION CHECK")
    # Divide domain into thirds along x-axis and compare mean velocity
    x_coords = xyz[:, 0]
    x_min, x_max = x_coords.min(), x_coords.max()
    x_third = (x_max - x_min) / 3.0
    r1 = vel_mag[x_coords < x_min + x_third]
    r2 = vel_mag[(x_coords >= x_min + x_third) & (x_coords < x_min + 2*x_third)]
    r3 = vel_mag[x_coords >= x_min + 2*x_third]
    print(f"      Mean vel, X-third 1 (inlet side) : {r1.mean():.4e} m/s  ({len(r1):,} nodes)")
    print(f"      Mean vel, X-third 2 (middle)     : {r2.mean():.4e} m/s  ({len(r2):,} nodes)")
    print(f"      Mean vel, X-third 3 (outlet side): {r3.mean():.4e} m/s  ({len(r3):,} nodes)")

    # ── 7. WSS stats ──────────────────────────────────────────────────────────
    print(f"\n  [7] WSS STATISTICS (from {WSS_PATH.name})")
    if WSS_PATH.exists():
        import csv
        wss_vals = []
        with open(WSS_PATH) as f:
            reader = csv.DictReader(f)
            for row in reader:
                wss_vals.append(float(row['wss_magnitude']))
        wss = np.array(wss_vals)
        print(f"      Wall-adjacent nodes  : {len(wss):,}")
        print(f"      WSS min (Pa)         : {wss.min():.4e}")
        print(f"      WSS max (Pa)         : {wss.max():.4e}")
        print(f"      WSS mean (Pa)        : {wss.mean():.4e}")
        print(f"      WSS median (Pa)      : {np.median(wss):.4e}")
        print(f"      WSS std (Pa)         : {wss.std():.4e}")
        print(f"      Nodes below 0.4 Pa   : {(wss < 0.4).sum():,}  ({100*(wss<0.4).mean():.1f}%)")
        print(f"      Nodes above 1.5 Pa   : {(wss > 1.5).sum():,}  ({100*(wss>1.5).mean():.1f}%)")
        nwss = wss / wss.mean()
        print(f"      Normalised WSS range : {nwss.min():.3f} to {nwss.max():.3f}")
        print(f"      (nWSS = WSS / mean_WSS, pressure-independent)")
    else:
        print(f"      WSS file not found. Run compute_wss.py first.")

    print(f"\n{DIVIDER}\n")


if __name__ == "__main__":
    main()
