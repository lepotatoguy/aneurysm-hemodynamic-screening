#!/usr/bin/env python3
"""
animate_flow.py

Animates flow development across all 10 extracted timesteps for a patient case.
Shows velocity magnitude on a 2D cross-section slice through the vessel centroid,
plus a 3D scatter view coloured by velocity magnitude.

Saves: <patient_id>_flow_animation.gif

Usage
-----
    python animate_flow.py          # default: C0001
    python animate_flow.py C0005    # specific case
    python animate_flow.py C0001 C0005  # side-by-side comparison
"""

import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as mcolors
from pathlib import Path

BASE_DIR    = Path(__file__).parent
OUT_DIR     = BASE_DIR / "data/outputs_csv"
NUM_PATTERN = re.compile(r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?')

SUBSAMPLE   = 20     # Plot every Nth node (reduce for denser plot, increase for speed)
SLICE_AXIS  = 1      # Axis to slice through: 0=X, 1=Y, 2=Z
SLICE_TOL   = 3      # Tolerance in voxel units for slice selection


# ── Parser ────────────────────────────────────────────────────────────────────

def parse_all_timesteps(csv_path):
    """
    Parse all timestep blocks from a HemeLB fluid data CSV.

    Returns
    -------
    origin      : np.ndarray (3,) metres
    voxel_size  : float metres
    timesteps   : list of int
    grids       : list of np.ndarray (N,3) int32
    velocities  : list of np.ndarray (N,3) float64
    pressures   : list of np.ndarray (N,) float64
    """
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
            ts_num = int(NUM_PATTERN.findall(line)[0])
            ts_indices.append((i, ts_num))

    print(f"  Found {len(ts_indices)} timesteps: "
          f"{[ts for _, ts in ts_indices]}")

    timesteps  = []
    grids      = []
    velocities = []
    pressures  = []

    for block_idx, (line_idx, ts_num) in enumerate(ts_indices):
        # Data runs from line_idx+1 to next timestep marker (or EOF)
        if block_idx + 1 < len(ts_indices):
            end_idx = ts_indices[block_idx + 1][0]
        else:
            end_idx = len(lines)

        data_lines = [
            l.strip() for l in lines[line_idx + 1:end_idx]
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

        timesteps.append(ts_num)
        grids.append(np.array(grid_list, dtype=np.int32))
        velocities.append(np.array(vel_list, dtype=np.float64))
        pressures.append(np.array(pres_list, dtype=np.float64))

        print(f"  Timestep {ts_num}: {len(grid_list)} sites")

    return origin, voxel_size, timesteps, grids, velocities, pressures


# ── Animation ─────────────────────────────────────────────────────────────────

def animate_case(patient_id, ax_3d, ax_slice, timesteps, grids,
                 velocities, origin, voxel_size, vmin, vmax):
    """
    Build animation frames for one patient case on the given axes.
    Returns a list of artist lists for FuncAnimation.
    """
    frames = []

    for ts_idx, (ts, grid, vel) in enumerate(zip(timesteps, grids, velocities)):
        xyz    = origin + grid.astype(np.float64) * voxel_size
        mag    = np.linalg.norm(vel, axis=1)

        # Subsample for 3D scatter
        idx_sub = np.arange(0, len(xyz), SUBSAMPLE)
        xyz_sub = xyz[idx_sub]
        mag_sub = mag[idx_sub]

        # 2D cross-section slice
        centre_val  = np.median(xyz[:, SLICE_AXIS])
        tol_m       = SLICE_TOL * voxel_size
        slice_mask  = np.abs(xyz[:, SLICE_AXIS] - centre_val) < tol_m

        ax_3d.cla()
        ax_slice.cla()

        # 3D scatter
        sc3 = ax_3d.scatter(
            xyz_sub[:, 0], xyz_sub[:, 1], xyz_sub[:, 2],
            c=mag_sub, cmap='plasma', s=1.5, alpha=0.6,
            vmin=vmin, vmax=vmax
        )
        ax_3d.set_title(
            f"{patient_id} | Step {ts}\n"
            f"|v| mean={mag.mean():.4f} m/s  max={mag.max():.4f} m/s",
            color='white', fontsize=9
        )
        ax_3d.set_xlabel("x", fontsize=7, color='#888')
        ax_3d.set_ylabel("y", fontsize=7, color='#888')
        ax_3d.set_zlabel("z", fontsize=7, color='#888')
        ax_3d.tick_params(labelsize=6, colors='#888')
        ax_3d.set_facecolor('#0f0f0f')
        ax_3d.xaxis.pane.fill = False
        ax_3d.yaxis.pane.fill = False
        ax_3d.zaxis.pane.fill = False

        # 2D slice
        axis_labels = ['x', 'y', 'z']
        h_ax = [a for a in [0, 1, 2] if a != SLICE_AXIS]
        if slice_mask.sum() > 10:
            sc2 = ax_slice.scatter(
                xyz[slice_mask, h_ax[0]],
                xyz[slice_mask, h_ax[1]],
                c=mag[slice_mask], cmap='plasma', s=4, alpha=0.8,
                vmin=vmin, vmax=vmax
            )
            ax_slice.set_title(
                f"Cross-section (slice {axis_labels[SLICE_AXIS]}={centre_val*1000:.1f} mm)\n"
                f"Step {ts} | {slice_mask.sum()} nodes",
                color='white', fontsize=9
            )
        else:
            ax_slice.text(0.5, 0.5, 'No nodes in slice',
                         ha='center', va='center', color='white',
                         transform=ax_slice.transAxes)

        ax_slice.set_xlabel(f"{axis_labels[h_ax[0]]} (m)", fontsize=7, color='#888')
        ax_slice.set_ylabel(f"{axis_labels[h_ax[1]]} (m)", fontsize=7, color='#888')
        ax_slice.tick_params(labelsize=6, colors='#888')
        ax_slice.set_facecolor('#1a1a2e')

        frames.append(ts_idx)

    return frames


def run(case_ids):
    print(f"\nLoading data for: {case_ids}")

    # Load all data
    all_data = {}
    for pid in case_ids:
        csv_matches = list((OUT_DIR / pid).glob(f"{pid}_fluid_data.csv"))
        if not csv_matches:
            raise FileNotFoundError(
                f"No fluid data CSV for {pid} in {OUT_DIR / pid}"
            )
        print(f"\n{pid}:")
        origin, voxel_size, timesteps, grids, velocities, pressures = \
            parse_all_timesteps(csv_matches[0])
        all_data[pid] = (origin, voxel_size, timesteps, grids, velocities)

    # Shared velocity scale across all cases and timesteps
    all_mags = []
    for pid in case_ids:
        _, _, _, _, velocities = all_data[pid]
        for vel in velocities:
            all_mags.append(np.linalg.norm(vel, axis=1))
    all_mags_flat = np.concatenate(all_mags)
    vmin = 0.0
    vmax = np.percentile(all_mags_flat, 97)
    print(f"\nShared velocity colormap: 0 to {vmax:.4f} m/s (97th percentile)")

    n_cases = len(case_ids)
    fig = plt.figure(figsize=(8 * n_cases, 9))
    fig.patch.set_facecolor('#0f0f0f')

    n_ts = len(all_data[case_ids[0]][2])   # number of timesteps

    def update(frame_idx):
        fig.clear()
        fig.patch.set_facecolor('#0f0f0f')

        for col, pid in enumerate(case_ids):
            origin, voxel_size, timesteps, grids, velocities = all_data[pid]
            ts    = timesteps[frame_idx]
            grid  = grids[frame_idx]
            vel   = velocities[frame_idx]
            mag   = np.linalg.norm(vel, axis=1)
            xyz   = origin + grid.astype(np.float64) * voxel_size

            # 3D scatter
            ax3 = fig.add_subplot(2, n_cases, col + 1, projection='3d')
            ax3.set_facecolor('#0f0f0f')
            ax3.xaxis.pane.fill = False
            ax3.yaxis.pane.fill = False
            ax3.zaxis.pane.fill = False

            idx_sub = np.arange(0, len(xyz), SUBSAMPLE)
            ax3.scatter(
                xyz[idx_sub, 0], xyz[idx_sub, 1], xyz[idx_sub, 2],
                c=mag[idx_sub], cmap='plasma', s=1.5, alpha=0.6,
                vmin=vmin, vmax=vmax
            )
            ax3.set_title(
                f"{pid} | Timestep {ts} ({ts * 1e-5:.3f}s)\n"
                f"mean={mag.mean():.4f} m/s  max={mag.max():.4f} m/s",
                color='white', fontsize=9
            )
            ax3.tick_params(labelsize=6, colors='#555')

            # 2D slice
            ax2 = fig.add_subplot(2, n_cases, n_cases + col + 1)
            ax2.set_facecolor('#0f0f1a')

            centre_val = np.median(xyz[:, SLICE_AXIS])
            tol_m      = SLICE_TOL * voxel_size
            mask       = np.abs(xyz[:, SLICE_AXIS] - centre_val) < tol_m
            h_ax       = [a for a in [0, 1, 2] if a != SLICE_AXIS]

            if mask.sum() > 10:
                sc = ax2.scatter(
                    xyz[mask, h_ax[0]], xyz[mask, h_ax[1]],
                    c=mag[mask], cmap='plasma', s=5, alpha=0.85,
                    vmin=vmin, vmax=vmax
                )
                ax2.set_title(
                    f"Cross-section slice | {mask.sum()} nodes",
                    color='white', fontsize=8
                )
            else:
                ax2.text(0.5, 0.5, 'No nodes in slice',
                        ha='center', va='center', color='white',
                        transform=ax2.transAxes)

            ax2.tick_params(labelsize=6, colors='#555')
            for spine in ax2.spines.values():
                spine.set_edgecolor('#333')

        # Colorbar
        sm = plt.cm.ScalarMappable(
            cmap='plasma',
            norm=mcolors.Normalize(vmin=vmin, vmax=vmax)
        )
        sm.set_array([])
        cbar_ax = fig.add_axes([0.92, 0.15, 0.012, 0.7])
        cbar = fig.colorbar(sm, cax=cbar_ax)
        cbar.set_label('Velocity magnitude (m/s)', color='white', fontsize=9)
        cbar.ax.yaxis.set_tick_params(color='white', labelcolor='white')

        fig.suptitle(
            f"Flow Development: Transient to Steady State  "
            f"[{' vs '.join(case_ids)}]",
            color='white', fontsize=11, y=0.99
        )

    ani = animation.FuncAnimation(
        fig, update,
        frames=n_ts,
        interval=600,
        repeat=True
    )

    out_name = f"{'_vs_'.join(case_ids)}_flow_animation.gif"
    print(f"\nSaving animation to {out_name} ...")
    ani.save(out_name, writer='pillow', fps=1.5,
             savefig_kwargs={'facecolor': '#0f0f0f'})
    print(f"Saved: {out_name}")
    plt.show()


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    case_ids = sys.argv[1:] if len(sys.argv) > 1 else ["C0001"]
    run(case_ids)
