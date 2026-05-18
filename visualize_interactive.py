#!/usr/bin/env python3
"""
visualize_interactive.py

Interactive 3D visualization of WSS and flow data using PyVista.
Full zoom, rotate, pan, clip planes, and timestep scrubbing.

Install dependencies:
    pip install pyvista pyvistaqt

Usage
-----
    python visualize_interactive.py C0001       # WSS view
    python visualize_interactive.py C0001 flow  # Flow animation with slider
    python visualize_interactive.py C0001 C0005 # Side-by-side WSS comparison

Controls (PyVista window)
-------------------------
    Left-click drag   : rotate
    Right-click drag  : zoom
    Middle-click drag : pan
    Scroll wheel      : zoom in/out
    r                 : reset camera
    s                 : surface rendering
    w                 : wireframe rendering
    q                 : quit
"""

import sys
import re
import numpy as np
import pandas as pd
import trimesh
import pyvista as pv
from scipy.spatial import KDTree
from pathlib import Path

BASE_DIR    = Path(__file__).parent
RAW_DIR     = BASE_DIR / "data/raw_meshes"
WSS_DIR     = BASE_DIR / "data/wss"
OUT_DIR     = BASE_DIR / "data/outputs_csv"
MANIFEST    = BASE_DIR / "data/dataset_manifest.csv"

NUM_PATTERN = re.compile(r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?')

# Clinical thresholds
LOW_WSS     = 0.4    # Pa
HIGH_WSS    = 2.5    # Pa

# PyVista theme
pv.set_plot_theme("dark")


# ── Data loaders ──────────────────────────────────────────────────────────────

def load_manifest():
    if not MANIFEST.exists():
        return {}
    import csv
    info = {}
    with open(MANIFEST, newline='') as f:
        for row in csv.DictReader(f):
            pid = row.get('patient_id', '')
            info[pid] = {
                'status'  : row.get('ruptureStatus', ''),
                'location': row.get('aneurysmLocation', ''),
                'age'     : row.get('age', ''),
                'sex'     : row.get('sex', ''),
            }
    return info


def load_wss_csv(patient_id):
    path = WSS_DIR / f"{patient_id}_wss.csv"
    if not path.exists():
        raise FileNotFoundError(f"WSS file not found: {path}. Run compute_wss.py first.")
    df = pd.read_csv(path)
    return df[['x', 'y', 'z']].values, df['wss_magnitude'].values


def load_stl_as_pyvista(patient_id):
    """Load STL and return as PyVista PolyData with vertices in metres."""
    stl_path = RAW_DIR / f"{patient_id}.stl"
    if not stl_path.exists():
        raise FileNotFoundError(f"STL not found: {stl_path}")
    mesh_tm   = trimesh.load(stl_path, force='mesh')
    verts_m   = mesh_tm.vertices * 0.001   # mm -> m
    faces     = mesh_tm.faces

    # PyVista face format: [3, v0, v1, v2, 3, v0, v1, v2, ...]
    pv_faces  = np.hstack([np.full((len(faces), 1), 3), faces]).ravel()
    mesh_pv   = pv.PolyData(verts_m, pv_faces)
    return mesh_pv


def map_wss_to_mesh(mesh_pv, wss_xyz, wss_values):
    """KDTree map WSS from fluid nodes to mesh vertices."""
    tree       = KDTree(wss_xyz)
    _, nn_idx  = tree.query(mesh_pv.points)
    vertex_wss = wss_values[nn_idx]
    mesh_pv.point_data['WSS_Pa']      = vertex_wss
    mesh_pv.point_data['WSS_log']     = np.log10(np.clip(vertex_wss, 1e-6, None))
    mesh_pv.point_data['LowWSS_mask'] = (vertex_wss < LOW_WSS).astype(float)
    return mesh_pv


def parse_all_timesteps(patient_id):
    """Parse all timestep blocks from the fluid CSV."""
    csv_matches = list((OUT_DIR / patient_id).glob(f"{patient_id}_fluid_data.csv"))
    if not csv_matches:
        raise FileNotFoundError(f"No CSV for {patient_id}")

    with open(csv_matches[0], 'r') as f:
        lines = f.readlines()

    origin = None
    voxel_size = None
    ts_indices = []

    for i, line in enumerate(lines):
        if line.startswith('# Geometry origin'):
            nums   = NUM_PATTERN.findall(line)
            origin = np.array([float(n) for n in nums[:3]])
        elif line.startswith('# Voxel size'):
            voxel_size = float(NUM_PATTERN.findall(line)[0])
        elif line.startswith('# Timestep'):
            ts_num = int(NUM_PATTERN.findall(line)[0])
            ts_indices.append((i, ts_num))

    result = []
    for block_idx, (line_idx, ts_num) in enumerate(ts_indices):
        end_idx    = ts_indices[block_idx + 1][0] if block_idx + 1 < len(ts_indices) else len(lines)
        data_lines = [
            l.strip() for l in lines[line_idx + 1:end_idx]
            if l.strip() and not l.strip().startswith('#')
        ]
        grid_list = []
        vel_list  = []
        pres_list = []
        for line in data_lines:
            nums = [float(x) for x in NUM_PATTERN.findall(line)]
            if len(nums) == 7:
                grid_list.append(nums[0:3])
                vel_list.append(nums[3:6])
                pres_list.append(nums[6])

        grid = np.array(grid_list, dtype=np.int32)
        xyz  = origin + grid.astype(np.float64) * voxel_size
        vel  = np.array(vel_list, dtype=np.float64)
        mag  = np.linalg.norm(vel, axis=1)
        pres = np.array(pres_list, dtype=np.float64)

        result.append({
            'timestep' : ts_num,
            'time_s'   : ts_num * 1e-5,
            'xyz'      : xyz,
            'velocity' : vel,
            'speed'    : mag,
            'pressure' : pres,
        })

    print(f"  {patient_id}: {len(result)} timesteps, "
          f"{len(result[0]['xyz']):,} fluid sites")
    return result, origin, voxel_size


# ── View modes ────────────────────────────────────────────────────────────────

def view_wss(patient_ids):
    """
    Side-by-side WSS surface view with clip plane widget.
    Fully interactive: zoom, rotate, clip plane for cross-sections.
    """
    manifest = load_manifest()
    n        = len(patient_ids)
    pl       = pv.Plotter(shape=(1, n), window_size=[900 * n, 750])

    for col, pid in enumerate(patient_ids):
        pl.subplot(0, col)

        info   = manifest.get(pid, {})
        status = 'RUPTURED' if info.get('status') == 'R' else 'Unruptured'
        loc    = info.get('location', '')
        age    = info.get('age', '')
        sex    = info.get('sex', '')

        print(f"\nLoading {pid} ({status}, {loc})...")
        wss_xyz, wss_vals = load_wss_csv(pid)
        mesh_pv           = load_stl_as_pyvista(pid)
        mesh_pv           = map_wss_to_mesh(mesh_pv, wss_xyz, wss_vals)

        # Smooth the mesh slightly for better visual
        mesh_pv = mesh_pv.smooth(n_iter=20, relaxation_factor=0.1)

        vmax = np.percentile(wss_vals, 97)

        pl.add_mesh(
            mesh_pv,
            scalars='WSS_Pa',
            cmap='RdBu_r',
            clim=[0, vmax],
            smooth_shading=True,
            show_scalar_bar=False,
        )

        # Scalar bar
        pl.add_scalar_bar(
            title='WSS (Pa)',
            n_labels=5,
            fmt='%.2f',
            position_x=0.05,
            position_y=0.05,
            width=0.4,
            height=0.08,
            title_font_size=12,
            label_font_size=10,
            color='white',
        )

        # Add threshold contour for low WSS zone
        low_wss_region = mesh_pv.threshold(
            value=[0, LOW_WSS],
            scalars='WSS_Pa',
            invert=False
        )
        if low_wss_region.n_points > 0:
            pl.add_mesh(
                low_wss_region,
                color='#00bfff',
                opacity=0.3,
                label=f'Low WSS < {LOW_WSS} Pa',
                show_edges=False,
            )

        # Stats text
        n_low  = (wss_vals < LOW_WSS).sum()
        pct    = 100 * n_low / len(wss_vals)
        pl.add_text(
            f"{pid}  [{status} | {loc}]\n"
            f"Age: {age}  Sex: {sex}\n"
            f"mean WSS: {wss_vals.mean():.2f} Pa\n"
            f"max WSS: {wss_vals.max():.2f} Pa\n"
            f"Low-WSS nodes: {pct:.1f}%",
            position='upper_left',
            font_size=9,
            color='white',
        )

        pl.add_axes(color='white')

    pl.link_views()   # sync camera across subplots
    pl.add_title(
        "WSS Distribution — Use scroll to zoom, drag to rotate, "
        "right-click drag to pan",
        font_size=9, color='white'
    )
    print("\nInteractive WSS viewer ready.")
    print("Tip: Press 'c' to activate clip plane for cross-section view.")
    pl.show(cpos='iso')


def view_flow(patient_id, subsample=30):
    """
    Interactive flow viewer with timestep slider.
    Shows velocity magnitude on fluid nodes with a scrub bar.
    """
    manifest = load_manifest()
    info     = manifest.get(patient_id, {})
    status   = 'RUPTURED' if info.get('status') == 'R' else 'Unruptured'
    loc      = info.get('location', '')

    print(f"\nLoading all timesteps for {patient_id}...")
    ts_data, origin, voxel_size = parse_all_timesteps(patient_id)

    # Load WSS for final timestep overlay
    try:
        wss_xyz, wss_vals = load_wss_csv(patient_id)
        mesh_pv           = load_stl_as_pyvista(patient_id)
        mesh_pv           = map_wss_to_mesh(mesh_pv, wss_xyz, wss_vals)
        has_wss           = True
    except FileNotFoundError:
        has_wss = False
        mesh_pv = None

    # Shared velocity scale
    all_speeds = np.concatenate([d['speed'] for d in ts_data])
    vmax       = np.percentile(all_speeds, 97)
    vmin       = 0.0

    pl = pv.Plotter(window_size=[1100, 800])

    # Add vessel surface (semi-transparent WSS overlay)
    if has_wss:
        wss_vmax = np.percentile(wss_vals, 97)
        pl.add_mesh(
            mesh_pv,
            scalars='WSS_Pa',
            cmap='RdBu_r',
            clim=[0, wss_vmax],
            opacity=0.25,
            smooth_shading=True,
            show_scalar_bar=False,
            name='vessel_surface',
        )

    # Initial fluid point cloud
    ts0   = ts_data[-1]   # Start at final (converged) timestep
    idx_s = np.arange(0, len(ts0['xyz']), subsample)
    pc0   = pv.PolyData(ts0['xyz'][idx_s])
    pc0['speed'] = ts0['speed'][idx_s]
    pc0['vx']    = ts0['velocity'][idx_s, 0]
    pc0['vy']    = ts0['velocity'][idx_s, 1]
    pc0['vz']    = ts0['velocity'][idx_s, 2]

    pl.add_mesh(
        pc0,
        scalars='speed',
        cmap='plasma',
        clim=[vmin, vmax],
        point_size=3,
        render_points_as_spheres=True,
        show_scalar_bar=True,
        scalar_bar_args={
            'title'         : 'Speed (m/s)',
            'color'         : 'white',
            'title_font_size': 12,
            'label_font_size': 10,
        },
        name='fluid_nodes',
    )

    # ── Timestep slider ───────────────────────────────────────────────────────
    ts_labels = [f"Step {d['timestep']}  t={d['time_s']:.4f}s"
                 for d in ts_data]

    info_text = pl.add_text(
        f"{patient_id} [{status} | {loc}] | "
        f"Step {ts0['timestep']}  t={ts0['time_s']:.4f}s\n"
        f"mean speed: {ts0['speed'].mean():.4f} m/s  "
        f"max speed: {ts0['speed'].max():.4f} m/s",
        position='upper_left',
        font_size=9,
        color='white',
        name='info_text',
    )

    def update_timestep(value):
        # Find nearest timestep index
        idx   = int(round(value))
        idx   = max(0, min(idx, len(ts_data) - 1))
        ts    = ts_data[idx]

        idx_s = np.arange(0, len(ts['xyz']), subsample)
        pc    = pv.PolyData(ts['xyz'][idx_s])
        pc['speed'] = ts['speed'][idx_s]

        pl.add_mesh(
            pc,
            scalars='speed',
            cmap='plasma',
            clim=[vmin, vmax],
            point_size=3,
            render_points_as_spheres=True,
            show_scalar_bar=False,
            name='fluid_nodes',
        )
        pl.add_text(
            f"{patient_id} [{status} | {loc}] | "
            f"{ts_labels[idx]}\n"
            f"mean speed: {ts['speed'].mean():.4f} m/s  "
            f"max speed: {ts['speed'].max():.4f} m/s",
            position='upper_left',
            font_size=9,
            color='white',
            name='info_text',
        )

    pl.add_slider_widget(
        callback=update_timestep,
        rng=[0, len(ts_data) - 1],
        value=len(ts_data) - 1,
        title='Timestep',
        pointa=(0.1, 0.05),
        pointb=(0.9, 0.05),
        style='modern',
        color='white',
        title_color='white',
        title_height=0.02,
        fmt='%0.0f',
        slider_width=0.02,
        tube_width=0.005,
    )

    pl.add_axes(color='white')
    pl.add_title(
        f"Flow Development: {patient_id} — Drag slider to scrub timesteps",
        font_size=9, color='white'
    )

    print("\nInteractive flow viewer ready.")
    print("Drag the slider at the bottom to scrub through timesteps.")
    print("Semi-transparent surface = WSS distribution (red=high, blue=low)")
    pl.show(cpos='iso')


def view_comparison_flow(pid1, pid2, subsample=30):
    """Side-by-side flow comparison with shared timestep slider."""
    manifest = load_manifest()

    print(f"\nLoading timesteps for {pid1} and {pid2}...")
    ts1, _, _ = parse_all_timesteps(pid1)
    ts2, _, _ = parse_all_timesteps(pid2)

    all_speeds = np.concatenate(
        [d['speed'] for d in ts1] + [d['speed'] for d in ts2]
    )
    vmax = np.percentile(all_speeds, 97)

    # Load WSS meshes
    meshes = {}
    for pid in [pid1, pid2]:
        try:
            wss_xyz, wss_vals = load_wss_csv(pid)
            mesh              = load_stl_as_pyvista(pid)
            mesh              = map_wss_to_mesh(mesh, wss_xyz, wss_vals)
            meshes[pid]       = (mesh, wss_vals)
        except FileNotFoundError:
            meshes[pid] = (None, None)

    pl = pv.Plotter(shape=(1, 2), window_size=[1800, 800])

    for col, (pid, ts_data) in enumerate([(pid1, ts1), (pid2, ts2)]):
        pl.subplot(0, col)
        info   = manifest.get(pid, {})
        status = 'RUPTURED' if info.get('status') == 'R' else 'Unruptured'
        loc    = info.get('location', '')

        mesh, wss_vals = meshes[pid]
        if mesh is not None:
            wss_vmax = np.percentile(wss_vals, 97)
            pl.add_mesh(
                mesh,
                scalars='WSS_Pa',
                cmap='RdBu_r',
                clim=[0, wss_vmax],
                opacity=0.2,
                smooth_shading=True,
                show_scalar_bar=False,
                name=f'surface_{pid}',
            )

        ts    = ts_data[-1]
        idx_s = np.arange(0, len(ts['xyz']), subsample)
        pc    = pv.PolyData(ts['xyz'][idx_s])
        pc['speed'] = ts['speed'][idx_s]

        pl.add_mesh(
            pc,
            scalars='speed',
            cmap='plasma',
            clim=[0, vmax],
            point_size=3,
            render_points_as_spheres=True,
            show_scalar_bar=(col == 1),
            scalar_bar_args={'title': 'Speed (m/s)', 'color': 'white'},
            name=f'fluid_{pid}',
        )
        pl.add_text(
            f"{pid}  [{status} | {loc}]",
            position='upper_left',
            font_size=10,
            color='white',
        )
        pl.add_axes(color='white')

    pl.link_views()

    def update_both(value):
        idx = int(round(value))
        idx = max(0, min(idx, len(ts1) - 1))
        for col, (pid, ts_data) in enumerate([(pid1, ts1), (pid2, ts2)]):
            pl.subplot(0, col)
            ts    = ts_data[idx]
            idx_s = np.arange(0, len(ts['xyz']), subsample)
            pc    = pv.PolyData(ts['xyz'][idx_s])
            pc['speed'] = ts['speed'][idx_s]
            pl.add_mesh(
                pc,
                scalars='speed',
                cmap='plasma',
                clim=[0, vmax],
                point_size=3,
                render_points_as_spheres=True,
                show_scalar_bar=False,
                name=f'fluid_{pid}',
            )

    pl.subplot(0, 0)
    pl.add_slider_widget(
        callback=update_both,
        rng=[0, len(ts1) - 1],
        value=len(ts1) - 1,
        title='Timestep (both)',
        pointa=(0.1, 0.04),
        pointb=(0.9, 0.04),
        style='modern',
        color='white',
        title_color='white',
        fmt='%0.0f',
    )

    pl.add_title(
        f"Flow Comparison: {pid1} vs {pid2} — Drag slider to scrub",
        font_size=9, color='white'
    )
    pl.show(cpos='iso')


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith('-')]

    if len(args) == 0:
        print("Usage:")
        print("  python visualize_interactive.py C0001            # WSS view")
        print("  python visualize_interactive.py C0001 flow       # Flow slider")
        print("  python visualize_interactive.py C0001 C0005      # WSS comparison")
        print("  python visualize_interactive.py C0001 C0005 flow # Flow comparison")
        sys.exit(0)

    mode = 'wss'
    pids = []
    for a in args:
        if a.lower() == 'flow':
            mode = 'flow'
        else:
            pids.append(a)

    if mode == 'flow':
        if len(pids) == 1:
            view_flow(pids[0])
        elif len(pids) == 2:
            view_comparison_flow(pids[0], pids[1])
        else:
            print("Specify 1 or 2 patient IDs for flow mode.")
    else:
        if len(pids) >= 1:
            view_wss(pids[:2])
        else:
            print("Specify at least one patient ID.")
