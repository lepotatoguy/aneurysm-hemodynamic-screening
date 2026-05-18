#!/usr/bin/env python3
"""
visualize_interactive.py

Research-standard interactive 3D visualization of WSS and hemodynamic flow data.
Built on PyVista (VTK/OpenGL) for hardware-accelerated rendering.

Install
-------
    pip install pyvista

Usage
-----
    python visualize_interactive.py C0001              # WSS surface view
    python visualize_interactive.py C0001 C0005        # Side-by-side WSS comparison
    python visualize_interactive.py C0001 flow         # Flow with timestep slider
    python visualize_interactive.py C0001 C0005 flow   # Side-by-side flow comparison

Controls (all modes)
--------------------
    Scroll wheel        : zoom in/out
    Left-click drag     : rotate
    Middle-click drag   : pan
    Right-click drag    : zoom (fine)
    r                   : reset camera
    c                   : clip plane (interactive cross-section)
    p                   : save screenshot to screenshots/
    f                   : focus on picked point
    q / Escape          : quit

WSS mode buttons (top-left)
----------------------------
    WSS (Pa)    : absolute WSS magnitude
    nWSS        : normalised WSS (WSS / case mean) - pressure independent
    log(WSS)    : log-scale WSS - reveals low-magnitude spatial structure
    Low-WSS     : binary mask of low-WSS zone (< 0.4 Pa)
    Focus Peak  : snap camera to peak WSS impingement point

Research outputs per case
--------------------------
    - WSS magnitude surface coloured blue (low) to red (high)
    - Normalised WSS (nWSS = WSS/mean): velocity-independent, cross-case comparable
    - Low-WSS zone overlay: chronic low-WSS region associated with rupture risk
    - Impingement marker: yellow sphere at peak WSS location
    - Full statistics panel: mean, median, max, std, low/high WSS fractions, nWSS range
    - Clinical threshold annotations: 0.4 Pa (low) and 2.5 Pa (high)
    - Zoom slider on right edge
    - Camera orientation cube
"""

import sys
import re
import numpy as np
import pandas as pd
import trimesh
import pyvista as pv
from scipy.spatial import KDTree
from pathlib import Path

BASE_DIR       = Path(__file__).parent
RAW_DIR        = BASE_DIR / "data/raw_meshes"
WSS_DIR        = BASE_DIR / "data/wss"
OUT_DIR        = BASE_DIR / "data/outputs_csv"
MANIFEST       = BASE_DIR / "data/dataset_manifest.csv"
SCREENSHOT_DIR = BASE_DIR / "screenshots"
SCREENSHOT_DIR.mkdir(exist_ok=True)

NUM_PATTERN = re.compile(r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?')

LOW_WSS  = 0.4
HIGH_WSS = 2.5
WSS_CMAP  = 'RdBu_r'
FLOW_CMAP = 'plasma'

pv.set_plot_theme("dark")


# ── Manifest ──────────────────────────────────────────────────────────────────

def load_manifest():
    if not MANIFEST.exists():
        return {}
    import csv as _csv
    info = {}
    with open(MANIFEST, newline='') as f:
        for row in _csv.DictReader(f):
            pid = row.get('patient_id', '')
            info[pid] = {
                'status'  : row.get('ruptureStatus', ''),
                'location': row.get('aneurysmLocation', ''),
                'age'     : row.get('age', ''),
                'sex'     : row.get('sex', ''),
            }
    return info


def patient_label(pid, info):
    status = 'RUPTURED' if info.get('status') == 'R' else 'Unruptured'
    return status, info.get('location',''), info.get('age',''), info.get('sex','')


# ── Data loaders ──────────────────────────────────────────────────────────────

def load_wss_csv(patient_id):
    path = WSS_DIR / f"{patient_id}_wss.csv"
    if not path.exists():
        raise FileNotFoundError(f"WSS not found: {path}. Run compute_wss.py first.")
    df = pd.read_csv(path)
    return df[['x', 'y', 'z']].values, df['wss_magnitude'].values


def load_stl_as_pyvista(patient_id):
    stl_path = RAW_DIR / f"{patient_id}.stl"
    if not stl_path.exists():
        raise FileNotFoundError(f"STL not found: {stl_path}")
    tm       = trimesh.load(stl_path, force='mesh')
    verts_m  = tm.vertices * 0.001
    pv_faces = np.hstack([np.full((len(tm.faces), 1), 3), tm.faces]).ravel()
    return pv.PolyData(verts_m, pv_faces)


def map_wss_to_mesh(mesh_pv, wss_xyz, wss_vals):
    tree      = KDTree(wss_xyz)
    _, nn_idx = tree.query(mesh_pv.points)
    vwss      = wss_vals[nn_idx]
    nwss      = vwss / (vwss.mean() + 1e-12)
    mesh_pv.point_data['WSS_Pa']      = vwss
    mesh_pv.point_data['nWSS']        = nwss
    mesh_pv.point_data['WSS_log']     = np.log10(np.clip(vwss, 1e-6, None))
    mesh_pv.point_data['LowWSS_mask'] = (vwss < LOW_WSS).astype(float)
    return mesh_pv


def parse_all_timesteps(patient_id):
    csv_list = list((OUT_DIR / patient_id).glob(f"{patient_id}_fluid_data.csv"))
    if not csv_list:
        raise FileNotFoundError(f"No CSV for {patient_id}")
    with open(csv_list[0], 'r') as f:
        lines = f.readlines()
    origin, voxel_size, ts_indices = None, None, []
    for i, line in enumerate(lines):
        if line.startswith('# Geometry origin'):
            origin = np.array([float(n) for n in NUM_PATTERN.findall(line)[:3]])
        elif line.startswith('# Voxel size'):
            voxel_size = float(NUM_PATTERN.findall(line)[0])
        elif line.startswith('# Timestep'):
            ts_indices.append((i, int(NUM_PATTERN.findall(line)[0])))
    result = []
    for bi, (li, ts_num) in enumerate(ts_indices):
        end  = ts_indices[bi+1][0] if bi+1 < len(ts_indices) else len(lines)
        rows = [l.strip() for l in lines[li+1:end]
                if l.strip() and not l.strip().startswith('#')]
        gl, vl, pl_ = [], [], []
        for row in rows:
            nums = [float(x) for x in NUM_PATTERN.findall(row)]
            if len(nums) == 7:
                gl.append(nums[:3]); vl.append(nums[3:6]); pl_.append(nums[6])
        g   = np.array(gl, dtype=np.int32)
        xyz = origin + g.astype(np.float64) * voxel_size
        vel = np.array(vl, dtype=np.float64)
        result.append({'timestep': ts_num, 'time_s': ts_num*1e-5,
                       'xyz': xyz, 'velocity': vel,
                       'speed': np.linalg.norm(vel, axis=1),
                       'pressure': np.array(pl_, dtype=np.float64)})
    print(f"  {patient_id}: {len(result)} timesteps, {len(result[0]['xyz']):,} sites")
    return result, origin, voxel_size


# ── Shared UI components ──────────────────────────────────────────────────────

def add_zoom_slider(pl):
    zoom_state = {'factor': 1.0}
    def zoom_cb(value):
        ratio = float(value) / zoom_state['factor']
        pl.camera.zoom(ratio)
        zoom_state['factor'] = float(value)
    pl.add_slider_widget(
        callback=zoom_cb, rng=[0.2, 8.0], value=1.0,
        title='Zoom',
        pointa=(0.935, 0.18), pointb=(0.935, 0.82),
        style='modern', color='#aaaaaa', title_color='white',
        title_height=0.014, fmt='%.1fx',
        slider_width=0.014, tube_width=0.004,
    )


def add_screenshot_key(pl, prefix='screenshot'):
    counter = {'n': 0}
    def cb():
        counter['n'] += 1
        path = SCREENSHOT_DIR / f"{prefix}_{counter['n']:03d}.png"
        pl.screenshot(str(path))
        print(f"  Screenshot: {path}")
    pl.add_key_event('p', cb)


def add_stats_text(pl, pid, wss_vals, status, loc, age, sex):
    nwss      = wss_vals / (wss_vals.mean() + 1e-12)
    pct_low   = 100 * (wss_vals < LOW_WSS).mean()
    pct_high  = 100 * (wss_vals > HIGH_WSS).mean()
    pl.add_text(
        f"{pid}  |  {status}  |  {loc}\n"
        f"Age: {age}   Sex: {sex}\n"
        f"{'─'*26}\n"
        f"WSS mean   : {wss_vals.mean():.3f} Pa\n"
        f"WSS median : {np.median(wss_vals):.3f} Pa\n"
        f"WSS max    : {wss_vals.max():.2f} Pa\n"
        f"WSS std    : {wss_vals.std():.3f} Pa\n"
        f"{'─'*26}\n"
        f"Low  WSS (<{LOW_WSS}Pa): {pct_low:.1f}%\n"
        f"High WSS (>{HIGH_WSS}Pa): {pct_high:.1f}%\n"
        f"nWSS range: {nwss.min():.2f} – {nwss.max():.2f}",
        position='upper_left', font_size=8, color='white', font='courier',
    )


def add_impingement_marker(pl, mesh_pv, wss_vals):
    peak_idx   = np.argmax(mesh_pv.point_data['WSS_Pa'])
    peak_coord = mesh_pv.points[peak_idx]
    sphere     = pv.Sphere(radius=0.0008, center=peak_coord)
    pl.add_mesh(sphere, color='yellow', opacity=1.0)
    pl.add_point_labels(
        [peak_coord], [f"Peak WSS\n{wss_vals.max():.2f} Pa"],
        font_size=9, text_color='yellow',
        shape_color='black', shape_opacity=0.5, show_points=False,
    )
    return peak_coord


# ── WSS view ──────────────────────────────────────────────────────────────────

def view_wss(patient_ids):
    manifest = load_manifest()
    n        = len(patient_ids)
    pl       = pv.Plotter(shape=(1, n), window_size=[1000*n, 840])

    for col, pid in enumerate(patient_ids):
        pl.subplot(0, col)
        info                  = manifest.get(pid, {})
        status, loc, age, sex = patient_label(pid, info)

        print(f"\nLoading {pid} ({status}, {loc})...")
        wss_xyz, wss_vals = load_wss_csv(pid)
        mesh_pv           = load_stl_as_pyvista(pid)
        mesh_pv           = map_wss_to_mesh(mesh_pv, wss_xyz, wss_vals)
        mesh_pv           = mesh_pv.smooth(n_iter=30, relaxation_factor=0.1)

        vmax_pa   = np.percentile(wss_vals, 97)
        vmax_nwss = np.percentile(mesh_pv.point_data['nWSS'], 97)

        # Primary WSS surface
        pl.add_mesh(
            mesh_pv, scalars='WSS_Pa', cmap=WSS_CMAP,
            clim=[0.0, vmax_pa], smooth_shading=True,
            show_scalar_bar=True,
            scalar_bar_args={
                'title': 'WSS (Pa)', 'color': 'white',
                'title_font_size': 13, 'label_font_size': 11,
                'n_labels': 6, 'fmt': '%.2f',
                'position_x': 0.05, 'position_y': 0.06,
                'width': 0.32, 'height': 0.055,
            },
            name='wss_surface',
        )

        # Low-WSS overlay
        low_region = mesh_pv.threshold(value=[0.0, LOW_WSS], scalars='WSS_Pa')
        if low_region.n_points > 0:
            pl.add_mesh(low_region, color='#00cfff', opacity=0.4,
                        show_edges=False, name='low_wss_overlay')

        # Impingement marker
        peak_coord = add_impingement_marker(pl, mesh_pv, wss_vals)

        # Stats
        add_stats_text(pl, pid, wss_vals, status, loc, age, sex)

        # Orientation cube + zoom slider + screenshot
        pl.add_camera_orientation_widget()
        add_zoom_slider(pl)
        add_screenshot_key(pl, prefix=f"wss_{pid}")
        pl.add_axes(color='white', xlabel='X', ylabel='Y', zlabel='Z')

        # ── Toggle buttons ────────────────────────────────────────────────────
        win_h = pl.window_size[1]
        y0    = 0.80
        bh    = 0.058
        bsize = 24

        def _make_scalar_cb(scalar, vmin, vmax, _pl=pl, _mesh=mesh_pv):
            def cb(_state):
                _pl.update_scalars(scalar, mesh=_mesh, render=False)
                _pl.update_scalar_bar_range([vmin, vmax])
                _pl.render()
            return cb

        log_min = float(mesh_pv.point_data['WSS_log'].min())
        log_max = float(mesh_pv.point_data['WSS_log'].max())

        buttons = [
            ('WSS (Pa)',  '#e74c3c', _make_scalar_cb('WSS_Pa',      0,       vmax_pa)),
            ('nWSS',      '#3498db', _make_scalar_cb('nWSS',         0,       vmax_nwss)),
            ('log(WSS)',  '#2ecc71', _make_scalar_cb('WSS_log',      log_min, log_max)),
            ('Low-WSS',   '#00cfff', _make_scalar_cb('LowWSS_mask',  0,       1)),
        ]

        for i, (label, colour, cb) in enumerate(buttons):
            ypos = int(win_h * (y0 - i * bh))
            pl.add_checkbox_button_widget(
                cb, value=(i == 0),
                position=(10, ypos), size=bsize,
                color_on=colour, color_off='#444444',
            )
            pl.add_text(label, position=(42, ypos + 4),
                        font_size=8, color='white')

        # Focus aneurysm button
        def focus_cb(_state, _pl=pl, _coord=peak_coord, _mesh=mesh_pv):
            b = _mesh.bounds
            span = max(b[1]-b[0], b[3]-b[2], b[5]-b[4])
            _pl.camera.focal_point = _coord.tolist()
            _pl.camera.position    = (_coord + np.array([0, 0, span*0.35])).tolist()
            _pl.camera.zoom(3.0)
            _pl.render()

        ypos_focus = int(win_h * (y0 - len(buttons) * bh))
        pl.add_checkbox_button_widget(
            focus_cb, value=False,
            position=(10, ypos_focus), size=bsize,
            color_on='#f39c12', color_off='#444444',
        )
        pl.add_text('Focus Peak', position=(42, ypos_focus + 4),
                    font_size=8, color='white')

    pl.link_views()
    pl.add_title(
        "WSS  |  Scroll=zoom  Drag=rotate  'c'=clip  'p'=screenshot  "
        "Buttons: scalar toggle  Right slider: zoom",
        font_size=8, color='#999999',
    )
    print("\nWSS viewer ready.")
    print("  Top-left buttons: WSS(Pa) | nWSS | log(WSS) | Low-WSS | Focus Peak")
    print("  Right edge: zoom slider")
    print("  'c' key: clip plane for cross-sections")
    print("  'p' key: save screenshot")
    pl.show(cpos='iso')


# ── Flow view ─────────────────────────────────────────────────────────────────

def view_flow(patient_id, subsample=25):
    manifest              = load_manifest()
    info                  = manifest.get(patient_id, {})
    status, loc, age, sex = patient_label(patient_id, info)

    print(f"\nLoading timesteps for {patient_id}...")
    ts_data, origin, voxel_size = parse_all_timesteps(patient_id)

    try:
        wss_xyz, wss_vals = load_wss_csv(patient_id)
        mesh_pv           = load_stl_as_pyvista(patient_id)
        mesh_pv           = map_wss_to_mesh(mesh_pv, wss_xyz, wss_vals)
        has_wss           = True
    except FileNotFoundError:
        has_wss = False
        mesh_pv = None

    all_spd  = np.concatenate([d['speed']    for d in ts_data])
    all_pres = np.concatenate([d['pressure'] for d in ts_data])
    vmax_spd  = np.percentile(all_spd, 97)
    vmin_pres = np.percentile(all_pres, 3)
    vmax_pres = np.percentile(all_pres, 97)

    pl = pv.Plotter(window_size=[1250, 860])

    if has_wss:
        wss_vmax = np.percentile(wss_vals, 97)
        pl.add_mesh(mesh_pv, scalars='WSS_Pa', cmap=WSS_CMAP,
                    clim=[0, wss_vmax], opacity=0.18,
                    smooth_shading=True, show_scalar_bar=False,
                    name='vessel_surface')

    ts0   = ts_data[-1]
    idx_s = np.arange(0, len(ts0['xyz']), subsample)
    pc0   = pv.PolyData(ts0['xyz'][idx_s])
    pc0['speed']    = ts0['speed'][idx_s]
    pc0['pressure'] = ts0['pressure'][idx_s]
    pc0['vx']       = ts0['velocity'][idx_s, 0]
    pc0['vy']       = ts0['velocity'][idx_s, 1]
    pc0['vz']       = ts0['velocity'][idx_s, 2]

    pl.add_mesh(pc0, scalars='speed', cmap=FLOW_CMAP,
                clim=[0, vmax_spd], point_size=3,
                render_points_as_spheres=True, show_scalar_bar=True,
                scalar_bar_args={'title': 'Speed (m/s)', 'color': 'white',
                                 'title_font_size': 12, 'label_font_size': 10,
                                 'position_x': 0.05, 'position_y': 0.06,
                                 'width': 0.32, 'height': 0.055},
                name='fluid_nodes')

    def info_str(ts):
        return (f"{patient_id}  [{status} | {loc}]  Age:{age}  Sex:{sex}\n"
                f"Step {ts['timestep']}  t={ts['time_s']:.4f}s\n"
                f"Speed  mean={ts['speed'].mean():.4f}  max={ts['speed'].max():.4f} m/s\n"
                f"Press  mean={ts['pressure'].mean():.4e} Pa")

    pl.add_text(info_str(ts0), position='upper_left',
                font_size=8, color='white', font='courier', name='info')

    show_state   = {'scalar': 'speed', 'arrows': False}
    arrows_actor = {'obj': None}

    def update_ts(value):
        idx   = max(0, min(int(round(value)), len(ts_data)-1))
        ts    = ts_data[idx]
        idx_s = np.arange(0, len(ts['xyz']), subsample)
        pc    = pv.PolyData(ts['xyz'][idx_s])
        pc['speed']    = ts['speed'][idx_s]
        pc['pressure'] = ts['pressure'][idx_s]
        sc    = show_state['scalar']
        clim  = [0, vmax_spd] if sc == 'speed' else [vmin_pres, vmax_pres]
        cmap  = FLOW_CMAP if sc == 'speed' else 'coolwarm'
        pl.add_mesh(pc, scalars=sc, cmap=cmap, clim=clim,
                    point_size=3, render_points_as_spheres=True,
                    show_scalar_bar=False, name='fluid_nodes')
        pl.add_text(info_str(ts), position='upper_left',
                    font_size=8, color='white', font='courier', name='info')
        if show_state['arrows']:
            _refresh_arrows(ts, idx_s)

    pl.add_slider_widget(
        callback=update_ts,
        rng=[0, len(ts_data)-1], value=len(ts_data)-1,
        title='Timestep',
        pointa=(0.10, 0.04), pointb=(0.88, 0.04),
        style='modern', color='white', title_color='white',
        title_height=0.018, fmt='%0.0f',
        slider_width=0.018, tube_width=0.004,
    )

    win_h = pl.window_size[1]
    y0    = 0.80
    bh    = 0.058
    bsize = 24

    def toggle_pressure(state):
        show_state['scalar'] = 'pressure' if state else 'speed'
        update_ts(len(ts_data)-1)

    pl.add_checkbox_button_widget(toggle_pressure, value=False,
        position=(10, int(win_h*y0)), size=bsize,
        color_on='#9b59b6', color_off='#444444')
    pl.add_text('Pressure', position=(42, int(win_h*y0)+4),
                font_size=8, color='white')

    def _refresh_arrows(ts, idx_s):
        if arrows_actor['obj'] is not None:
            pl.remove_actor(arrows_actor['obj'])
        pc_a = pv.PolyData(ts['xyz'][idx_s])
        pc_a['velocity'] = ts['velocity'][idx_s]
        glyphs = pc_a.glyph(orient='velocity', scale='velocity',
                             factor=0.002, geom=pv.Arrow(), tolerance=0.05)
        arrows_actor['obj'] = pl.add_mesh(
            glyphs, color='#ffff55', opacity=0.75, name='vel_arrows')

    def toggle_arrows(state):
        show_state['arrows'] = state
        if not state and arrows_actor['obj'] is not None:
            pl.remove_actor(arrows_actor['obj'])
            arrows_actor['obj'] = None
        elif state:
            ts    = ts_data[-1]
            idx_s = np.arange(0, len(ts['xyz']), subsample)
            _refresh_arrows(ts, idx_s)

    pl.add_checkbox_button_widget(toggle_arrows, value=False,
        position=(10, int(win_h*(y0-bh))), size=bsize,
        color_on='#f1c40f', color_off='#444444')
    pl.add_text('Vel arrows', position=(42, int(win_h*(y0-bh))+4),
                font_size=8, color='white')

    pl.add_camera_orientation_widget()
    add_zoom_slider(pl)
    add_screenshot_key(pl, prefix=f"flow_{patient_id}")
    pl.add_axes(color='white')
    pl.add_title(
        f"Flow: {patient_id}  |  Slider=timestep  Buttons=Pressure/Arrows  "
        f"Right=zoom  'p'=screenshot  'c'=clip",
        font_size=8, color='#999999')

    print("\nFlow viewer ready.")
    print("  Bottom slider: scrub through timesteps")
    print("  Buttons: Pressure toggle | Velocity arrows")
    print("  Right slider: zoom")
    pl.show(cpos='iso')


# ── Flow comparison ───────────────────────────────────────────────────────────

def view_comparison_flow(pid1, pid2, subsample=25):
    manifest = load_manifest()
    print(f"\nLoading {pid1} and {pid2}...")
    ts1, _, _ = parse_all_timesteps(pid1)
    ts2, _, _ = parse_all_timesteps(pid2)

    all_spd = np.concatenate([d['speed'] for d in ts1+ts2])
    vmax    = np.percentile(all_spd, 97)

    meshes = {}
    for pid in [pid1, pid2]:
        try:
            wss_xyz, wss_vals = load_wss_csv(pid)
            mesh              = load_stl_as_pyvista(pid)
            mesh              = map_wss_to_mesh(mesh, wss_xyz, wss_vals)
            meshes[pid]       = (mesh, wss_vals)
        except FileNotFoundError:
            meshes[pid]       = (None, None)

    pl = pv.Plotter(shape=(1, 2), window_size=[1960, 860])

    for col, (pid, ts_data) in enumerate([(pid1, ts1), (pid2, ts2)]):
        pl.subplot(0, col)
        info                  = manifest.get(pid, {})
        status, loc, age, sex = patient_label(pid, info)
        mesh, wss_vals        = meshes[pid]

        if mesh is not None:
            wss_vmax = np.percentile(wss_vals, 97)
            pl.add_mesh(mesh, scalars='WSS_Pa', cmap=WSS_CMAP,
                        clim=[0, wss_vmax], opacity=0.18,
                        smooth_shading=True, show_scalar_bar=False,
                        name=f'surface_{pid}')

        ts    = ts_data[-1]
        idx_s = np.arange(0, len(ts['xyz']), subsample)
        pc    = pv.PolyData(ts['xyz'][idx_s])
        pc['speed'] = ts['speed'][idx_s]

        pl.add_mesh(pc, scalars='speed', cmap=FLOW_CMAP, clim=[0, vmax],
                    point_size=3, render_points_as_spheres=True,
                    show_scalar_bar=(col == 1),
                    scalar_bar_args={'title': 'Speed (m/s)', 'color': 'white'},
                    name=f'fluid_{pid}')

        pl.add_text(f"{pid}  [{status} | {loc}]\nAge:{age}  Sex:{sex}",
                    position='upper_left', font_size=9, color='white')
        pl.add_camera_orientation_widget()
        add_zoom_slider(pl)
        pl.add_axes(color='white')

    pl.link_views()
    add_screenshot_key(pl, prefix=f"cmp_{pid1}_{pid2}")

    def update_both(value):
        idx = max(0, min(int(round(value)), len(ts1)-1))
        for col, (pid, ts_data) in enumerate([(pid1, ts1), (pid2, ts2)]):
            pl.subplot(0, col)
            ts    = ts_data[idx]
            idx_s = np.arange(0, len(ts['xyz']), subsample)
            pc    = pv.PolyData(ts['xyz'][idx_s])
            pc['speed'] = ts['speed'][idx_s]
            pl.add_mesh(pc, scalars='speed', cmap=FLOW_CMAP, clim=[0, vmax],
                        point_size=3, render_points_as_spheres=True,
                        show_scalar_bar=False, name=f'fluid_{pid}')

    pl.subplot(0, 0)
    pl.add_slider_widget(update_both,
        rng=[0, len(ts1)-1], value=len(ts1)-1,
        title='Timestep (synced)',
        pointa=(0.10, 0.04), pointb=(0.88, 0.04),
        style='modern', color='white', title_color='white', fmt='%0.0f')

    pl.add_title(f"Flow: {pid1} vs {pid2}  |  Slider=timestep  'p'=screenshot",
                 font_size=8, color='#999999')
    pl.show(cpos='iso')


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith('-')]

    if len(args) == 0:
        print(__doc__)
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
        elif len(pids) >= 2:
            view_comparison_flow(pids[0], pids[1])
        else:
            print("Specify 1 or 2 patient IDs.")
    else:
        if len(pids) >= 1:
            view_wss(pids[:2])
        else:
            print("Specify at least one patient ID.")
