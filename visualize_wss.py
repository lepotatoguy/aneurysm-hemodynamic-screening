#!/usr/bin/env python3
"""
visualize_wss.py

Side-by-side WSS visualization comparing ruptured vs unruptured aneurysm cases.
Colors the STL surface mesh by WSS magnitude mapped from the nearest fluid lattice node.

Usage
-----
    python visualize_wss.py                    # default: C0001 vs C0005
    python visualize_wss.py C0001 C0005        # explicit case IDs
    python visualize_wss.py C0001              # single case

Clinical reference thresholds
------------------------------
    Low WSS  : < 0.4 Pa  (associated with aneurysm growth and rupture risk)
    High WSS : > 2.5 Pa  (associated with wall degradation at neck)
"""

import sys
import numpy as np
import pandas as pd
import trimesh
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import KDTree
from pathlib import Path

# ── Configuration ─────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).parent
RAW_DIR  = BASE_DIR / "data/raw_meshes"
WSS_DIR  = BASE_DIR / "data/wss"

LOW_WSS_THRESHOLD  = 0.4    # Pa - clinical low-WSS threshold
HIGH_WSS_THRESHOLD = 2.5    # Pa - clinical high-WSS threshold

# Clinical colormap: deep blue (low WSS) -> white -> deep red (high WSS)
# Matches the convention used in published hemodynamics papers
WSS_CMAP = LinearSegmentedColormap.from_list(
    "wss_clinical",
    ["#053061", "#2166ac", "#4393c3", "#92c5de", "#f7f7f7",
     "#f4a582", "#d6604d", "#b2182b", "#67001f"],
    N=512
)


# ── WSS mapping ───────────────────────────────────────────────────────────────

def load_wss(patient_id):
    """Load WSS CSV. Returns (xyz array in metres, wss array in Pa)."""
    path = WSS_DIR / f"{patient_id}_wss.csv"
    if not path.exists():
        raise FileNotFoundError(f"WSS file not found: {path}. Run compute_wss.py first.")
    df = pd.read_csv(path)
    xyz = df[['x', 'y', 'z']].values
    wss = df['wss_magnitude'].values
    return xyz, wss


def map_wss_to_stl(stl_path, wss_xyz, wss_values):
    """
    Map WSS values from fluid lattice nodes onto STL surface vertices
    using nearest-neighbour KDTree query.

    STL vertices are in mm; WSS coordinates are in metres.
    Converts STL to metres before matching.

    Returns
    -------
    mesh        : trimesh.Trimesh
    vertex_wss  : np.ndarray (V,)  WSS at each STL vertex in Pa
    face_wss    : np.ndarray (F,)  WSS at each face (mean of 3 vertices) in Pa
    """
    mesh = trimesh.load(stl_path, force='mesh')
    vertices_m = mesh.vertices * 0.001   # mm -> m

    tree = KDTree(wss_xyz)
    _, nn_idx = tree.query(vertices_m)

    vertex_wss = wss_values[nn_idx]
    face_wss   = vertex_wss[mesh.faces].mean(axis=1)

    return mesh, vertex_wss, face_wss


# ── Plotting ──────────────────────────────────────────────────────────────────

def plot_case(ax, patient_id, mesh, face_wss, manifest_label, vmin, vmax):
    """
    Render one case on the given 3D axes.
    Face colour encodes WSS magnitude. Low-WSS faces are highlighted.
    """
    vertices  = mesh.vertices * 0.001   # mm -> m
    triangles = vertices[mesh.faces]    # (F, 3, 3)

    # Normalise face WSS to [0,1] for colormap using shared scale
    norm      = mcolors.Normalize(vmin=vmin, vmax=vmax)
    face_rgba = WSS_CMAP(norm(face_wss))

    # Highlight low-WSS faces with extra transparency to make them pop
    low_wss_mask = face_wss < LOW_WSS_THRESHOLD
    face_rgba[low_wss_mask, 3]  = 1.0    # fully opaque for low-WSS
    face_rgba[~low_wss_mask, 3] = 0.85

    poly = Poly3DCollection(
        triangles,
        facecolor=face_rgba,
        edgecolor='none',
        linewidth=0,
        shade=False
    )
    ax.add_collection3d(poly)

    # Set axis limits from mesh bounds
    mins = vertices.min(axis=0)
    maxs = vertices.max(axis=0)
    centre = (mins + maxs) / 2
    half_range = (maxs - mins).max() / 2 * 1.1

    ax.set_xlim(centre[0] - half_range, centre[0] + half_range)
    ax.set_ylim(centre[1] - half_range, centre[1] + half_range)
    ax.set_zlim(centre[2] - half_range, centre[2] + half_range)

    # Labels
    n_low  = low_wss_mask.sum()
    pct    = 100 * n_low / len(face_wss)
    mean_w = face_wss.mean()
    max_w  = face_wss.max()

    ax.set_title(
        f"{patient_id}  [{manifest_label}]\n"
        f"mean={mean_w:.2f} Pa  max={max_w:.1f} Pa\n"
        f"low-WSS faces (<{LOW_WSS_THRESHOLD} Pa): {pct:.1f}%",
        fontsize=10, pad=8
    )
    ax.set_xlabel("x (m)", fontsize=7)
    ax.set_ylabel("y (m)", fontsize=7)
    ax.set_zlabel("z (m)", fontsize=7)
    ax.tick_params(labelsize=6)
    ax.view_init(elev=25, azim=-60)

    return norm


def visualize(case_ids, labels=None):
    """
    Main visualization function. case_ids is a list of patient ID strings.
    labels is an optional list of rupture status strings.
    """
    n = len(case_ids)
    if labels is None:
        labels = [''] * n

    # Load WSS and build mesh data for all cases
    all_face_wss = []
    meshes       = []
    all_wss_vals = []

    for pid in case_ids:
        stl_path = RAW_DIR / f"{pid}.stl"
        if not stl_path.exists():
            raise FileNotFoundError(f"STL not found: {stl_path}")

        wss_xyz, wss_vals = load_wss(pid)
        mesh, _, face_wss = map_wss_to_stl(stl_path, wss_xyz, wss_vals)

        meshes.append(mesh)
        all_face_wss.append(face_wss)
        all_wss_vals.append(wss_vals)
        print(f"{pid} [{labels[case_ids.index(pid)]}]: "
              f"mean={wss_vals.mean():.3f} Pa  "
              f"max={wss_vals.max():.2f} Pa  "
              f"low-WSS={100*(wss_vals < LOW_WSS_THRESHOLD).mean():.1f}%")

    # Shared colormap scale across all cases for direct comparison
    vmin = 0.0
    vmax = np.percentile(np.concatenate(all_wss_vals), 97)   # 97th percentile cap
    print(f"\nColormap range: 0 to {vmax:.2f} Pa (97th percentile)")

    # Figure layout
    fig = plt.figure(figsize=(7 * n, 8))
    fig.patch.set_facecolor('#0f0f0f')

    axes = []
    for i in range(n):
        ax = fig.add_subplot(1, n, i + 1, projection='3d')
        ax.set_facecolor('#0f0f0f')
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.xaxis.pane.set_edgecolor('#333333')
        ax.yaxis.pane.set_edgecolor('#333333')
        ax.zaxis.pane.set_edgecolor('#333333')
        ax.tick_params(colors='#888888')
        ax.xaxis.label.set_color('#888888')
        ax.yaxis.label.set_color('#888888')
        ax.zaxis.label.set_color('#888888')
        axes.append(ax)

    for i, (pid, label, mesh, face_wss) in enumerate(
            zip(case_ids, labels, meshes, all_face_wss)):
        norm = plot_case(axes[i], pid, mesh, face_wss, label, vmin, vmax)
        axes[i].title.set_color('white')

    # Shared colorbar
    sm  = cm.ScalarMappable(cmap=WSS_CMAP, norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    cbar    = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label('WSS (Pa)', color='white', fontsize=11)
    cbar.ax.yaxis.set_tick_params(color='white', labelcolor='white')

    # Threshold lines on colorbar
    norm_obj = mcolors.Normalize(vmin=vmin, vmax=vmax)
    for threshold, label_text, color in [
        (LOW_WSS_THRESHOLD,  f'Low WSS\n{LOW_WSS_THRESHOLD} Pa',  '#00bfff'),
        (HIGH_WSS_THRESHOLD, f'High WSS\n{HIGH_WSS_THRESHOLD} Pa', '#ff6b35'),
    ]:
        if vmin <= threshold <= vmax:
            y_pos = norm_obj(threshold)
            cbar_ax.axhline(y=y_pos, color=color, linewidth=1.5, linestyle='--')
            cbar_ax.text(1.3, y_pos, label_text, color=color,
                        fontsize=7, va='center', transform=cbar_ax.transAxes)

    fig.suptitle(
        'Wall Shear Stress Distribution — Ruptured vs Unruptured Cerebral Aneurysm\n'
        'Blue = low WSS (rupture risk zone)  |  Red = high WSS (neck jet impingement)',
        color='white', fontsize=11, y=0.98
    )

    plt.tight_layout(rect=[0, 0, 0.91, 0.96])
    plt.savefig('wss_comparison.png', dpi=150, bbox_inches='tight',
                facecolor='#0f0f0f')
    print("\nSaved: wss_comparison.png")
    plt.show()


# ── CLI ───────────────────────────────────────────────────────────────────────

def load_rupture_status(case_ids):
    """Try to read rupture status from dataset_manifest.csv."""
    manifest_path = BASE_DIR / "data/dataset_manifest.csv"
    if not manifest_path.exists():
        return [''] * len(case_ids)

    import csv
    status_map = {}
    with open(manifest_path, newline='') as f:
        for row in csv.DictReader(f):
            pid = row.get('patient_id', '')
            rs  = row.get('ruptureStatus', '')
            loc = row.get('aneurysmLocation', '')
            if rs == 'R':
                status_map[pid] = f'RUPTURED | {loc}'
            elif rs == 'U':
                status_map[pid] = f'Unruptured | {loc}'
            else:
                status_map[pid] = rs

    return [status_map.get(pid, '') for pid in case_ids]


if __name__ == "__main__":
    if len(sys.argv) >= 2:
        case_ids = sys.argv[1:]
    else:
        case_ids = ["C0001", "C0005"]

    labels = load_rupture_status(case_ids)
    visualize(case_ids, labels)
