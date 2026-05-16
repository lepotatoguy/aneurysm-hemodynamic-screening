# import os
# import re
# import numpy as np
# import matplotlib.pyplot as plt
# from pathlib import Path

# # ==================================================
# # CONFIGURATION
# # ==================================================
# CSV_PATH = Path("data/outputs_csv/model/model_fluid_data.csv")
# SUBSAMPLE_RATE = 50  # Plot every N-th point for Matplotlib to stay responsive

# def load_and_parse_hlb(filepath):
#     """Parses custom HemeLB output format with bracketed strings and comments."""
#     print(f"Reading HemeLB raw data matrix from: {filepath}")
    
#     grids = []
#     velocities = []
#     pressures = []
#     timesteps = []
#     current_ts = 0
    
#     with open(filepath, 'r') as f:
#         for line in f:
#             line = line.strip()
#             if not line:
#                 continue
#             if line.startswith('#'):
#                 if "Timestep" in line:
#                     current_ts = int(line.split()[-1])
#                 continue
            
#             # Row template: [ 29  31 182], [ 1.015e-12 -1.176e-11  2.288e-12], -2.502e-11
#             try:
#                 parts = line.split("],")
#                 if len(parts) == 3:
#                     # Parse grid positions [x y z]
#                     g_str = parts[0].replace("[", "").strip()
#                     grid = [int(x) for x in g_str.split()]
                    
#                     # Parse velocity vectors [vx vy vz]
#                     v_str = parts[1].replace("[", "").strip()
#                     vel = [float(x) for x in v_str.split()]
                    
#                     # Parse pressure scalar
#                     p = float(parts[2].strip())
                    
#                     grids.append(grid)
#                     velocities.append(vel)
#                     pressures.append(p)
#                     timesteps.append(current_ts)
#             except Exception:
#                 continue
                
#     print(f"Successfully loaded {len(grids)} total lattice records.")
#     return np.array(grids), np.array(velocities), np.array(pressures), np.array(timesteps)

# def plot_with_matplotlib(X, V, P, mag):
#     """Generates a fast static 3D spatial check plot using Matplotlib."""
#     print(f"Generating Matplotlib 3D visualization (subsampled 1/{SUBSAMPLE_RATE})...")
    
#     # Subsample data array to keep matplotlib interactively smooth
#     X_sub = X[::SUBSAMPLE_RATE]
#     mag_sub = mag[::SUBSAMPLE_RATE]
    
#     fig = plt.figure(figsize=(10, 8))
#     ax = fig.add_subplot(111, projection='3d')
    
#     # Scatter grid nodes colored by local velocity scalar magnitude
#     sc = ax.scatter(X_sub[:, 0], X_sub[:, 1], X_sub[:, 2], 
#                     c=mag_sub, cmap='jet', s=2, alpha=0.6)
    
#     ax.set_title("Vascular Fluid Grid Node Distribution (Colored by Velocity Magnitude)")
#     ax.set_xlabel("Lattice X")
#     ax.set_ylabel("Lattice Y")
#     ax.set_zlabel("Lattice Z")
#     fig.colorbar(sc, ax=ax, label='Velocity Magnitude (Lattice Units)')
#     plt.show()

# def plot_with_pyvista(X, V, P, mag):
#     """Generates a high-fidelity interactive 3D render window using PyVista."""
#     try:
#         import pyvista as pv
#     except ImportError:
#         print("PyVista is not installed. Skipping high-fidelity rendering.")
#         return

#     print("Generating hardware-accelerated interactive 3D PyVista workspace...")
#     # Instantiate unstructured point cloud object geometry
#     point_cloud = pv.PolyData(X)
#     point_cloud["Velocity Magnitude"] = mag
#     point_cloud["Pressure"] = P
#     point_cloud["Vectors"] = V

#     # Setup rendering environment layout
#     plotter = pv.Plotter(window_size=[1000, 800])
#     plotter.background_color = 'white'
    
#     # Add scalar mesh surface nodes
#     plotter.add_mesh(point_cloud, scalars="Velocity Magnitude", cmap="jet", 
#                      point_size=3.0, render_points_as_spheres=True, opacity=0.8)
    
#     # Generate vector arrows tracking local velocity direction paths
#     arrows = point_cloud.glyph(orient="Vectors", scale="Velocity Magnitude", factor=5.0e10)
#     plotter.add_mesh(arrows, color="black", opacity=0.4)
    
#     plotter.add_title("Steady-State Vascular Hemodynamics Profile", color="black", font_size=12)
#     plotter.add_axes()
#     plotter.show()

# if __name__ == "__main__":
#     if not CSV_PATH.exists():
#         print(f"Error: Target data matrix not found at: {CSV_PATH}")
#         exit(1)
        
#     grids, vels, pressures, timesteps = load_and_parse_hlb(CSV_PATH)
    
#     # Isolate only the final simulated iteration frame (converged state)
#     final_ts = np.max(timesteps)
#     print(f"Filtering dataset parameters for the final simulated step: frame {final_ts}")
#     idx = (timesteps == final_ts)
    
#     X = grids[idx]
#     V = vels[idx]
#     P = pressures[idx]
    
#     # Calculate scalar velocity vector norms across nodes
#     vel_magnitude = np.linalg.norm(V, axis=1)
    
#     # Run visualizations
#     plot_with_matplotlib(X, V, P, vel_magnitude)
#     plot_with_pyvista(X, V, P, vel_magnitude)



#######################################################

# import os
# import xml.etree.ElementTree as ET
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.colors import to_rgba
# from pathlib import Path
# import trimesh
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# # ==================================================
# # CONFIGURATION - TUNE THESE FOR YOUR VISUAL PREFERENCE
# # ==================================================
# CSV_PATH = Path("data/outputs_csv/model/model_fluid_data.csv")
# XML_PATH = Path("data/processed_gmy/model_input.xml")
# STL_PATH = Path("data/raw_meshes/model.stl")

# SUBSAMPLE_RATE = 50       # Higher = faster responsiveness; Lower = denser fluid points
# MESH_ALPHA = 0.07         # Face transparency of the vessel walls (kept low so inside nodes pop)
# MESH_COLOR = "gainsboro"  # Color of the vessel body shell

# # --- NEW BORDER CONTROLS ---
# MESH_EDGE_COLOR = "black" # Color of the vessel contour borders
# MESH_EDGE_ALPHA = 0.15    # Transparency of the border lines (0.0 to 1.0)
# MESH_LINEWIDTH = 0.2      # Thickness of the border lines

# def load_and_parse_hlb(filepath):
#     """Parses custom HemeLB output format with bracketed strings and comments."""
#     print(f"Reading HemeLB raw data matrix from: {filepath}")
    
#     grids = []
#     velocities = []
#     pressures = []
#     timesteps = []
#     current_ts = 0
    
#     with open(filepath, 'r') as f:
#         for line in f:
#             line = line.strip()
#             if not line:
#                 continue
#             if line.startswith('#'):
#                 if "Timestep" in line:
#                     current_ts = int(line.split()[-1])
#                 continue
            
#             try:
#                 parts = line.split("],")
#                 if len(parts) == 3:
#                     # Parse grid positions [x y z]
#                     g_str = parts[0].replace("[", "").strip()
#                     grid = [int(x) for x in g_str.split()]
                    
#                     # Parse velocity vectors [vx vy vz]
#                     v_str = parts[1].replace("[", "").strip()
#                     vel = [float(x) for x in v_str.split()]
                    
#                     # Parse pressure scalar
#                     p = float(parts[2].strip())
                    
#                     grids.append(grid)
#                     velocities.append(vel)
#                     pressures.append(p)
#                     timesteps.append(current_ts)
#             except Exception:
#                 continue
                
#     print(f"Successfully loaded {len(grids)} total lattice records.")
#     return np.array(grids), np.array(velocities), np.array(pressures), np.array(timesteps)

# def parse_xml_metadata(xml_filepath):
#     """Reads simulation metadata parameters (origin, voxel size, iolets)."""
#     if not xml_filepath.exists():
#         print(f"Warning: XML file not found at {xml_filepath}. Missing outer hull alignment markers.")
#         return None, None, []
        
#     print(f"Extracting spatial scale parameters from: {xml_filepath.name}")
#     tree = ET.parse(xml_filepath)
#     root = tree.getroot()
    
#     # Parse Origin: "(val, val, val)" -> numpy array
#     origin_str = root.find(".//simulation/origin").get("value").strip("()")
#     origin = np.array([float(x) for x in origin_str.split(",")])
    
#     # Parse Voxel size
#     voxel_size = float(root.find(".//simulation/voxel_size").get("value"))
    
#     iolet_markers = []
#     # Process Inlets
#     for idx, inlet in enumerate(root.findall(".//inlets/inlet")):
#         pos_str = inlet.find("position").get("value").strip("()")
#         pos_meters = np.array([float(x) for x in pos_str.split(",")])
#         pos_lattice = (pos_meters - origin) / voxel_size
#         iolet_markers.append({"name": f"Inlet {idx+1}", "coord": pos_lattice, "type": "inlet"})
        
#     # Process Outlets
#     for idx, outlet in enumerate(root.findall(".//outlets/outlet")):
#         pos_str = outlet.find("position").get("value").strip("()")
#         pos_meters = np.array([float(x) for x in pos_str.split(",")])
#         pos_lattice = (pos_meters - origin) / voxel_size
#         iolet_markers.append({"name": f"Outlet {idx+1}", "coord": pos_lattice, "type": "outlet"})
        
#     return origin, voxel_size, iolet_markers

# def plot_with_matplotlib(X, V, P, mag, origin, voxel_size, iolet_markers):
#     """Generates an advanced 3D visualization combining translucent STL walls, fluid nodes, and iolets."""
#     print(f"Generating Matplotlib 3D visualization (subsampled 1/{SUBSAMPLE_RATE})...")
    
#     X_sub = X[::SUBSAMPLE_RATE]
#     mag_sub = mag[::SUBSAMPLE_RATE]
    
#     fig = plt.figure(figsize=(12, 9))
#     ax = fig.add_subplot(111, projection='3d')
    
#     # 1. LOAD AND OVERLAY TRANSLUCENT STL MESH HULL WITH DEFINED BORDERS
#     if STL_PATH.exists() and origin is not None and voxel_size is not None:
#         print(f"Loading outer STL vascular boundary shell from: {STL_PATH.name}")
#         mesh = trimesh.load(STL_PATH)
        
#         # Convert vertices from mm to meters
#         vertices_meters = mesh.vertices * 0.001
        
#         # Project 3D coordinates perfectly into identical Lattice unit space
#         vertices_lattice = (vertices_meters - origin) / voxel_size
#         triangles_lattice = vertices_lattice[mesh.faces]
        
#         # Generate an independent soft RGBA color mapping for the borders
#         edge_rgba = to_rgba(MESH_EDGE_COLOR, alpha=MESH_EDGE_ALPHA)
        
#         # Add the mesh collection with explicit edge styling
#         mesh_col = Poly3DCollection(triangles_lattice, alpha=MESH_ALPHA, 
#                                     facecolor=MESH_COLOR, edgecolor=edge_rgba, 
#                                     linewidth=MESH_LINEWIDTH, zorder=1)
#         ax.add_collection3d(mesh_col)
    
#     # 2. PLOT THE COLOR-CODED INTERNAL FLUID SITES
#     sc = ax.scatter(X_sub[:, 0], X_sub[:, 1], X_sub[:, 2], 
#                     c=mag_sub, cmap='jet', s=2.5, alpha=0.5, label='Internal Fluid Domain', zorder=2)
    
#     # 3. INJECT THE INLET & OUTLET POSITION MARKERS
#     for iolet in iolet_markers:
#         color = 'red' if iolet['type'] == 'inlet' else 'blue'
#         marker_shape = '^' if iolet['type'] == 'inlet' else 'v'
        
#         # Draw high-visibility spatial coordinate markers
#         ax.scatter(iolet['coord'][0], iolet['coord'][1], iolet['coord'][2],
#                    color=color, marker=marker_shape, s=300, edgecolor='black', 
#                    linewidth=2.0, label=iolet['name'], zorder=10)
        
#         # Offset text string label for clarity
#         ax.text(iolet['coord'][0] + 6, iolet['coord'][1] + 6, iolet['coord'][2] + 6, 
#                 iolet['name'], color='black', fontsize=11, weight='bold', zorder=11)
    
#     # Title & Label settings
#     ax.set_title("Converged Fluid Profile with Defined Outer Vascular Borders", fontsize=14, weight='bold')
#     ax.set_xlabel("Lattice X")
#     ax.set_ylabel("Lattice Y")
#     ax.set_zlabel("Lattice Z")
    
#     # Filter legend handles to avoid duplicate entries
#     handles, labels = ax.get_legend_handles_labels()
#     by_label = dict(zip(labels, handles))
#     ax.legend(by_label.values(), by_label.keys(), loc='upper left', markerscale=1.4)
    
#     fig.colorbar(sc, ax=ax, shrink=0.6, label='Velocity Magnitude (Lattice Units)')
#     plt.show()

# if __name__ == "__main__":
#     if not CSV_PATH.exists():
#         print(f"Error: Target data matrix not found at: {CSV_PATH}")
#         exit(1)
        
#     grids, vels, pressures, timesteps = load_and_parse_hlb(CSV_PATH)
#     origin, voxel_size, iolet_markers = parse_xml_metadata(XML_PATH)
    
#     # Filter for final steady-state frame
#     final_ts = np.max(timesteps)
#     print(f"Filtering dataset parameters for the final simulated step: frame {final_ts}")
#     idx = (timesteps == final_ts)
    
#     X = grids[idx]
#     V = vels[idx]
#     P = pressures[idx]
    
#     vel_magnitude = np.linalg.norm(V, axis=1)
    
#     # Launch visualization
#     plot_with_matplotlib(X, V, P, vel_magnitude, origin, voxel_size, iolet_markers)


import os
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from pathlib import Path
import trimesh
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ==================================================================
# CONFIGURATION - TUNE THESE FOR YOUR SURROGATE TRAINING ANALYSIS
# ==================================================================
CSV_PATH = Path("data/outputs_csv/model/model_fluid_data.csv")
XML_PATH = Path("data/processed_gmy/model_input.xml")
STL_PATH = Path("data/raw_meshes/model.stl")

SUBSAMPLE_RATE = 40       # Denser sampling (40) makes fluid patterns easier to trace
ENABLE_CUTAWAY = True     # Slices the vessel in half lengthwise to reveal the interior

# Cutaway plane orientation: 'Y' cuts horizontally, 'X' cuts vertically, 'Z' cuts cross-section
CLIPPING_AXIS = 'Y'       

# Back-Shell Visual Styling (Forms the background cradle)
MESH_COLOR = "lightgray"  # Opaque backing color to catch light/shadow depth
MESH_ALPHA = 0.0         # Transparency of the back vessel wall
MESH_EDGE_COLOR = "black" # Defines the internal geometric ridges and contours
MESH_EDGE_ALPHA = 0.35    # High contrast outline to make vessel borders razor sharp
MESH_LINEWIDTH = 0.3      # Fineness of the mesh triangle lines

def load_and_parse_hlb(filepath):
    """Parses custom HemeLB output format with bracketed strings and comments."""
    print(f"Reading HemeLB raw data matrix from: {filepath}")
    
    grids, velocities, pressures, timesteps = [], [], [], []
    current_ts = 0
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                if "Timestep" in line:
                    current_ts = int(line.split()[-1])
                continue
            
            try:
                parts = line.split("],")
                if len(parts) == 3:
                    grid = [int(x) for x in parts[0].replace("[", "").split()]
                    vel = [float(x) for x in parts[1].replace("[", "").split()]
                    p = float(parts[2].strip())
                    
                    grids.append(grid)
                    velocities.append(vel)
                    pressures.append(p)
                    timesteps.append(current_ts)
            except Exception:
                continue
                
    print(f"Successfully loaded {len(grids)} total lattice records.")
    return np.array(grids), np.array(velocities), np.array(pressures), np.array(timesteps)

def parse_xml_metadata(xml_filepath):
    """Reads simulation metadata parameters (origin, voxel size, iolets)."""
    if not xml_filepath.exists():
        return None, None, []
        
    tree = ET.parse(xml_filepath)
    root = tree.getroot()
    
    origin_str = root.find(".//simulation/origin").get("value").strip("()")
    origin = np.array([float(x) for x in origin_str.split(",")])
    voxel_size = float(root.find(".//simulation/voxel_size").get("value"))
    
    iolet_markers = []
    for idx, inlet in enumerate(root.findall(".//inlets/inlet")):
        pos_str = inlet.find("position").get("value").strip("()")
        pos_lattice = (np.array([float(x) for x in pos_str.split(",")]) - origin) / voxel_size
        iolet_markers.append({"name": f"Inlet {idx+1}", "coord": pos_lattice, "type": "inlet"})
        
    for idx, outlet in enumerate(root.findall(".//outlets/outlet")):
        pos_str = outlet.find("position").get("value").strip("()")
        pos_lattice = (np.array([float(x) for x in pos_str.split(",")]) - origin) / voxel_size
        iolet_markers.append({"name": f"Outlet {idx+1}", "coord": pos_lattice, "type": "outlet"})
        
    return origin, voxel_size, iolet_markers

def plot_with_matplotlib(X, V, P, mag, origin, voxel_size, iolet_markers):
    """Generates an advanced 3D cutaway visualization separating outer walls from inner fluid fields."""
    print("Generating Cutaway Longitudinal 3D Plot Layout...")
    
    # Filter final converged timestep
    X_sub = X[::SUBSAMPLE_RATE]
    mag_sub = mag[::SUBSAMPLE_RATE]
    
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # 1. COMPUTE AND RENDER CLIPPED OPAQUE HALF-SHELL BOUNDARY
    if STL_PATH.exists() and origin is not None and voxel_size is not None:
        mesh = trimesh.load(STL_PATH)
        vertices_lattice = (mesh.vertices * 0.001 - origin) / voxel_size
        triangles_lattice = vertices_lattice[mesh.faces]
        
        if ENABLE_CUTAWAY:
            # Calculate the mid-point center of the mesh to establish our cutting plane
            face_centroids = triangles_lattice.mean(axis=1)
            mesh_center = vertices_lattice.mean(axis=0)
            
            # Select only triangles that reside on the back-half of the selected axis
            axis_idx = {'X': 0, 'Y': 1, 'Z': 2}[CLIPPING_AXIS]
            cut_mask = face_centroids[:, axis_idx] < mesh_center[axis_idx]
            triangles_to_render = triangles_lattice[cut_mask]
            title_suffix = f" (Longitudinal {CLIPPING_AXIS}-Cutaway View)"
        else:
            triangles_to_render = triangles_lattice
            title_suffix = " (Closed Translucent Shell View)"
            
        edge_rgba = to_rgba(MESH_EDGE_COLOR, alpha=MESH_EDGE_ALPHA)
        mesh_col = Poly3DCollection(triangles_to_render, alpha=MESH_ALPHA, 
                                    facecolor=MESH_COLOR, edgecolor=edge_rgba, 
                                    linewidth=MESH_LINEWIDTH, zorder=1)
        ax.add_collection3d(mesh_col)
        print(f"Rendered {len(triangles_to_render)} background contour panels safely.")

    # 2. PLOT THE COLOR-CODED FLUID PATTERNS INSIDE THE CHANNEL
    sc = ax.scatter(X_sub[:, 0], X_sub[:, 1], X_sub[:, 2], 
                    c=mag_sub, cmap='jet', s=3.5, alpha=0.7, 
                    label='Internal Lattice Fluid Nodes', zorder=2)
    
    # 3. INJECT HIGHLIGHTED BOUNDARY PORTS (INLETS/OUTLETS)
    for iolet in iolet_markers:
        color = 'red' if iolet['type'] == 'inlet' else 'blue'
        marker_shape = '^' if iolet['type'] == 'inlet' else 'v'
        
        ax.scatter(iolet['coord'][0], iolet['coord'][1], iolet['coord'][2],
                   color=color, marker=marker_shape, s=350, edgecolor='black', 
                   linewidth=2.5, label=iolet['name'], zorder=10)
        
        ax.text(iolet['coord'][0] + 5, iolet['coord'][1] + 5, iolet['coord'][2] + 5, 
                iolet['name'], color='black', fontsize=12, weight='bold', zorder=11)
    
    ax.set_title(f"Hemodynamic Dataset Sanity Check{title_suffix}", fontsize=14, weight='bold')
    ax.set_xlabel("Lattice X (Voxels)")
    ax.set_ylabel("Lattice Y (Voxels)")
    ax.set_zlabel("Lattice Z (Voxels)")
    
    # Clean legend entries
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper left', markerscale=1.3)
    
    # Auto-adjust camera view to point directly into the open cutaway channel
    if ENABLE_CUTAWAY and CLIPPING_AXIS == 'Y':
        ax.view_init(elev=25, azim=-45)
        
    fig.colorbar(sc, ax=ax, shrink=0.55, label='Velocity Magnitude (Lattice Units)')
    plt.show()

if __name__ == "__main__":
    if not CSV_PATH.exists():
        print(f"Error: Target data matrix not found at: {CSV_PATH}")
        exit(1)
        
    grids, vels, pressures, timesteps = load_and_parse_hlb(CSV_PATH)
    origin, voxel_size, iolet_markers = parse_xml_metadata(XML_PATH)
    
    final_ts = np.max(timesteps)
    print(f"Filtering variables for the final step: frame {final_ts}")
    idx = (timesteps == final_ts)
    
    X, V, P = grids[idx], vels[idx], pressures[idx]
    vel_magnitude = np.linalg.norm(V, axis=1)
    
    plot_with_matplotlib(X, V, P, vel_magnitude, origin, voxel_size, iolet_markers)