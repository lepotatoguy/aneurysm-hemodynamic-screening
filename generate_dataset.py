# import os
# import subprocess
# import logging
# import traceback
# import numpy as np
# import trimesh
# import networkx as nx
# from pathlib import Path

# # ==========================================
# # CONFIGURATION - UPDATE THESE PATHS
# # ==========================================
# HEMELB_BIN = "hemelb"                 # Documented binary path name (cite: 313)
# SETUP_TOOL_BIN = "hlb-gmy-cli"         # Documented CLI geometry tool (cite: 535)
# HLB_DUMP_BIN = "hlb-dump-extracted-properties" # Documented extraction tool (cite: 344)
# MPI_CORES = 4                         # Adjust based on your local hardware (cite: 312)

# # Directory Setup
# BASE_DIR = Path(__file__).parent
# RAW_DIR = BASE_DIR / "data/raw_meshes"
# GMY_DIR = BASE_DIR / "data/processed_gmy"
# OUT_DIR = BASE_DIR / "data/outputs_csv"
# TEMPLATE_XML = BASE_DIR / "templates/base_input.xml"

# logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')

# def to_hex_float(val):
#     """Converts a standard Python float to an exact hexadecimal float string required by HemeLB (cite: 397)."""
#     return float.hex(float(val))

# def auto_generate_pr2(mesh_path, pr2_path):
#     """Parses an open STL mesh, extracts boundaries, and generates a valid .pr2 profile[cite: 395]."""
#     logging.info(f"Analyzing mesh topology for: {mesh_path.name}")
    
#     # Load the vascular geometry surface
#     mesh = trimesh.load(mesh_path)
    
#     # Isolate boundary edges: rows in mesh.edges_sorted that appear exactly once
#     edges = mesh.edges_sorted
#     unique_edges, counts = np.unique(edges, axis=0, return_counts=True)
#     boundary_edges = unique_edges[counts == 1]
    
#     if len(boundary_edges) == 0:
#         raise ValueError(
#             f"Mesh {mesh_path.name} has no open boundaries! It is completely watertight. "
#             f"HemeLB requires open holes where inlets and outlets can be placed[cite: 540]."
#         )

#     # Clean group boundary edges into discrete closed loops using networkx
#     g = nx.Graph()
#     g.add_edges_from(boundary_edges)
#     loops = list(nx.connected_components(g))
    
#     if len(loops) < 2:
#         raise ValueError(f"Mesh {mesh_path.name} must have at least an inlet and an outlet boundary hole!")

#     # Establish an interior fluid seed point using the mesh center of mass [cite: 449]
#     centroid = mesh.center_mass
    
#     pr2_content = []
#     pr2_content.append(f"Duration Seconds: {to_hex_float(6.0)}") # cite: 473
#     pr2_content.append("Iolets:")
    
#     # Iterate through each detected boundary hole to calculate iolet properties
#     for idx, loop in enumerate(loops):
#         loop_nodes = list(loop)
#         loop_vertices = mesh.vertices[loop_nodes]
        
#         # Calculate geometric center of the boundary disc [cite: 436]
#         center = np.mean(loop_vertices, axis=0)
        
#         # Calculate approximate disc radius [cite: 438]
#         radius = np.max(np.linalg.norm(loop_vertices - center, axis=1)) + 0.5 # Add padding to intersect walls [cite: 276]
        
#         # Determine the principal plane normal of the loop using SVD
#         _, _, vh = np.linalg.svd(loop_vertices - center)
#         normal = vh[2, :]
        
#         # Ensure normal vector points inward toward the fluid domain center of mass [cite: 437]
#         if np.dot(normal, centroid - center) < 0:
#             normal = -normal
            
#         # Classify the primary largest hole as the Inlet, others as Outlets [cite: 434]
#         iolet_type = "Inlet" if idx == 0 else "Outlet" # cite: 434
#         iolet_name = f"{iolet_type}{idx+1}" # cite: 424
        
#         # Set boundary conditions: non-zero target pressure for inlet, zero for outlet baseline [cite: 485, 500]
#         pressure_x = 16.0 if iolet_type == "Inlet" else 0.0
        
#         pr2_content.append(f"  - Centre:") # cite: 420
#         pr2_content.append(f"      x: {to_hex_float(center[0])}")
#         pr2_content.append(f"      y: {to_hex_float(center[1])}")
#         pr2_content.append(f"      z: {to_hex_float(center[2])}")
#         pr2_content.append(f"    Name: {iolet_name}") # cite: 424
#         pr2_content.append(f"    Normal:") # cite: 424
#         pr2_content.append(f"      x: {to_hex_float(normal[0])}")
#         pr2_content.append(f"      y: {to_hex_float(normal[1])}")
#         pr2_content.append(f"      z: {to_hex_float(normal[2])}")
#         pr2_content.append(f"    Pressure:") # cite: 428
#         pr2_content.append(f"      x: {to_hex_float(pressure_x)}")
#         pr2_content.append(f"      y: {to_hex_float(0.0)}")
#         pr2_content.append(f"      z: {to_hex_float(1.0)}")
#         pr2_content.append(f"    Radius: {to_hex_float(radius)}") # cite: 433
#         pr2_content.append(f"    Type: {iolet_type}") # cite: 434

#     # Add downstream file configurations using clean f-strings [cite: 403-412]
#     gmy_out_name = f"{mesh_path.with_suffix('').name}.gmy"
#     xml_out_name = f"{mesh_path.with_suffix('').name}_input.xml"
    
#     pr2_content.append(f"OutputGeometryFile: {gmy_out_name}") # cite: 403
#     pr2_content.append(f"OutputXmlFile: {xml_out_name}") # cite: 404
#     pr2_content.append(f"SeedPoint:") # cite: 405
#     pr2_content.append(f"  x: {to_hex_float(centroid[0])}")
#     pr2_content.append(f"  y: {to_hex_float(centroid[1])}")
#     pr2_content.append(f"  z: {to_hex_float(centroid[2])}")
#     pr2_content.append(f"StlFile: {mesh_path.name}") # cite: 409
#     pr2_content.append(f"StlFileUnitId: 1") # 1 = Millimeters [cite: 410, 512]
#     pr2_content.append(f"TimeStepSeconds: {to_hex_float(1e-5)}") # cite: 411, 513
#     pr2_content.append(f"VoxelSize: {to_hex_float(0.05)}") # Default high resolution [cite: 412, 514]

#     # Save out the complete YAML string profile
#     with open(pr2_path, "w") as f:
#         f.write("\n".join(pr2_content))
#     logging.info(f"Successfully auto-generated profile layout at: {pr2_path.name}")

# def run_cmd(cmd, step_name):
#     """Executes background terminal shell steps safely."""
#     logging.info(f"Starting: {step_name}")
#     try:
#         subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         logging.info(f"Success: {step_name}")
#     except subprocess.CalledProcessError as e:
#         logging.error(f"FAILED: {step_name}\nError details:\n{e.stderr.decode('utf-8')}")
#         raise

# def process_patient(mesh_file):
#     patient_id = mesh_file.stem
#     logging.info(f"============ Processing Case: {patient_id} ============")
    
#     # Path initializations sharing the identical execution directory folder context (cite: 536)
#     patient_pr2 = RAW_DIR / f"{patient_id}.pr2"
#     patient_gmy = GMY_DIR / f"{patient_id}.gmy"
#     patient_xml = GMY_DIR / f"{patient_id}_input.xml"
#     patient_out_dir = OUT_DIR / patient_id
    
#     # gmy_out_name = mesh_file.with_suffix(".gmy").name
#     # xml_out_name = mesh_file.with_suffix("_input.xml").name
#     gmy_out_name = f"{mesh_file.stem}.gmy"
#     xml_out_name = f"{mesh_file.stem}_input.xml"
    
#     patient_out_dir.mkdir(parents=True, exist_ok=True)

#     # STEP 1: Topological Feature Extraction & Profile Auto-Writing
#     auto_generate_pr2(mesh_file, patient_pr2)

#     # STEP 2: Non-GUI Command-Line Voxelization
#     # The setup tool expects to run directly inside the folder where the files sit (cite: 536)
#     setup_cmd = f"cd {RAW_DIR} && {SETUP_TOOL_BIN} {patient_pr2.name}" # cite: 535
#     run_cmd(setup_cmd, f"Headless Voxelization Loop ({patient_id})")

#     # Move generated assets into their organized tracking directories
#     os.rename(RAW_DIR / gmy_out_name, patient_gmy)
#     os.rename(RAW_DIR / xml_out_name, patient_xml)

#     # STEP 3: Run Simulation via Compiled Parallel Binary
#     hemelb_cmd = f"mpirun -n {MPI_CORES} {HEMELB_BIN} in {patient_xml} -out {patient_out_dir}" # cite: 312
#     run_cmd(hemelb_cmd, f"HemeLB Core Simulation ({patient_id})")

#     # STEP 4: High-Throughput Matrix Feature Extraction to CSV
#     xtr_path = patient_out_dir / "whole.xtr" # cite: 302
#     csv_path = patient_out_dir / f"{patient_id}_fluid_data.csv"
    
#     extract_cmd = f"{HLB_DUMP_BIN} {xtr_path} > {csv_path}" # cite: 344
#     run_cmd(extract_cmd, f"Compiling Output Features Matrix ({patient_id})")
    
#     logging.info(f"============ Case {patient_id} Successfully Generated! ============\n")

# if __name__ == "__main__":
#     GMY_DIR.mkdir(parents=True, exist_ok=True)
#     OUT_DIR.mkdir(parents=True, exist_ok=True)

#     # Support all native medical file representations interchangeably (cite: 357)
#     mesh_files = list(RAW_DIR.glob("*.stl")) + list(RAW_DIR.glob("*.vtp"))
    
#     if not mesh_files:
#         logging.warning(f"No surface files found in {RAW_DIR}. Drop your downloaded case file there.")
    
#     for mesh in mesh_files:
#         try:
#             process_patient(mesh)
#         except Exception as e:
#             logging.error(f"Skipping case tracking loop for {mesh.name} due to unexpected execution errors.")
#             logging.error(traceback.format_exc()) # Print full traceback log to target issues precisely
#             continue




###############################################

# import os
# import subprocess
# import logging
# import traceback
# import numpy as np
# import trimesh
# import networkx as nx
# from pathlib import Path

# # ==========================================
# # CONFIGURATION - UPDATE THESE PATHS
# # ==========================================
# HEMELB_BIN = "hemelb"                 # Documented binary path name (cite: 313)
# SETUP_TOOL_BIN = "hlb-gmy-cli"         # Documented CLI geometry tool (cite: 535)
# HLB_DUMP_BIN = "hlb-dump-extracted-properties" # Documented extraction tool (cite: 344)
# MPI_CORES = 4                         # Adjust based on your local hardware (cite: 312)

# # Directory Setup
# BASE_DIR = Path(__file__).parent
# RAW_DIR = BASE_DIR / "data/raw_meshes"
# GMY_DIR = BASE_DIR / "data/processed_gmy"
# OUT_DIR = BASE_DIR / "data/outputs_csv"
# TEMPLATE_XML = BASE_DIR / "templates/base_input.xml"

# logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')

# def to_hex_float(val):
#     """Converts a standard Python float to an exact hexadecimal float string required by HemeLB (cite: 397)."""
#     return float.hex(float(val))

# def auto_generate_pr2(mesh_path, pr2_path):
#     """Parses an open STL mesh, extracts boundaries, and generates a valid .pr2 profile (cite: 395)."""
#     logging.info(f"Analyzing mesh topology for: {mesh_path.name}")
    
#     # Load the vascular geometry surface
#     mesh = trimesh.load(mesh_path)
    
#     # Isolate boundary edges: rows in mesh.edges_sorted that appear exactly once
#     edges = mesh.edges_sorted
#     unique_edges, counts = np.unique(edges, axis=0, return_counts=True)
#     boundary_edges = unique_edges[counts == 1]
    
#     if len(boundary_edges) == 0:
#         raise ValueError(
#             f"Mesh {mesh_path.name} has no open boundaries! It is completely watertight. "
#             f"HemeLB requires open holes where inlets and outlets can be placed."
#         )

#     # Clean group boundary edges into discrete closed loops using networkx
#     g = nx.Graph()
#     g.add_edges_from(boundary_edges)
#     loops = list(nx.connected_components(g))
    
#     if len(loops) < 2:
#         raise ValueError(f"Mesh {mesh_path.name} must have at least an inlet and an outlet boundary hole!")

#     # Establish an interior fluid seed point using the mesh center of mass (cite: 449)
#     centroid = mesh.center_mass
    
#     pr2_content = []
#     pr2_content.append(f"Duration Seconds: {to_hex_float(6.0)}") # cite: 473
#     pr2_content.append("Iolets:") # Changing lowercase 'lolets' to uppercase 'Iolets' fixes the KeyError
    
#     # Iterate through each detected boundary hole to calculate iolet properties
#     for idx, loop in enumerate(loops):
#         loop_nodes = list(loop)
#         loop_vertices = mesh.vertices[loop_nodes]
        
#         # Calculate geometric center of the boundary disc (cite: 436)
#         center = np.mean(loop_vertices, axis=0)
        
#         # Calculate approximate disc radius (cite: 438)
#         radius = np.max(np.linalg.norm(loop_vertices - center, axis=1)) + 0.5 # Add padding to intersect walls (cite: 276)
        
#         # Determine the principal plane normal of the loop using SVD
#         _, _, vh = np.linalg.svd(loop_vertices - center)
#         normal = vh[2, :]
        
#         # Ensure normal vector points inward toward the fluid domain center of mass (cite: 437)
#         if np.dot(normal, centroid - center) < 0:
#             normal = -normal
            
#         # Classify the primary largest hole as the Inlet, others as Outlets
#         iolet_type = "Inlet" if idx == 0 else "Outlet" # cite: 434
#         iolet_name = f"{iolet_type}{idx+1}" # cite: 424
        
#         # Set boundary conditions: non-zero target pressure for inlet, zero for outlet baseline (cite: 485, 500)
#         pressure_x = 16.0 if iolet_type == "Inlet" else 0.0
        
#         pr2_content.append(f"  - Centre:") # cite: 420
#         pr2_content.append(f"      x: {to_hex_float(center[0])}")
#         pr2_content.append(f"      y: {to_hex_float(center[1])}")
#         pr2_content.append(f"      z: {to_hex_float(center[2])}")
#         pr2_content.append(f"    Name: {iolet_name}") # cite: 424
#         pr2_content.append(f"    Normal:") # cite: 424
#         pr2_content.append(f"      x: {to_hex_float(normal[0])}")
#         pr2_content.append(f"      y: {to_hex_float(normal[1])}")
#         pr2_content.append(f"      z: {to_hex_float(normal[2])}")
#         pr2_content.append(f"    Pressure:") # cite: 428
#         pr2_content.append(f"      x: {to_hex_float(pressure_x)}")
#         pr2_content.append(f"      y: {to_hex_float(0.0)}")
#         pr2_content.append(f"      z: {to_hex_float(1.0)}")
#         pr2_content.append(f"    Radius: {to_hex_float(radius)}") # cite: 433
#         pr2_content.append(f"    Type: {iolet_type}") # cite: 434

#     # Add downstream file configurations using clean f-strings (cite: 403-412)
#     gmy_out_name = f"{mesh_path.with_suffix('').name}.gmy"
#     xml_out_name = f"{mesh_path.with_suffix('').name}_input.xml"
    
#     pr2_content.append(f"OutputGeometryFile: {gmy_out_name}") # cite: 403
#     pr2_content.append(f"OutputXmlFile: {xml_out_name}") # cite: 404
#     pr2_content.append(f"SeedPoint:") # cite: 405
#     pr2_content.append(f"  x: {to_hex_float(centroid[0])}")
#     pr2_content.append(f"  y: {to_hex_float(centroid[1])}")
#     pr2_content.append(f"  z: {to_hex_float(centroid[2])}")
#     pr2_content.append(f"StlFile: {mesh_path.name}") # cite: 409
#     pr2_content.append(f"StlFileUnitId: 1") # 1 = Millimeters (cite: 512)
#     pr2_content.append(f"TimeStepSeconds: {to_hex_float(1e-5)}") # cite: 411, 513
#     pr2_content.append(f"VoxelSize: {to_hex_float(0.05)}") # Default high resolution (cite: 514)

#     # Save out the complete YAML string profile
#     with open(pr2_path, "w") as f:
#         f.write("\n".join(pr2_content))
#     logging.info(f"Successfully auto-generated profile layout at: {pr2_path.name}")

# def run_cmd(cmd, step_name):
#     """Executes background terminal shell steps with robust error dumping."""
#     logging.info(f"Starting: {step_name}")
#     try:
#         subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         logging.info(f"Success: {step_name}")
#     except subprocess.CalledProcessError as e:
#         # Crucial fix: capturing BOTH stdout and stderr to view internal HemeLB logging errors
#         stdout_msg = e.stdout.decode('utf-8') if e.stdout else "None"
#         stderr_msg = e.stderr.decode('utf-8') if e.stderr else "None"
#         logging.error(f"FAILED: {step_name}\n\n[STANDARD ERROR]:\n{stderr_msg}\n\n[STANDARD OUTPUT]:\n{stdout_msg}")
#         raise

# def process_patient(mesh_file):
#     patient_id = mesh_file.stem
#     logging.info(f"============ Processing Case: {patient_id} ============")
    
#     # Path initializations sharing the identical execution directory folder context (cite: 536)
#     patient_pr2 = RAW_DIR / f"{patient_id}.pr2"
#     patient_gmy = GMY_DIR / f"{patient_id}.gmy"
#     patient_xml = GMY_DIR / f"{patient_id}_input.xml"
#     patient_out_dir = OUT_DIR / patient_id
    
#     gmy_out_name = f"{mesh_file.stem}.gmy"
#     xml_out_name = f"{mesh_file.stem}_input.xml"
    
#     import shutil
#     if patient_out_dir.exists():
#         logging.info(f"Purging stale target output folder directory: {patient_out_dir.name}")
#         shutil.rmtree(patient_out_dir)
#     patient_out_dir.mkdir(parents=True, exist_ok=True)

#     # STEP 1: Topological Feature Extraction & Profile Auto-Writing
#     auto_generate_pr2(mesh_file, patient_pr2)

#     # STEP 2: Non-GUI Command-Line Voxelization (cite: 536)
#     setup_cmd = f"cd {RAW_DIR} && {SETUP_TOOL_BIN} {patient_pr2.name}" # cite: 535
#     run_cmd(setup_cmd, f"Headless Voxelization Loop ({patient_id})")

#     # Move the voxelized grid geometry file to its proper tracking directory
#     os.rename(RAW_DIR / gmy_out_name, patient_gmy)
    
#     # Clean up the bare-bones XML generated by the setup tool
#     if (RAW_DIR / xml_out_name).exists():
#         os.remove(RAW_DIR / xml_out_name)

#     # STEP 3: Inject custom properties block and ABSOLUTE path into the XML template
#     with open(TEMPLATE_XML, 'r') as file:
#         xml_content = file.read()
    
#     # Using .resolve() forces an absolute path so HemeLB can discover the .gmy file from any directory (cite: 292)
#     xml_content = xml_content.format(gmy_path=str(patient_gmy.resolve()))
    
#     with open(patient_xml, 'w') as file:
#         file.write(xml_content)
#     logging.info(f"Generated custom production XML configuration for {patient_id}")

#     # STEP 4: Run Simulation via Compiled Parallel Binary
#     hemelb_cmd = f"mpirun -n {MPI_CORES} {HEMELB_BIN} -in {patient_xml} -out {patient_out_dir}" # cite: 312
#     run_cmd(hemelb_cmd, f"HemeLB Core Simulation ({patient_id})")

#     # STEP 5: High-Throughput Matrix Feature Extraction to CSV
#     xtr_path = patient_out_dir / "whole.xtr" # cite: 302
#     csv_path = patient_out_dir / f"{patient_id}_fluid_data.csv"
    
#     extract_cmd = f"{HLB_DUMP_BIN} {xtr_path} > {csv_path}" # cite: 344
#     run_cmd(extract_cmd, f"Compiling Output Features Matrix ({patient_id})")
    
#     logging.info(f"============ Case {patient_id} Successfully Generated! ============\n")

# if __name__ == "__main__":
#     GMY_DIR.mkdir(parents=True, exist_ok=True)
#     OUT_DIR.mkdir(parents=True, exist_ok=True)

#     # Support all native medical file representations interchangeably (cite: 357)
#     mesh_files = list(RAW_DIR.glob("*.stl")) + list(RAW_DIR.glob("*.vtp"))
    
#     if not mesh_files:
#         logging.warning(f"No surface files found in {RAW_DIR}. Drop your downloaded case file there.")
    
#     for mesh in mesh_files:
#         try:
#             process_patient(mesh)
#         except Exception as e:
#             logging.error(f"Skipping case tracking loop for {mesh.name} due to unexpected execution errors.")
#             logging.error(traceback.format_exc())
#             continue

import os
import subprocess
import logging
import traceback
import shutil
import xml.etree.ElementTree as ET
import numpy as np
import trimesh
import networkx as nx
from pathlib import Path

# ==========================================
# CONFIGURATION - UPDATE THESE PATHS
# ==========================================
HEMELB_BIN = "hemelb"                 # Documented binary path name [cite: 703, 884]
SETUP_TOOL_BIN = "hlb-gmy-cli"         # Documented CLI geometry tool [cite: 795]
HLB_DUMP_BIN = "hlb-dump-extracted-properties" # Documented extraction tool [cite: 915]
MPI_CORES = 4                         # Adjust based on your local hardware [cite: 862, 883]

# Directory Setup
BASE_DIR = Path(__file__).parent
RAW_DIR = BASE_DIR / "data/raw_meshes"
GMY_DIR = BASE_DIR / "data/processed_gmy"
OUT_DIR = BASE_DIR / "data/outputs_csv"

logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')

def auto_generate_pr2(mesh_path, pr2_path):
    """Parses an open STL mesh, extracts boundaries, filters out topological noise, and generates a valid .pr2 profile."""
    logging.info(f"Withdrawing mesh topology markers for: {mesh_path.name}")
    
    # Load the vascular geometry surface
    mesh = trimesh.load(mesh_path)
    
    # Isolate boundary edges: rows in mesh.edges_sorted that appear exactly once
    edges = mesh.edges_sorted
    unique_edges, counts = np.unique(edges, axis=0, return_counts=True)
    boundary_edges = unique_edges[counts == 1]
    
    if len(boundary_edges) == 0:
        raise ValueError(f"Mesh {mesh_path.name} has no open boundaries! HemeLB requires open holes [cite: 1110-1112].")

    # Clean group boundary edges into discrete closed loops using networkx
    g = nx.Graph()
    g.add_edges_from(boundary_edges)
    all_loops = list(nx.connected_components(g))
    
    # TOPO FILTERING: Discard microscopic cracks/artifacts by requiring a minimum vertex count per loop
    loops = [l for l in all_loops if len(l) >= 20]
    
    if len(loops) < 2:
        logging.warning("Topological filtering was too aggressive. Falling back to the largest available components.")
        all_loops.sort(key=len, reverse=True)
        loops = all_loops[:2]

    logging.info(f"Identified {len(loops)} genuine fluid boundaries after filtering out topological mesh noise.")

    # Establish an interior fluid seed point using the mesh center of mass [cite: 848, 1020-1021]
    centroid = mesh.center_mass
    
    pr2_content = []
    pr2_content.append("DurationSeconds: 6.0") # Spaceless key matching verified profile format [cite: 1044]
    pr2_content.append("Iolets:") # Uppercase key matching verified profile format [cite: 1045]
    
    inlet_counter = 1
    outlet_counter = 1
    
    # Sort remaining valid loops by size (vertex count) to easily isolate the major inlet opening
    loops.sort(key=len, reverse=True)
    
    # Iterate through each verified boundary hole to calculate iolet properties
    for idx, loop in enumerate(loops):
        loop_nodes = list(loop)
        loop_vertices = mesh.vertices[loop_nodes]
        
        # Calculate geometric center of the boundary disc [cite: 1007]
        center = np.mean(loop_vertices, axis=0)
        
        # Padding clears wall boundary voxels smoothly at a 0.1 mm grid size [cite: 847, 1086]
        radius = np.max(np.linalg.norm(loop_vertices - center, axis=1)) + 0.1
        
        # Determine the principal plane normal of the loop using SVD
        _, _, vh = np.linalg.svd(loop_vertices - center)
        normal = vh[2, :]
        
        # Ensure normal vector points inward toward the fluid domain center of mass [cite: 1008]
        if np.dot(normal, centroid - center) < 0:
            normal = -normal
            
        # Classify boundaries sequentially starting from index 1 [cite: 1050, 1065]
        if idx == 0:
            iolet_type = "Inlet"
            iolet_name = f"Inlet{inlet_counter}"
            inlet_counter += 1
            # Lowering this from 16.0 to 0.1 or 1.0 mmHg ensures LBM numerical stability 
            # in complex meshes while testing at a coarse prototyping resolution.
            pressure_x = 1.0
        else:
            iolet_type = "Outlet"
            iolet_name = f"Outlet{outlet_counter}"
            outlet_counter += 1
            pressure_x = 0.0   # Outlets ground to zero baseline [cite: 1071]
        
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
        pr2_content.append("    z: 0.0")  # Phase offset initialization
        pr2_content.append(f"  Radius: {radius:.10f}")
        pr2_content.append(f"  Type: {iolet_type}")

    gmy_out_name = f"{mesh_path.with_suffix('').name}.gmy"
    xml_out_name = f"{mesh_path.with_suffix('').name}_input.xml"
    
    pr2_content.append(f"OutputGeometryFile: {gmy_out_name}") #[cite: 1076]
    pr2_content.append(f"OutputXmlFile: {xml_out_name}") #[cite: 1077]
    pr2_content.append("SeedPoint:") #[cite: 1078]
    pr2_content.append(f"  x: {centroid[0]:.10f}")
    pr2_content.append(f"  y: {centroid[1]:.10f}")
    pr2_content.append(f"  z: {centroid[2]:.10f}")
    pr2_content.append(f"StlFile: {mesh_path.name}") #[cite: 1082]
    pr2_content.append("StlFileUnitId: 1") # 1 = Millimeters [cite: 1083]
    pr2_content.append("TimeStepSeconds: 0.00001") #[cite: 1084]
    pr2_content.append("VoxelSize: 0.1") # High-speed stable grid optimization [cite: 1086]

    with open(pr2_path, "w") as f:
        f.write("\n".join(pr2_content))
    logging.info(f"Successfully auto-generated profile layout at: {pr2_path.name}")

def run_cmd(cmd, step_name):
    """Executes background terminal shell steps with robust error dumping."""
    logging.info(f"Starting: {step_name}")
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.info(f"Success: {step_name}")
    except subprocess.CalledProcessError as e:
        stdout_msg = e.stdout.decode('utf-8') if e.stdout else "None"
        stderr_msg = e.stderr.decode('utf-8') if e.stderr else "None"
        logging.error(f"FAILED: {step_name}\n\n[STANDARD ERROR]:\n{stderr_msg}\n\n[STANDARD OUTPUT]:\n{stdout_msg}")
        raise

def process_patient(mesh_file):
    patient_id = mesh_file.stem
    logging.info(f"============ Processing Case: {patient_id} ============")
    
    patient_pr2 = RAW_DIR / f"{patient_id}.pr2"
    patient_gmy = GMY_DIR / f"{patient_id}.gmy"
    patient_xml = GMY_DIR / f"{patient_id}_input.xml"
    patient_out_dir = OUT_DIR / patient_id
    
    xml_out_name = f"{mesh_file.stem}_input.xml"
    gmy_out_name = f"{mesh_file.stem}.gmy"

    # STEP 1: Topological Feature Extraction & Profile Auto-Writing
    auto_generate_pr2(mesh_file, patient_pr2)

    # STEP 2: Non-GUI Command-Line Voxelization [cite: 795]
    setup_cmd = f"cd {RAW_DIR} && {SETUP_TOOL_BIN} {patient_pr2.name}"
    run_cmd(setup_cmd, f"Headless Voxelization Loop ({patient_id})")

    # Move the natively voxelized grid geometry file to its proper tracking directory
    os.rename(RAW_DIR / gmy_out_name, patient_gmy)
    
    # STEP 3: Patch and Modify Native Version 5 XML Layout directly [cite: 930]
    generated_xml_path = RAW_DIR / xml_out_name
    logging.info(f"Patching machine learning tracking blocks directly into native Version 5 XML")
    
    tree = ET.parse(generated_xml_path)
    root = tree.getroot()
    
    # Update datafile path to use absolute resolved pathway [cite: 863]
    datafile = root.find(".//geometry/datafile")
    if datafile is not None:
        datafile.set("path", str(patient_gmy.resolve()))
        
    # Accelerated prototyping constraint: 1000 iteration run-time limit
    steps_element = root.find(".//simulation/steps")
    if steps_element is not None:
        steps_element.set("value", "1000")
        
    # Append custom machine learning properties extraction block [cite: 872-878]
    properties = root.find("properties")
    if properties is None:
        properties = ET.SubElement(root, "properties")
        prop_output = ET.SubElement(properties, "propertyoutput", {"file": "whole.xtr", "period": "100"})
        ET.SubElement(prop_output, "geometry", {"type": "whole"})
        ET.SubElement(prop_output, "field", {"type": "velocity"})
        ET.SubElement(prop_output, "field", {"type": "pressure"})
        
    # Save modified tree to processing path
    tree.write(patient_xml, encoding="utf-8", xml_declaration=True)
    
    # Clean up temporary version 5 XML asset file
    if generated_xml_path.exists():
        os.remove(generated_xml_path)

    # STEP 4: Flush old output directory paths so the v5 engine doesn't encounter file collisions [cite: 864]
    if patient_out_dir.exists():
        logging.info(f"Purging pre-existing output directory to clear execution checks: {patient_out_dir.name}")
        shutil.rmtree(patient_out_dir)

    # DIAGNOSTIC STEP: Print out the complete generated XML file to console right before running
    try:
        with open(patient_xml, 'r') as f:
            logging.info(f"\n==================================================\n--- DIAGNOSTIC: GENERATED XML CONTENT ---\n==================================================\n{f.read()}\n==================================================")
    except Exception as e:
        logging.error(f"Could not read generated XML for diagnostic printing: {e}")

    # STEP 5: Run Simulation via Parallel Engine [cite: 860]
    hemelb_cmd = f"mpirun -n {MPI_CORES} {HEMELB_BIN} -in {patient_xml} -out {patient_out_dir}"
    run_cmd(hemelb_cmd, f"HemeLB Core Simulation ({patient_id})")

    # # STEP 6: High-Throughput Matrix Feature Extraction to CSV [cite: 915]
    # xtr_path = patient_out_dir / "whole.xtr"
    # csv_path = patient_out_dir / f"{patient_id}_fluid_data.csv"
    
    # extract_cmd = f"{HLB_DUMP_BIN} {xtr_path} > {csv_path}"
    # run_cmd(extract_cmd, f"Compiling Output Features Matrix ({patient_id})")
    
    # logging.info(f"============ Case {patient_id} Successfully Generated! ============\n")

    # STEP 6: High-Throughput Matrix Feature Extraction to CSV (cite: 344)
    # Dynamic search handles HemeLB's 'Extracted/' subfolder layout automatically
    xtr_matches = list(patient_out_dir.glob("**/whole.xtr"))
    
    if not xtr_matches:
        # Diagnostic fallback prints files if something goes wrong
        existing_files = [str(p.relative_to(BASE_DIR)) for p in patient_out_dir.rglob("*") if p.is_file()]
        raise FileNotFoundError(
            f"Could not find 'whole.xtr' inside {patient_out_dir}.\n"
            f"Files actually generated by HemeLB:\n" + "\n".join(existing_files)
        )
        
    xtr_path = xtr_matches[0]
    csv_path = patient_out_dir / f"{patient_id}_fluid_data.csv"
    
    extract_cmd = f"{HLB_DUMP_BIN} {xtr_path} > {csv_path}" # cite: 344
    run_cmd(extract_cmd, f"Compiling Output Features Matrix ({patient_id})")
    
    logging.info(f"============ Case {patient_id} Successfully Generated! ============\n")

if __name__ == "__main__":
    GMY_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Support all native medical file representations interchangeably [cite: 882]
    mesh_files = list(RAW_DIR.glob("*.stl")) + list(RAW_DIR.glob("*.vtp"))
    
    if not mesh_files:
        logging.warning(f"No surface files found in {RAW_DIR}. Drop your downloaded case file there.")
    
    for mesh in mesh_files:
        try:
            process_patient(mesh)
        except Exception as e:
            logging.error(f"Skipping case tracking loop for {mesh.name} due to unexpected execution errors.")
            logging.error(traceback.format_exc())
            continue