#!/usr/bin/env python3


# import os
# import subprocess
# import logging
# import traceback
# import shutil
# import xml.etree.ElementTree as ET
# import numpy as np
# import trimesh
# import networkx as nx
# from pathlib import Path

# # ==========================================
# # CONFIGURATION - UPDATE THESE PATHS
# # ==========================================
# HEMELB_BIN = "hemelb"                 # Documented binary path name [cite: 703, 884]
# SETUP_TOOL_BIN = "hlb-gmy-cli"         # Documented CLI geometry tool [cite: 795]
# HLB_DUMP_BIN = "hlb-dump-extracted-properties" # Documented extraction tool [cite: 915]
# MPI_CORES = 4                         # Adjust based on your local hardware [cite: 862, 883]

# # Directory Setup
# BASE_DIR = Path(__file__).parent
# RAW_DIR = BASE_DIR / "data/raw_meshes"
# GMY_DIR = BASE_DIR / "data/processed_gmy"
# OUT_DIR = BASE_DIR / "data/outputs_csv"

# logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')

# def auto_generate_pr2(mesh_path, pr2_path):
#     """Parses an open STL mesh, extracts boundaries, filters out topological noise, and generates a valid .pr2 profile."""
#     logging.info(f"Withdrawing mesh topology markers for: {mesh_path.name}")
    
#     # Load the vascular geometry surface
#     mesh = trimesh.load(mesh_path)
    
#     # Isolate boundary edges: rows in mesh.edges_sorted that appear exactly once
#     edges = mesh.edges_sorted
#     unique_edges, counts = np.unique(edges, axis=0, return_counts=True)
#     boundary_edges = unique_edges[counts == 1]
    
#     if len(boundary_edges) == 0:
#         raise ValueError(f"Mesh {mesh_path.name} has no open boundaries! HemeLB requires open holes [cite: 1110-1112].")

#     # Clean group boundary edges into discrete closed loops using networkx
#     g = nx.Graph()
#     g.add_edges_from(boundary_edges)
#     all_loops = list(nx.connected_components(g))
    
#     # TOPO FILTERING: Discard microscopic cracks/artifacts by requiring a minimum vertex count per loop
#     loops = [l for l in all_loops if len(l) >= 20]
    
#     if len(loops) < 2:
#         logging.warning("Topological filtering was too aggressive. Falling back to the largest available components.")
#         all_loops.sort(key=len, reverse=True)
#         loops = all_loops[:2]

#     logging.info(f"Identified {len(loops)} genuine fluid boundaries after filtering out topological mesh noise.")

#     # Establish an interior fluid seed point using the mesh center of mass [cite: 848, 1020-1021]
#     centroid = mesh.center_mass
    
#     pr2_content = []
#     pr2_content.append("DurationSeconds: 6.0") # Spaceless key matching verified profile format [cite: 1044]
#     pr2_content.append("Iolets:") # Uppercase key matching verified profile format [cite: 1045]
    
#     inlet_counter = 1
#     outlet_counter = 1
    
#     # Sort remaining valid loops by size (vertex count) to easily isolate the major inlet opening
#     loops.sort(key=len, reverse=True)
    
#     # Iterate through each verified boundary hole to calculate iolet properties
#     for idx, loop in enumerate(loops):
#         loop_nodes = list(loop)
#         loop_vertices = mesh.vertices[loop_nodes]
        
#         # Calculate geometric center of the boundary disc [cite: 1007]
#         center = np.mean(loop_vertices, axis=0)
        
#         # Padding clears wall boundary voxels smoothly at a 0.1 mm grid size [cite: 847, 1086]
#         radius = np.max(np.linalg.norm(loop_vertices - center, axis=1)) + 0.1
        
#         # Determine the principal plane normal of the loop using SVD
#         _, _, vh = np.linalg.svd(loop_vertices - center)
#         normal = vh[2, :]
        
#         # Ensure normal vector points inward toward the fluid domain center of mass [cite: 1008]
#         if np.dot(normal, centroid - center) < 0:
#             normal = -normal
            
#         # Classify boundaries sequentially starting from index 1 [cite: 1050, 1065]
#         if idx == 0:
#             iolet_type = "Inlet"
#             iolet_name = f"Inlet{inlet_counter}"
#             inlet_counter += 1
#             # Lowering this from 16.0 to 0.1 or 1.0 mmHg ensures LBM numerical stability 
#             # in complex meshes while testing at a coarse prototyping resolution.
#             pressure_x = 1.0
#         else:
#             iolet_type = "Outlet"
#             iolet_name = f"Outlet{outlet_counter}"
#             outlet_counter += 1
#             pressure_x = 0.0   # Outlets ground to zero baseline [cite: 1071]
        
#         pr2_content.append("- Centre:")
#         pr2_content.append(f"    x: {center[0]:.10f}")
#         pr2_content.append(f"    y: {center[1]:.10f}")
#         pr2_content.append(f"    z: {center[2]:.10f}")
#         pr2_content.append(f"  Name: {iolet_name}")
#         pr2_content.append("  Normal:")
#         pr2_content.append(f"    x: {normal[0]:.10f}")
#         pr2_content.append(f"    y: {normal[1]:.10f}")
#         pr2_content.append(f"    z: {normal[2]:.10f}")
#         pr2_content.append("  Pressure:")
#         pr2_content.append(f"    x: {pressure_x:.10f}")
#         pr2_content.append("    y: 0.0")
#         pr2_content.append("    z: 0.0")  # Phase offset initialization
#         pr2_content.append(f"  Radius: {radius:.10f}")
#         pr2_content.append(f"  Type: {iolet_type}")

#     gmy_out_name = f"{mesh_path.with_suffix('').name}.gmy"
#     xml_out_name = f"{mesh_path.with_suffix('').name}_input.xml"
    
#     pr2_content.append(f"OutputGeometryFile: {gmy_out_name}") #[cite: 1076]
#     pr2_content.append(f"OutputXmlFile: {xml_out_name}") #[cite: 1077]
#     pr2_content.append("SeedPoint:") #[cite: 1078]
#     pr2_content.append(f"  x: {centroid[0]:.10f}")
#     pr2_content.append(f"  y: {centroid[1]:.10f}")
#     pr2_content.append(f"  z: {centroid[2]:.10f}")
#     pr2_content.append(f"StlFile: {mesh_path.name}") #[cite: 1082]
#     pr2_content.append("StlFileUnitId: 1") # 1 = Millimeters [cite: 1083]
#     pr2_content.append("TimeStepSeconds: 0.00001") #[cite: 1084]
#     pr2_content.append("VoxelSize: 0.1") # High-speed stable grid optimization [cite: 1086]

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
#         stdout_msg = e.stdout.decode('utf-8') if e.stdout else "None"
#         stderr_msg = e.stderr.decode('utf-8') if e.stderr else "None"
#         logging.error(f"FAILED: {step_name}\n\n[STANDARD ERROR]:\n{stderr_msg}\n\n[STANDARD OUTPUT]:\n{stdout_msg}")
#         raise

# def process_patient(mesh_file):
#     patient_id = mesh_file.stem
#     logging.info(f"============ Processing Case: {patient_id} ============")
    
#     patient_pr2 = RAW_DIR / f"{patient_id}.pr2"
#     patient_gmy = GMY_DIR / f"{patient_id}.gmy"
#     patient_xml = GMY_DIR / f"{patient_id}_input.xml"
#     patient_out_dir = OUT_DIR / patient_id
    
#     xml_out_name = f"{mesh_file.stem}_input.xml"
#     gmy_out_name = f"{mesh_file.stem}.gmy"

#     # STEP 1: Topological Feature Extraction & Profile Auto-Writing
#     auto_generate_pr2(mesh_file, patient_pr2)

#     # STEP 2: Non-GUI Command-Line Voxelization [cite: 795]
#     setup_cmd = f"cd {RAW_DIR} && {SETUP_TOOL_BIN} {patient_pr2.name}"
#     run_cmd(setup_cmd, f"Headless Voxelization Loop ({patient_id})")

#     # Move the natively voxelized grid geometry file to its proper tracking directory
#     os.rename(RAW_DIR / gmy_out_name, patient_gmy)
    
#     # STEP 3: Patch and Modify Native Version 5 XML Layout directly [cite: 930]
#     generated_xml_path = RAW_DIR / xml_out_name
#     logging.info(f"Patching machine learning tracking blocks directly into native Version 5 XML")
    
#     tree = ET.parse(generated_xml_path)
#     root = tree.getroot()
    
#     # Update datafile path to use absolute resolved pathway [cite: 863]
#     datafile = root.find(".//geometry/datafile")
#     if datafile is not None:
#         datafile.set("path", str(patient_gmy.resolve()))
        
#     # Accelerated prototyping constraint: 1000 iteration run-time limit
#     steps_element = root.find(".//simulation/steps")
#     if steps_element is not None:
#         steps_element.set("value", "1000")
        
#     # Append custom machine learning properties extraction block [cite: 872-878]
#     properties = root.find("properties")
#     if properties is None:
#         properties = ET.SubElement(root, "properties")
#         prop_output = ET.SubElement(properties, "propertyoutput", {"file": "whole.xtr", "period": "100"})
#         ET.SubElement(prop_output, "geometry", {"type": "whole"})
#         ET.SubElement(prop_output, "field", {"type": "velocity"})
#         ET.SubElement(prop_output, "field", {"type": "pressure"})
        
#     # Save modified tree to processing path
#     tree.write(patient_xml, encoding="utf-8", xml_declaration=True)
    
#     # Clean up temporary version 5 XML asset file
#     if generated_xml_path.exists():
#         os.remove(generated_xml_path)

#     # STEP 4: Flush old output directory paths so the v5 engine doesn't encounter file collisions [cite: 864]
#     if patient_out_dir.exists():
#         logging.info(f"Purging pre-existing output directory to clear execution checks: {patient_out_dir.name}")
#         shutil.rmtree(patient_out_dir)

#     # DIAGNOSTIC STEP: Print out the complete generated XML file to console right before running
#     try:
#         with open(patient_xml, 'r') as f:
#             logging.info(f"\n==================================================\n--- DIAGNOSTIC: GENERATED XML CONTENT ---\n==================================================\n{f.read()}\n==================================================")
#     except Exception as e:
#         logging.error(f"Could not read generated XML for diagnostic printing: {e}")

#     # STEP 5: Run Simulation via Parallel Engine [cite: 860]
#     hemelb_cmd = f"mpirun -n {MPI_CORES} {HEMELB_BIN} -in {patient_xml} -out {patient_out_dir}"
#     run_cmd(hemelb_cmd, f"HemeLB Core Simulation ({patient_id})")

#     # # STEP 6: High-Throughput Matrix Feature Extraction to CSV [cite: 915]
#     # xtr_path = patient_out_dir / "whole.xtr"
#     # csv_path = patient_out_dir / f"{patient_id}_fluid_data.csv"
    
#     # extract_cmd = f"{HLB_DUMP_BIN} {xtr_path} > {csv_path}"
#     # run_cmd(extract_cmd, f"Compiling Output Features Matrix ({patient_id})")
    
#     # logging.info(f"============ Case {patient_id} Successfully Generated! ============\n")

#     # STEP 6: High-Throughput Matrix Feature Extraction to CSV (cite: 344)
#     # Dynamic search handles HemeLB's 'Extracted/' subfolder layout automatically
#     xtr_matches = list(patient_out_dir.glob("**/whole.xtr"))
    
#     if not xtr_matches:
#         # Diagnostic fallback prints files if something goes wrong
#         existing_files = [str(p.relative_to(BASE_DIR)) for p in patient_out_dir.rglob("*") if p.is_file()]
#         raise FileNotFoundError(
#             f"Could not find 'whole.xtr' inside {patient_out_dir}.\n"
#             f"Files actually generated by HemeLB:\n" + "\n".join(existing_files)
#         )
        
#     xtr_path = xtr_matches[0]
#     csv_path = patient_out_dir / f"{patient_id}_fluid_data.csv"
    
#     extract_cmd = f"{HLB_DUMP_BIN} {xtr_path} > {csv_path}" # cite: 344
#     run_cmd(extract_cmd, f"Compiling Output Features Matrix ({patient_id})")
    
#     logging.info(f"============ Case {patient_id} Successfully Generated! ============\n")

# if __name__ == "__main__":
#     GMY_DIR.mkdir(parents=True, exist_ok=True)
#     OUT_DIR.mkdir(parents=True, exist_ok=True)

#     # Support all native medical file representations interchangeably [cite: 882]
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


###################################################

# """
# generate_dataset.py

# Automated pipeline for generating HemeLB CFD ground truth data from AneuRisk STL meshes.

# Steps per patient:
#   1. Auto-generate .pr2 boundary profile from STL topology
#   2. Voxelize STL to .gmy lattice geometry via hlb-gmy-cli
#   3. Patch HemeLB XML (absolute .gmy path, step count, properties extraction block)
#   4. Purge stale output directory to prevent HemeLB file collision errors
#   5. Run HemeLB simulation via mpirun
#   6. Extract whole.xtr to CSV via hlb-dump-extracted-properties
#   7. Convergence check: relative velocity L2 change between last two timesteps < 1%
#   8. Mass conservation check: iolet flux imbalance < 2%

# Cases failing checks 7 or 8 are logged to data/excluded_cases.csv and excluded
# from training. Their CSV output is preserved for inspection.

# Boundary conditions (fixed throughout project, do not change):
#   Inlet  pressure: 1.0 mmHg
#   Outlet pressure: 0.0 mmHg
#   VoxelSize: 0.1 mm | TimeStepSeconds: 1e-5 s | Steps: 1000 (prototype)
# """

# import os
# import re
# import csv
# import subprocess
# import logging
# import traceback
# import shutil
# import xml.etree.ElementTree as ET
# from datetime import datetime
# import numpy as np
# import trimesh
# import networkx as nx
# from pathlib import Path

# # ── Configuration ─────────────────────────────────────────────────────────────
# HEMELB_BIN     = "hemelb"
# SETUP_TOOL_BIN = "hlb-gmy-cli"
# HLB_DUMP_BIN   = "hlb-dump-extracted-properties"
# MPI_CORES      = 4

# BASE_DIR = Path(__file__).parent
# RAW_DIR  = BASE_DIR / "data/raw_meshes"
# GMY_DIR  = BASE_DIR / "data/processed_gmy"
# OUT_DIR  = BASE_DIR / "data/outputs_csv"

# # Quality check thresholds
# CONVERGENCE_THRESHOLD       = 0.01   # Max relative velocity L2 change between last two timesteps
# MASS_CONSERVATION_THRESHOLD = 0.02   # Max relative flux imbalance across all iolets

# # Exclusion log path
# EXCLUSION_LOG = BASE_DIR / "data" / "excluded_cases.csv"

# # Regex for extracting floats/ints from text lines
# NUM_PATTERN = re.compile(r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?')

# logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')


# # ── Custom exception ──────────────────────────────────────────────────────────

# class SimulationQualityError(Exception):
#     """Raised when a simulation fails a quality check. Patient is excluded from training."""
#     pass


# # ── PR2 generation (unchanged) ────────────────────────────────────────────────

# def auto_generate_pr2(mesh_path, pr2_path):
#     """Parses an open STL mesh, extracts boundaries, filters topological noise,
#     and generates a valid .pr2 profile."""
#     logging.info(f"Withdrawing mesh topology markers for: {mesh_path.name}")

#     mesh = trimesh.load(mesh_path)

#     edges = mesh.edges_sorted
#     unique_edges, counts = np.unique(edges, axis=0, return_counts=True)
#     boundary_edges = unique_edges[counts == 1]

#     if len(boundary_edges) == 0:
#         raise ValueError(
#             f"Mesh {mesh_path.name} has no open boundaries. "
#             f"HemeLB requires open holes at inlets and outlets."
#         )

#     g = nx.Graph()
#     g.add_edges_from(boundary_edges)
#     all_loops = list(nx.connected_components(g))

#     loops = [l for l in all_loops if len(l) >= 20]

#     if len(loops) < 2:
#         logging.warning("Topological filtering was too aggressive. Falling back to largest two components.")
#         all_loops.sort(key=len, reverse=True)
#         loops = all_loops[:2]

#     logging.info(f"Identified {len(loops)} genuine fluid boundaries after topological filtering.")

#     centroid = mesh.center_mass

#     pr2_content = []
#     pr2_content.append("DurationSeconds: 6.0")
#     pr2_content.append("Iolets:")

#     inlet_counter  = 1
#     outlet_counter = 1

#     loops.sort(key=len, reverse=True)

#     for idx, loop in enumerate(loops):
#         loop_nodes    = list(loop)
#         loop_vertices = mesh.vertices[loop_nodes]

#         center = np.mean(loop_vertices, axis=0)
#         radius = np.max(np.linalg.norm(loop_vertices - center, axis=1)) + 0.1

#         _, _, vh = np.linalg.svd(loop_vertices - center)
#         normal = vh[2, :]

#         if np.dot(normal, centroid - center) < 0:
#             normal = -normal

#         if idx == 0:
#             iolet_type = "Inlet"
#             iolet_name = f"Inlet{inlet_counter}"
#             inlet_counter += 1
#             # pressure_x = 1.0   # 1.0 mmHg inlet BC - do not change mid-project
#             pressure_x = 0.1   # 0.1 mmHg inlet BC - stable across all AneuRisk geometries
#         else:
#             iolet_type = "Outlet"
#             iolet_name = f"Outlet{outlet_counter}"
#             outlet_counter += 1
#             pressure_x = 0.0

#         pr2_content.append("- Centre:")
#         pr2_content.append(f"    x: {center[0]:.10f}")
#         pr2_content.append(f"    y: {center[1]:.10f}")
#         pr2_content.append(f"    z: {center[2]:.10f}")
#         pr2_content.append(f"  Name: {iolet_name}")
#         pr2_content.append("  Normal:")
#         pr2_content.append(f"    x: {normal[0]:.10f}")
#         pr2_content.append(f"    y: {normal[1]:.10f}")
#         pr2_content.append(f"    z: {normal[2]:.10f}")
#         pr2_content.append("  Pressure:")
#         pr2_content.append(f"    x: {pressure_x:.10f}")
#         pr2_content.append("    y: 0.0")
#         pr2_content.append("    z: 0.0")
#         pr2_content.append(f"  Radius: {radius:.10f}")
#         pr2_content.append(f"  Type: {iolet_type}")

#     gmy_out_name = f"{mesh_path.with_suffix('').name}.gmy"
#     xml_out_name = f"{mesh_path.with_suffix('').name}_input.xml"

#     pr2_content.append(f"OutputGeometryFile: {gmy_out_name}")
#     pr2_content.append(f"OutputXmlFile: {xml_out_name}")
#     pr2_content.append("SeedPoint:")
#     pr2_content.append(f"  x: {centroid[0]:.10f}")
#     pr2_content.append(f"  y: {centroid[1]:.10f}")
#     pr2_content.append(f"  z: {centroid[2]:.10f}")
#     pr2_content.append(f"StlFile: {mesh_path.name}")
#     pr2_content.append("StlFileUnitId: 1")
#     pr2_content.append("TimeStepSeconds: 0.00005")
#     pr2_content.append("VoxelSize: 0.1")

#     with open(pr2_path, "w") as f:
#         f.write("\n".join(pr2_content))
#     logging.info(f"Successfully auto-generated profile: {pr2_path.name}")


# # ── Shell execution (unchanged) ───────────────────────────────────────────────

# def run_cmd(cmd, step_name):
#     """Executes a shell command with full stdout/stderr capture on failure."""
#     logging.info(f"Starting: {step_name}")
#     try:
#         subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         logging.info(f"Success: {step_name}")
#     except subprocess.CalledProcessError as e:
#         stdout_msg = e.stdout.decode('utf-8') if e.stdout else "None"
#         stderr_msg = e.stderr.decode('utf-8') if e.stderr else "None"
#         logging.error(
#             f"FAILED: {step_name}\n\n"
#             f"[STANDARD ERROR]:\n{stderr_msg}\n\n"
#             f"[STANDARD OUTPUT]:\n{stdout_msg}"
#         )
#         raise


# # ── PR2 parser (new) ──────────────────────────────────────────────────────────

# def parse_pr2_iolets(pr2_path):
#     """
#     Parse iolet definitions from a .pr2 profile file.

#     Returns
#     -------
#     list of dict, each with keys:
#         'centre'  : [x, y, z] in mm  (STL coordinate units)
#         'normal'  : [nx, ny, nz]     (unit vector pointing into fluid domain)
#         'radius'  : float in mm
#         'type'    : 'Inlet' or 'Outlet'
#     """
#     TOP_LEVEL_EXIT_KEYS = {
#         'OutputGeometryFile', 'OutputXmlFile', 'SeedPoint',
#         'StlFile', 'StlFileUnitId', 'TimeStepSeconds', 'VoxelSize', 'DurationSeconds'
#     }

#     with open(pr2_path, 'r') as f:
#         lines = f.readlines()

#     iolets          = []
#     current         = None
#     current_sub     = None   # 'centre' | 'normal' | None
#     in_iolets_block = False

#     for line in lines:
#         stripped = line.strip()
#         if not stripped:
#             continue

#         # Enter the Iolets block
#         if stripped == 'Iolets:':
#             in_iolets_block = True
#             continue

#         if not in_iolets_block:
#             continue

#         # Detect top-level keys that signal the end of the Iolets block
#         key_candidate = stripped.split(':')[0].strip()
#         if key_candidate in TOP_LEVEL_EXIT_KEYS:
#             if current is not None:
#                 iolets.append(current)
#                 current = None
#             in_iolets_block = False
#             continue

#         # New iolet entry starts with "- Centre:"
#         if stripped == '- Centre:':
#             if current is not None:
#                 iolets.append(current)
#             current     = {'centre': [0.0, 0.0, 0.0], 'normal': [0.0, 0.0, 0.0],
#                            'radius': 0.0, 'type': None}
#             current_sub = 'centre'
#             continue

#         if stripped == 'Normal:':
#             current_sub = 'normal'
#             continue

#         # Keys that end a coordinate sub-block
#         if stripped in ('Pressure:', 'Centre:') or stripped.startswith('Name:'):
#             current_sub = None
#             continue

#         if stripped.startswith('Radius:') and current is not None:
#             current['radius'] = float(stripped.split(':', 1)[1].strip())
#             current_sub = None
#             continue

#         if stripped.startswith('Type:') and current is not None:
#             current['type'] = stripped.split(':', 1)[1].strip()
#             current_sub = None
#             continue

#         # Parse x/y/z coordinates for current sub-block
#         if current_sub in ('centre', 'normal') and current is not None:
#             if stripped.startswith('x:'):
#                 current[current_sub][0] = float(stripped[2:].strip())
#             elif stripped.startswith('y:'):
#                 current[current_sub][1] = float(stripped[2:].strip())
#             elif stripped.startswith('z:'):
#                 current[current_sub][2] = float(stripped[2:].strip())

#     if current is not None:
#         iolets.append(current)

#     if not iolets:
#         raise ValueError(f"No iolets parsed from {pr2_path}. Check the file format.")

#     return iolets


# # ── CSV parsing utilities (new) ───────────────────────────────────────────────

# def _parse_csv_header(lines):
#     """
#     Extract geometry origin (m) and voxel size (m) from the CSV header comment block.
#     Returns (origin: np.ndarray shape (3,), voxel_size: float).
#     """
#     origin     = None
#     voxel_size = None

#     for line in lines:
#         if line.startswith('# Geometry origin'):
#             nums = NUM_PATTERN.findall(line)
#             if len(nums) < 3:
#                 raise ValueError(f"Cannot parse origin from: {line.strip()}")
#             origin = np.array([float(n) for n in nums[:3]])
#         elif line.startswith('# Voxel size'):
#             nums = NUM_PATTERN.findall(line)
#             if not nums:
#                 raise ValueError(f"Cannot parse voxel size from: {line.strip()}")
#             voxel_size = float(nums[0])

#     if origin is None:
#         raise ValueError("'# Geometry origin' not found in CSV header.")
#     if voxel_size is None:
#         raise ValueError("'# Voxel size' not found in CSV header.")

#     return origin, voxel_size


# def _parse_data_block(lines, start_idx):
#     """
#     Parse data rows from lines[start_idx+1:] until the next comment line or EOF.
#     Returns (grid: np.ndarray (N,3) int32, velocity: np.ndarray (N,3) float64).
#     """
#     grid_list = []
#     vel_list  = []

#     for line in lines[start_idx + 1:]:
#         stripped = line.strip()
#         if not stripped or stripped.startswith('#'):
#             break
#         nums = [float(x) for x in NUM_PATTERN.findall(stripped)]
#         if len(nums) != 7:
#             continue
#         grid_list.append(nums[0:3])
#         vel_list.append(nums[3:6])

#     if not grid_list:
#         raise ValueError(f"No valid data rows found after line {start_idx}.")

#     return (np.array(grid_list, dtype=np.int32),
#             np.array(vel_list,  dtype=np.float64))


# def parse_final_two_timesteps(csv_path):
#     """
#     Parse the last two timestep blocks from a HemeLB extracted properties CSV.

#     Returns
#     -------
#     origin      : np.ndarray (3,)  in metres
#     voxel_size  : float            in metres
#     grid_a      : np.ndarray (N,3) second-to-last timestep grid indices
#     vel_a       : np.ndarray (N,3) second-to-last timestep velocities (m/s)
#     grid_b      : np.ndarray (N,3) final timestep grid indices
#     vel_b       : np.ndarray (N,3) final timestep velocities (m/s)

#     Raises ValueError if fewer than two timestep blocks are found.
#     """
#     with open(csv_path, 'r') as f:
#         lines = f.readlines()

#     origin, voxel_size = _parse_csv_header(lines)

#     ts_line_indices = [i for i, l in enumerate(lines) if l.startswith('# Timestep')]

#     if len(ts_line_indices) < 2:
#         raise ValueError(
#             f"Need at least 2 timestep blocks for convergence check; "
#             f"found {len(ts_line_indices)} in {csv_path}."
#         )

#     idx_a = ts_line_indices[-2]   # Second-to-last timestep
#     idx_b = ts_line_indices[-1]   # Final timestep

#     logging.info(f"  Convergence: comparing {lines[idx_a].strip()} vs {lines[idx_b].strip()}")

#     grid_a, vel_a = _parse_data_block(lines, idx_a)
#     grid_b, vel_b = _parse_data_block(lines, idx_b)

#     return origin, voxel_size, grid_a, vel_a, grid_b, vel_b


# def parse_final_timestep(csv_path):
#     """
#     Parse only the final timestep block. Used by the mass conservation check.

#     Returns
#     -------
#     origin     : np.ndarray (3,)  in metres
#     voxel_size : float            in metres
#     grid       : np.ndarray (N,3) int32
#     velocity   : np.ndarray (N,3) float64 in m/s
#     """
#     with open(csv_path, 'r') as f:
#         lines = f.readlines()

#     origin, voxel_size = _parse_csv_header(lines)

#     ts_line_indices = [i for i, l in enumerate(lines) if l.startswith('# Timestep')]

#     if not ts_line_indices:
#         raise ValueError(f"No timestep blocks found in {csv_path}.")

#     # Collect all non-comment, non-empty lines after the last timestep marker
#     last_idx  = ts_line_indices[-1]
#     data_lines_raw = [
#         line for line in lines[last_idx + 1:]
#         if line.strip() and not line.strip().startswith('#')
#     ]

#     grid_list = []
#     vel_list  = []

#     for line in data_lines_raw:
#         nums = [float(x) for x in NUM_PATTERN.findall(line.strip())]
#         if len(nums) != 7:
#             continue
#         grid_list.append(nums[0:3])
#         vel_list.append(nums[3:6])

#     return (origin, voxel_size,
#             np.array(grid_list, dtype=np.int32),
#             np.array(vel_list,  dtype=np.float64))


# # ── Quality checks (new) ──────────────────────────────────────────────────────

# def check_convergence(csv_path, threshold=CONVERGENCE_THRESHOLD):
#     """
#     Compare the velocity field between the last two extraction timesteps.

#     Metric: relative L2 norm of the velocity difference.
#         rel_change = ||V_last - V_second_last||_F / ||V_second_last||_F

#     Sites are matched by grid index; both blocks must cover the same set of sites.

#     Parameters
#     ----------
#     csv_path  : Path to the HemeLB fluid data CSV.
#     threshold : Maximum allowed relative change (default 0.01 = 1%).

#     Returns
#     -------
#     passed     : bool
#     rel_change : float
#     """
#     _, _, grid_a, vel_a, grid_b, vel_b = parse_final_two_timesteps(csv_path)

#     # Match sites by grid index (sort both arrays by grid tuple)
#     def sort_key(grid):
#         return np.lexsort((grid[:, 2], grid[:, 1], grid[:, 0]))

#     order_a = sort_key(grid_a)
#     order_b = sort_key(grid_b)

#     grid_a_sorted = grid_a[order_a]
#     grid_b_sorted = grid_b[order_b]

#     if not np.array_equal(grid_a_sorted, grid_b_sorted):
#         raise ValueError(
#             "Grid indices differ between the last two timestep blocks. "
#             "Cannot compute convergence metric."
#         )

#     vel_a_sorted = vel_a[order_a]
#     vel_b_sorted = vel_b[order_b]

#     diff_norm    = np.linalg.norm(vel_b_sorted - vel_a_sorted)
#     ref_norm     = np.linalg.norm(vel_a_sorted)

#     if ref_norm < 1e-30:
#         logging.warning("  Convergence check: reference velocity norm is effectively zero. Marking as failed.")
#         return False, float('inf')

#     rel_change = diff_norm / ref_norm
#     passed     = rel_change < threshold

#     logging.info(
#         f"  Convergence check: rel_change = {rel_change:.4e} "
#         f"(threshold {threshold:.2e}) -> {'PASS' if passed else 'FAIL'}"
#     )

#     return passed, rel_change


# def check_mass_conservation(csv_path, pr2_path, threshold=MASS_CONSERVATION_THRESHOLD):
#     """
#     Verify that total inflow flux equals total outflow flux at the iolet planes.

#     Method:
#     - For each iolet, identify fluid nodes within 1.5 voxel widths of the iolet plane
#       and within the iolet radius of the iolet centre.
#     - Compute signed flux: F_i = sum(v . n_hat_i) * voxel_size^2
#       where n_hat_i is the iolet inward normal (pointing into fluid domain).
#     - For an incompressible fluid the net flux across all boundaries must be zero:
#         sum_i(F_i) = 0
#     - Relative imbalance: |sum_i(F_i)| / |F_inlet|

#     Notes on coordinate systems:
#     - pr2 centre/radius are in mm (STL file units, StlFileUnitId = 1).
#     - CSV physical coordinates are in metres.
#     - Conversion: divide pr2 values by 1000.

#     Parameters
#     ----------
#     csv_path  : Path to the HemeLB fluid data CSV.
#     pr2_path  : Path to the patient .pr2 profile file.
#     threshold : Maximum allowed relative flux imbalance (default 0.02 = 2%).

#     Returns
#     -------
#     passed   : bool
#     rel_error: float
#     """
#     iolets                     = parse_pr2_iolets(pr2_path)
#     origin, voxel_size, grid, velocity = parse_final_timestep(csv_path)

#     # Physical coordinates of all fluid nodes (metres)
#     xyz = origin + grid.astype(np.float64) * voxel_size   # (N, 3)

#     # Search radius: 1.5 voxel widths normal to the iolet plane
#     plane_dist_threshold = 1.5 * voxel_size

#     inlet_flux  = None
#     outlet_flux = 0.0
#     fluxes      = []

#     for iolet in iolets:
#         centre_m = np.array(iolet['centre']) / 1000.0   # mm -> m
#         radius_m = iolet['radius']            / 1000.0  # mm -> m
#         normal   = np.array(iolet['normal'],  dtype=np.float64)
#         n_hat    = normal / np.linalg.norm(normal)

#         # Signed distance from each node to the iolet plane
#         disp          = xyz - centre_m[np.newaxis, :]           # (N, 3)
#         plane_dist    = np.abs(disp @ n_hat)                    # (N,)

#         # In-plane distance from iolet centre
#         proj          = (disp @ n_hat)[:, np.newaxis] * n_hat   # (N, 3)
#         in_plane_dist = np.linalg.norm(disp - proj, axis=1)     # (N,)

#         mask = (plane_dist < plane_dist_threshold) & (in_plane_dist < radius_m)
#         n_nodes = mask.sum()

#         if n_nodes == 0:
#             logging.warning(
#                 f"  Mass conservation: no fluid nodes found near iolet "
#                 f"'{iolet['type']}' (centre={centre_m}, radius={radius_m:.4f} m). "
#                 f"Skipping this iolet."
#             )
#             continue

#         # Flux = sum(v . n_hat) * dx^2
#         # Positive = flow in the inward-normal direction (into fluid)
#         flux = float(np.sum(velocity[mask] @ n_hat)) * voxel_size ** 2
#         fluxes.append(flux)

#         logging.info(
#             f"  Mass conservation: {iolet['type']} | "
#             f"nodes={n_nodes} | flux={flux:.4e} m^3/s"
#         )

#         if iolet['type'] == 'Inlet':
#             inlet_flux = flux
#         else:
#             outlet_flux += flux

#     if inlet_flux is None:
#         logging.warning("  Mass conservation: inlet flux not computed. Cannot perform check.")
#         return True, 0.0   # Cannot check; do not exclude

#     if abs(inlet_flux) < 1e-30:
#         logging.warning("  Mass conservation: inlet flux is effectively zero. Marking as failed.")
#         return False, float('inf')

#     net_flux  = inlet_flux + outlet_flux
#     rel_error = abs(net_flux) / abs(inlet_flux)
#     passed    = rel_error < threshold

#     logging.info(
#         f"  Mass conservation: net_flux={net_flux:.4e} m^3/s | "
#         f"rel_error={rel_error:.4e} (threshold {threshold:.2e}) "
#         f"-> {'PASS' if passed else 'FAIL'}"
#     )

#     return passed, rel_error


# # ── Exclusion logger (new) ────────────────────────────────────────────────────

# def log_exclusion(patient_id, reason, metric_name, metric_value):
#     """
#     Append an exclusion record to the exclusion log CSV.

#     Columns: timestamp, patient_id, reason, metric_name, metric_value
#     Creates the file and writes a header row if it does not exist.
#     """
#     EXCLUSION_LOG.parent.mkdir(parents=True, exist_ok=True)
#     write_header = not EXCLUSION_LOG.exists()

#     with open(EXCLUSION_LOG, 'a', newline='') as f:
#         writer = csv.writer(f)
#         if write_header:
#             writer.writerow(['timestamp', 'patient_id', 'reason', 'metric_name', 'metric_value'])
#         writer.writerow([
#             datetime.now().isoformat(timespec='seconds'),
#             patient_id,
#             reason,
#             metric_name,
#             f"{metric_value:.6e}"
#         ])

#     logging.warning(
#         f"  EXCLUDED: {patient_id} | {reason} | {metric_name}={metric_value:.4e} "
#         f"| Logged to {EXCLUSION_LOG.name}"
#     )


# # ── Patient processing (steps 1-6 unchanged, 7-8 added) ──────────────────────

# def process_patient(mesh_file):
#     patient_id = mesh_file.stem
#     logging.info(f"============ Processing Case: {patient_id} ============")

#     patient_pr2     = RAW_DIR / f"{patient_id}.pr2"
#     patient_gmy     = GMY_DIR / f"{patient_id}.gmy"
#     patient_xml     = GMY_DIR / f"{patient_id}_input.xml"
#     patient_out_dir = OUT_DIR / patient_id

#     xml_out_name = f"{mesh_file.stem}_input.xml"
#     gmy_out_name = f"{mesh_file.stem}.gmy"

#     # STEP 1: Topological feature extraction and PR2 generation
#     auto_generate_pr2(mesh_file, patient_pr2)

#     # STEP 2: Voxelization
#     setup_cmd = f"cd {RAW_DIR} && {SETUP_TOOL_BIN} {patient_pr2.name}"
#     run_cmd(setup_cmd, f"Headless Voxelization ({patient_id})")
#     os.rename(RAW_DIR / gmy_out_name, patient_gmy)

#     # STEP 3: Patch XML
#     generated_xml_path = RAW_DIR / xml_out_name
#     logging.info("Patching HemeLB XML with absolute .gmy path and properties extraction block.")

#     tree = ET.parse(generated_xml_path)
#     root = tree.getroot()

#     datafile = root.find(".//geometry/datafile")
#     if datafile is not None:
#         datafile.set("path", str(patient_gmy.resolve()))

#     steps_element = root.find(".//simulation/steps")
#     if steps_element is not None:
#         # steps_element.set("value", "1000")
#         steps_element.set("value", "5000")

#     properties = root.find("properties")
#     if properties is None:
#         properties  = ET.SubElement(root, "properties")
#         # prop_output = ET.SubElement(properties, "propertyoutput", {"file": "whole.xtr", "period": "100"})
#         prop_output = ET.SubElement(properties, "propertyoutput", {"file": "whole.xtr", "period": "500"})
#         ET.SubElement(prop_output, "geometry", {"type": "whole"})
#         ET.SubElement(prop_output, "field",    {"type": "velocity"})
#         ET.SubElement(prop_output, "field",    {"type": "pressure"})

#     tree.write(patient_xml, encoding="utf-8", xml_declaration=True)

#     if generated_xml_path.exists():
#         os.remove(generated_xml_path)

#     # STEP 4: Purge stale output directory
#     if patient_out_dir.exists():
#         logging.info(f"Purging pre-existing output directory: {patient_out_dir.name}")
#         shutil.rmtree(patient_out_dir)

#     # DIAGNOSTIC: Print XML before running
#     try:
#         with open(patient_xml, 'r') as f:
#             logging.info(
#                 f"\n{'='*50}\n--- DIAGNOSTIC: GENERATED XML ---\n{'='*50}\n"
#                 f"{f.read()}\n{'='*50}"
#             )
#     except Exception as e:
#         logging.error(f"Could not read XML for diagnostic: {e}")

#     # STEP 5: Run HemeLB simulation
#     hemelb_cmd = f"mpirun -n {MPI_CORES} {HEMELB_BIN} -in {patient_xml} -out {patient_out_dir}"
#     run_cmd(hemelb_cmd, f"HemeLB Simulation ({patient_id})")

#     # STEP 6: Extract whole.xtr to CSV
#     xtr_matches = list(patient_out_dir.glob("**/whole.xtr"))

#     if not xtr_matches:
#         existing = [str(p.relative_to(BASE_DIR)) for p in patient_out_dir.rglob("*") if p.is_file()]
#         raise FileNotFoundError(
#             f"whole.xtr not found in {patient_out_dir}.\n"
#             f"Files present:\n" + "\n".join(existing)
#         )

#     xtr_path = xtr_matches[0]
#     csv_path = patient_out_dir / f"{patient_id}_fluid_data.csv"

#     extract_cmd = f"{HLB_DUMP_BIN} {xtr_path} > {csv_path}"
#     run_cmd(extract_cmd, f"CSV Extraction ({patient_id})")

#     # STEP 7: Convergence check
#     logging.info(f"--- Quality Check: Convergence ({patient_id}) ---")
#     try:
#         conv_passed, rel_change = check_convergence(csv_path)
#     except Exception as e:
#         logging.error(f"  Convergence check raised an exception: {e}")
#         log_exclusion(patient_id, "convergence_check_error", "exception", float('nan'))
#         raise SimulationQualityError(f"Convergence check failed with exception for {patient_id}.") from e

#     if not conv_passed:
#         log_exclusion(patient_id, "not_converged", "rel_velocity_change", rel_change)
#         raise SimulationQualityError(
#             f"{patient_id} excluded: velocity not converged "
#             f"(rel_change={rel_change:.4e} > threshold={CONVERGENCE_THRESHOLD})."
#         )

#     # STEP 8: Mass conservation check
#     logging.info(f"--- Quality Check: Mass Conservation ({patient_id}) ---")
#     try:
#         mass_passed, rel_error = check_mass_conservation(csv_path, patient_pr2)
#     except Exception as e:
#         logging.error(f"  Mass conservation check raised an exception: {e}")
#         log_exclusion(patient_id, "mass_conservation_check_error", "exception", float('nan'))
#         raise SimulationQualityError(f"Mass conservation check failed with exception for {patient_id}.") from e

#     if not mass_passed:
#         log_exclusion(patient_id, "mass_not_conserved", "rel_flux_imbalance", rel_error)
#         raise SimulationQualityError(
#             f"{patient_id} excluded: mass not conserved "
#             f"(rel_error={rel_error:.4e} > threshold={MASS_CONSERVATION_THRESHOLD})."
#         )

#     logging.info(f"============ Case {patient_id}: All checks PASSED ============\n")


# # ── Entry point ───────────────────────────────────────────────────────────────

# if __name__ == "__main__":
#     GMY_DIR.mkdir(parents=True, exist_ok=True)
#     OUT_DIR.mkdir(parents=True, exist_ok=True)

#     mesh_files = sorted(list(RAW_DIR.glob("*.stl")) + list(RAW_DIR.glob("*.vtp")))

#     if not mesh_files:
#         logging.warning(f"No surface files found in {RAW_DIR}.")

#     passed_cases  = []
#     excluded_cases = []
#     failed_cases   = []

#     for mesh in mesh_files:
#         try:
#             process_patient(mesh)
#             passed_cases.append(mesh.stem)

#         except SimulationQualityError as e:
#             # Expected exclusion: quality check failed. Already logged to excluded_cases.csv.
#             logging.warning(f"Quality exclusion: {e}")
#             excluded_cases.append(mesh.stem)

#         except Exception:
#             # Unexpected pipeline error (voxelization crash, HemeLB crash, etc.)
#             logging.error(
#                 f"Pipeline error for {mesh.name}:\n{traceback.format_exc()}"
#             )
#             failed_cases.append(mesh.stem)

#     # Final summary
#     logging.info(
#         f"\n{'='*60}\n"
#         f"  BATCH COMPLETE\n"
#         f"  Passed    : {len(passed_cases)}\n"
#         f"  Excluded  : {len(excluded_cases)}  (see {EXCLUSION_LOG.name})\n"
#         f"  Errors    : {len(failed_cases)}\n"
#         f"  Total     : {len(mesh_files)}\n"
#         f"{'='*60}"
#     )

#     if excluded_cases:
#         logging.info(f"  Excluded cases : {', '.join(excluded_cases)}")
#     if failed_cases:
#         logging.info(f"  Error cases    : {', '.join(failed_cases)}")


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
  Inlet  pressure: 0.1 mmHg
  Outlet pressure: 0.0 mmHg
  VoxelSize: 0.1 mm | TimeStepSeconds: 1e-5 s | Steps: 5000 (benchmark)

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
HEMELB_BIN     = "hemelb"
SETUP_TOOL_BIN = "hlb-gmy-cli"
HLB_DUMP_BIN   = "hlb-dump-extracted-properties"
MPI_CORES      = 4

BASE_DIR = Path(__file__).parent
RAW_DIR  = BASE_DIR / "data/raw_meshes"
GMY_DIR  = BASE_DIR / "data/processed_gmy"
OUT_DIR    = BASE_DIR / "data/outputs_csv"
MODELS_DIR = BASE_DIR / "AneuriskDatabase" / "models"

# Quality check thresholds
CONVERGENCE_THRESHOLD       = 0.05   # Max relative velocity L2 change between last two timesteps
MASS_CONSERVATION_THRESHOLD = 0.02   # Max relative flux imbalance across all iolets

# Exclusion log path
EXCLUSION_LOG = BASE_DIR / "data" / "excluded_cases.csv"

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

    Boundary conditions (fixed for the entire project):
      Inlet  pressure: 0.1 mmHg
      Outlet pressure: 0.0 mmHg
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
    """Executes a shell command with full stdout/stderr capture on failure."""
    logging.info(f"Starting: {step_name}")
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.info(f"Success: {step_name}")
    except subprocess.CalledProcessError as e:
        stdout_msg = e.stdout.decode('utf-8') if e.stdout else "None"
        stderr_msg = e.stderr.decode('utf-8') if e.stderr else "None"
        logging.error(
            f"FAILED: {step_name}\n\n"
            f"[STANDARD ERROR]:\n{stderr_msg}\n\n"
            f"[STANDARD OUTPUT]:\n{stdout_msg}"
        )
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
    plane_dist_threshold = 3 * voxel_size

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

def process_patient(mesh_file):
    patient_id = mesh_file.stem
    logging.info(f"============ Processing Case: {patient_id} ============")

    patient_pr2     = RAW_DIR / f"{patient_id}.pr2"
    patient_gmy     = GMY_DIR / f"{patient_id}.gmy"
    patient_xml     = GMY_DIR / f"{patient_id}_input.xml"
    patient_out_dir = OUT_DIR / patient_id

    xml_out_name = f"{mesh_file.stem}_input.xml"
    gmy_out_name = f"{mesh_file.stem}.gmy"

    # STEP 1: Topological feature extraction and PR2 generation
    auto_generate_pr2(mesh_file, patient_pr2)

    # STEP 2: Voxelization
    setup_cmd = f"cd {RAW_DIR} && {SETUP_TOOL_BIN} {patient_pr2.name}"
    run_cmd(setup_cmd, f"Headless Voxelization ({patient_id})")
    os.rename(RAW_DIR / gmy_out_name, patient_gmy)

    # STEP 3: Patch XML
    generated_xml_path = RAW_DIR / xml_out_name
    logging.info("Patching HemeLB XML with absolute .gmy path and properties extraction block.")

    tree = ET.parse(generated_xml_path)
    root = tree.getroot()

    datafile = root.find(".//geometry/datafile")
    if datafile is not None:
        datafile.set("path", str(patient_gmy.resolve()))

    steps_element = root.find(".//simulation/steps")
    if steps_element is not None:
        steps_element.set("value", "20000")

    properties = root.find("properties")
    if properties is None:
        properties  = ET.SubElement(root, "properties")
        prop_output = ET.SubElement(properties, "propertyoutput", {"file": "whole.xtr", "period": "2000"})
        ET.SubElement(prop_output, "geometry", {"type": "whole"})
        ET.SubElement(prop_output, "field",    {"type": "velocity"})
        ET.SubElement(prop_output, "field",    {"type": "pressure"})

    tree.write(patient_xml, encoding="utf-8", xml_declaration=True)

    if generated_xml_path.exists():
        os.remove(generated_xml_path)

    # STEP 4: Purge stale output directory
    if patient_out_dir.exists():
        logging.info(f"Purging pre-existing output directory: {patient_out_dir.name}")
        shutil.rmtree(patient_out_dir)

    # DIAGNOSTIC: Print XML before running
    try:
        with open(patient_xml, 'r') as f:
            logging.info(
                f"\n{'='*50}\n--- DIAGNOSTIC: GENERATED XML ---\n{'='*50}\n"
                f"{f.read()}\n{'='*50}"
            )
    except Exception as e:
        logging.error(f"Could not read XML for diagnostic: {e}")

    # STEP 5: Run HemeLB simulation
    hemelb_cmd = f"mpirun -n {MPI_CORES} {HEMELB_BIN} -in {patient_xml} -out {patient_out_dir}"
    run_cmd(hemelb_cmd, f"HemeLB Simulation ({patient_id})")

    # STEP 6: Extract whole.xtr to CSV
    xtr_matches = list(patient_out_dir.glob("**/whole.xtr"))

    if not xtr_matches:
        existing = [str(p.relative_to(BASE_DIR)) for p in patient_out_dir.rglob("*") if p.is_file()]
        raise FileNotFoundError(
            f"whole.xtr not found in {patient_out_dir}.\n"
            f"Files present:\n" + "\n".join(existing)
        )

    xtr_path = xtr_matches[0]
    csv_path = patient_out_dir / f"{patient_id}_fluid_data.csv"

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

    passed_cases  = []
    excluded_cases = []
    failed_cases   = []

    for mesh in mesh_files:
        try:
            process_patient(mesh)
            passed_cases.append(mesh.stem)

        except SimulationQualityError as e:
            # Expected exclusion: quality check failed. Already logged to excluded_cases.csv.
            logging.warning(f"Quality exclusion: {e}")
            excluded_cases.append(mesh.stem)

        except Exception:
            # Unexpected pipeline error (voxelization crash, HemeLB crash, etc.)
            logging.error(
                f"Pipeline error for {mesh.name}:\n{traceback.format_exc()}"
            )
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
