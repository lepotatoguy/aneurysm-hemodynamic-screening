import os
import subprocess
import logging
from pathlib import Path

# # ==========================================
# # CONFIGURATION - UPDATE THESE PATHS
# # ==========================================
# HEMELB_BIN = "/usr/local/bin/hemelb"
# SETUP_TOOL_BIN = "hemelb-setup-main" # or path to the setup tool python script
# HEMEXTRACT_BIN = "/path/to/hemelb/build/hemeXtract"
# MPI_CORES = 4 # Adjust based on your local machine
# import os
# import subprocess
# import logging
# from pathlib import Path

# ==========================================
# CONFIGURATION - UPDATE THESE PATHS
# ==========================================
HEMELB_BIN = "hemelb"             # Accessible if installed in binary path (cite: 313)
SETUP_TOOL_BIN = "hlb-gmy-cli"     # Using the documented CLI option (cite: 285)
HLB_DUMP_BIN = "hlb-dump-extracted-properties" # The correct extraction tool (cite: 344)
MPI_CORES = 4                     # Adjust based on your system (cite: 312)

# Directory Setup
BASE_DIR = Path(__file__).parent
RAW_DIR = BASE_DIR / "data/raw_meshes"
GMY_DIR = BASE_DIR / "data/processed_gmy"
OUT_DIR = BASE_DIR / "data/outputs_csv"   # Changed to hold clean ML-ready CSVs
TEMPLATE_XML = BASE_DIR / "templates/base_input.xml"

logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')

def run_cmd(cmd, step_name):
    """Executes shell commands and logs outputs cleanly."""
    logging.info(f"Starting: {step_name}")
    try:
        # Using shell execution to allow for the '>' output redirection operator
        result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.info(f"Success: {step_name}")
    except subprocess.CalledProcessError as e:
        logging.error(f"FAILED: {step_name}\nError output:\n{e.stderr.decode('utf-8')}")
        raise

def process_patient(mesh_file):
    patient_id = mesh_file.stem  # e.g., "C0099"
    logging.info(f"--- Processing Patient: {patient_id} ---")
    
    # Path construction
    patient_gmy = GMY_DIR / f"{patient_id}.gmy"
    patient_xml = GMY_DIR / f"{patient_id}_input.xml"
    patient_out_dir = OUT_DIR / patient_id
    
    # Ensure patient-specific output directory exists
    patient_out_dir.mkdir(parents=True, exist_ok=True)

    # Note: Ensure your input mesh is pre-capped (capping=True) before running (cite: 539)
    
    # STEP 1: Geometry & Config Generation via Setup Tool
    # Custom adjustments can be handled here or pre-generated via your .pr2 profile method (cite: 220, 535)
    setup_cmd = f"{SETUP_TOOL_BIN} -i {mesh_file} -o {patient_gmy} -v 0.05" 
    # If using pre-generated configurations directly from the GUI, you can bypass this step.

    # STEP 2: Inject configuration path dynamically into base template
    with open(TEMPLATE_XML, 'r') as file:
        xml_content = file.read()
    
    xml_content = xml_content.format(gmy_path=str(patient_gmy))
    with open(patient_xml, 'w') as file:
        file.write(xml_content)
    logging.info(f"Generated runtime XML configuration for {patient_id}")

    # STEP 3: Run HemeLB Simulation via MPI
    # Documented format: mpirun -n N <exec> in <xml> -out <dir> (cite: 289)
    hemelb_cmd = f"mpirun -n {MPI_CORES} {HEMELB_BIN} in {patient_xml} -out {patient_out_dir}"
    run_cmd(hemelb_cmd, f"HemeLB Simulation ({patient_id})")

    # STEP 4: Extract Information using hlb-dump-extracted-properties
    # Based on your documentation, HemeLB creates 'whole.xtr' inside the specified output directory (cite: 302, 312)
    xtr_path = patient_out_dir / "whole.xtr"
    csv_path = patient_out_dir / f"{patient_id}_fluid_data.csv"
    
    # Command format: hlb-dump-extracted-properties input.xtr > output.csv (cite: 344)
    extract_cmd = f"{HLB_DUMP_BIN} {xtr_path} > {csv_path}"
    run_cmd(extract_cmd, f"Extracting XTR to CSV ({patient_id})")
    
    logging.info(f"--- Finished {patient_id}. Dataset ready at: {csv_path} ---\n")

if __name__ == "__main__":
    GMY_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Automatically detects both .stl and .vtp files from your dataset (cite: 357)
    mesh_files = list(RAW_DIR.glob("*.stl")) + list(RAW_DIR.glob("*.vtp"))
    
    if not mesh_files:
        logging.warning(f"No vascular surface meshes found in {RAW_DIR}. Please place your unzipped case files there.")
    
    for mesh in mesh_files:
        try:
            process_patient(mesh)
        except Exception as e:
            logging.error(f"Skipping case {mesh.name} due to execution errors.")
            continue
            
    logging.info("Batch processing complete. Data generated seamlessly.")