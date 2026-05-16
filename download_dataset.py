#!/usr/bin/env python3
"""
download_aneurisk.py

Clones the AneuriskDatabase GitHub repository and copies all patient STL files
into the pipeline's raw_meshes directory, named by patient ID (e.g. C0001.stl).

Usage
-----
    python download_aneurisk.py

What it does
------------
1. Git-clones https://github.com/hkjeldsberg/AneuriskDatabase into a local
   'AneuriskDatabase' folder next to this script (skips clone if already present).
2. Walks the repo's 'models/' directory. Each subdirectory is a patient case
   (e.g. C0001, C0002, ...).
3. Inside each patient folder, looks for the STL file inside the 'surface/'
   subdirectory. Expected path pattern:
       AneuriskDatabase/models/<PATIENT_ID>/surface/<anything>.stl
4. Copies each found STL to:
       data/raw_meshes/<PATIENT_ID>.stl
5. Reports how many STL files were found, copied, and skipped.

Requirements
------------
Git must be available on PATH. No extra Python packages required.
"""

import shutil
import logging
import subprocess
import traceback
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - [%(levelname)s] - %(message)s'
)

# ── Configuration ─────────────────────────────────────────────────────────────

BASE_DIR   = Path(__file__).parent
REPO_URL   = "https://github.com/hkjeldsberg/AneuriskDatabase"
REPO_DIR   = BASE_DIR / "AneuriskDatabase"
MODELS_DIR = REPO_DIR / "models"
RAW_DIR    = BASE_DIR / "data" / "raw_meshes"


# ── Clone ─────────────────────────────────────────────────────────────────────

def clone_repo():
    if REPO_DIR.exists():
        logging.info(f"Repo already present at {REPO_DIR}. Skipping clone.")
        return

    logging.info(f"Cloning {REPO_URL} ...")
    result = subprocess.run(
        ["git", "clone", "--depth", "1", REPO_URL, str(REPO_DIR)],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"git clone failed:\n"
            f"STDOUT: {result.stdout}\n"
            f"STDERR: {result.stderr}"
        )
    logging.info("Clone complete.")


# ── STL discovery ─────────────────────────────────────────────────────────────

def collect_stl_files():
    """
    Walk models/ and return a list of (patient_id, stl_path) tuples.
    Expects exactly one STL per patient inside their surface/ subfolder.
    Warns if zero or more than one STL is found for a patient.
    """
    if not MODELS_DIR.exists():
        raise FileNotFoundError(
            f"models/ directory not found at {MODELS_DIR}. "
            f"Check that the repo cloned correctly."
        )

    patient_dirs = sorted([d for d in MODELS_DIR.iterdir() if d.is_dir()])
    logging.info(f"Found {len(patient_dirs)} patient directories in models/.")

    entries = []

    for patient_dir in patient_dirs:
        patient_id  = patient_dir.name        # e.g. C0001
        surface_dir = patient_dir / "surface"

        if not surface_dir.exists():
            logging.warning(f"  {patient_id}: no 'surface/' subfolder. Skipping.")
            continue

        # Case-insensitive glob for .stl
        stl_files = list(surface_dir.glob("*.stl")) + list(surface_dir.glob("*.STL"))

        if len(stl_files) == 0:
            logging.warning(f"  {patient_id}: no STL in surface/. Skipping.")
            continue

        if len(stl_files) > 1:
            logging.warning(
                f"  {patient_id}: {len(stl_files)} STL files found. "
                f"Using: {stl_files[0].name}"
            )

        entries.append((patient_id, stl_files[0]))

    return entries


# ── Copy ──────────────────────────────────────────────────────────────────────

def copy_stl_files(entries):
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    copied  = 0
    skipped = 0
    failed  = 0

    for patient_id, src_path in entries:
        dst_path = RAW_DIR / f"{patient_id}.stl"

        if dst_path.exists():
            logging.info(f"  {patient_id}: destination already exists. Skipping.")
            skipped += 1
            continue

        try:
            shutil.copy2(src_path, dst_path)
            logging.info(f"  {patient_id}: {src_path.name} -> {dst_path.name}")
            copied += 1
        except Exception:
            logging.error(f"  {patient_id}: copy failed.\n{traceback.format_exc()}")
            failed += 1

    return copied, skipped, failed


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    try:
        clone_repo()
    except Exception:
        logging.error(f"Repo clone failed:\n{traceback.format_exc()}")
        raise SystemExit(1)

    try:
        entries = collect_stl_files()
    except Exception:
        logging.error(f"STL collection failed:\n{traceback.format_exc()}")
        raise SystemExit(1)

    if not entries:
        logging.warning("No STL files found. Nothing written to data/raw_meshes/.")
        raise SystemExit(0)

    logging.info(f"Copying {len(entries)} STL files to {RAW_DIR} ...")
    copied, skipped, failed = copy_stl_files(entries)

    logging.info(
        f"\n{'='*50}\n"
        f"  Total patients found : {len(entries)}\n"
        f"  Copied               : {copied}\n"
        f"  Skipped (exist)      : {skipped}\n"
        f"  Failed               : {failed}\n"
        f"  Destination          : {RAW_DIR}\n"
        f"{'='*50}"
    )
