#!/usr/bin/env python3
"""
build_dataset_manifest.py

Aggregates per-patient manifest.csv files from the AneuriskDatabase into a single
master manifest at data/dataset_manifest.csv.

Also cross-references the STL files present in data/raw_meshes/ to identify:
  - Cases with STL but no manifest (metadata missing)
  - Cases with manifest but no STL (not yet downloaded or failed)

Output columns (in addition to all manifest.csv fields):
  patient_id         : e.g. C0001
  has_stl            : bool
  has_manifest       : bool
  has_centerlines    : bool
  stl_path           : absolute path or empty
  centerlines_csv    : absolute path or empty

Usage
-----
    python build_dataset_manifest.py
"""

import csv
import logging
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - [%(levelname)s] - %(message)s'
)

BASE_DIR     = Path(__file__).parent
MODELS_DIR   = BASE_DIR / "AneuriskDatabase" / "models"
RAW_DIR      = BASE_DIR / "data" / "raw_meshes"
DATA_DIR     = BASE_DIR / "data"
OUTPUT_PATH  = DATA_DIR / "dataset_manifest.csv"


def build_manifest():
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    # All patient directories in the repo
    patient_dirs = sorted([d for d in MODELS_DIR.iterdir() if d.is_dir()])
    # All STL files present in raw_meshes
    stl_stems = {p.stem for p in RAW_DIR.glob("*.stl")}

    all_rows        = []
    missing_manifest = []
    missing_stl      = []
    missing_centerlines = []

    # Collect all manifest column names first (they should be identical across cases
    # but we union them defensively)
    all_manifest_keys = []
    for patient_dir in patient_dirs:
        mf = patient_dir / "manifest.csv"
        if mf.exists():
            with open(mf, newline='') as f:
                reader = csv.DictReader(f)
                for key in (reader.fieldnames or []):
                    if key not in all_manifest_keys:
                        all_manifest_keys.append(key)
            break  # One file is enough to get the column schema

    # Extra columns we add
    extra_keys = ['patient_id', 'has_stl', 'has_manifest',
                  'has_centerlines', 'stl_path', 'centerlines_csv']

    fieldnames = extra_keys + all_manifest_keys

    for patient_dir in patient_dirs:
        patient_id   = patient_dir.name
        manifest_path = patient_dir / "manifest.csv"
        stl_path      = RAW_DIR / f"{patient_id}.stl"
        cl_csv_path   = patient_dir / "morphology" / "centerlines.csv"

        has_stl          = stl_path.exists()
        has_manifest     = manifest_path.exists()
        has_centerlines  = cl_csv_path.exists()

        if not has_manifest:
            missing_manifest.append(patient_id)
        if not has_stl:
            missing_stl.append(patient_id)
        if not has_centerlines:
            missing_centerlines.append(patient_id)

        # Read manifest fields if available
        manifest_data = {k: "" for k in all_manifest_keys}
        if has_manifest:
            with open(manifest_path, newline='') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    for key in all_manifest_keys:
                        manifest_data[key] = row.get(key, "")
                    break  # Each manifest has exactly one data row

        row = {
            'patient_id'      : patient_id,
            'has_stl'         : has_stl,
            'has_manifest'    : has_manifest,
            'has_centerlines' : has_centerlines,
            'stl_path'        : str(stl_path) if has_stl else "",
            'centerlines_csv' : str(cl_csv_path) if has_centerlines else "",
        }
        row.update(manifest_data)
        all_rows.append(row)

    # Write master manifest
    with open(OUTPUT_PATH, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)

    # Report
    logging.info(f"\n{'='*60}")
    logging.info(f"  Total patient directories : {len(patient_dirs)}")
    logging.info(f"  STL files in raw_meshes   : {len(stl_stems)}")
    logging.info(f"  Cases with manifest       : {len(patient_dirs) - len(missing_manifest)}")
    logging.info(f"  Cases with centerlines    : {len(patient_dirs) - len(missing_centerlines)}")
    logging.info(f"  Master manifest written   : {OUTPUT_PATH}")
    logging.info(f"{'='*60}")

    if missing_manifest:
        logging.warning(
            f"  {len(missing_manifest)} cases missing manifest.csv "
            f"(no rupture status or metadata):\n  " +
            ", ".join(missing_manifest)
        )
    if missing_stl:
        logging.warning(
            f"  {len(missing_stl)} cases missing STL in raw_meshes:\n  " +
            ", ".join(missing_stl)
        )
    if missing_centerlines:
        logging.warning(
            f"  {len(missing_centerlines)} cases missing centerlines.csv:\n  " +
            ", ".join(missing_centerlines)
        )

    return all_rows


if __name__ == "__main__":
    build_manifest()
