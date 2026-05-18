# HemeLB Surrogate Pipeline

A full pipeline for generating CFD-derived hemodynamic ground truth data from patient-specific cerebral aneurysm geometries (AneuRisk dataset) and training a Graph Neural Network (GNN) surrogate model to predict Wall Shear Stress (WSS) distributions for rupture risk screening.

---

## Project Status

| Step | Script | Status |
|------|--------|--------|
| 0. Download data | `download_aneurisk.py` | Done |
| 1. Build manifest | `build_dataset_manifest.py` | Done |
| 2. CFD simulation (generate dataset) | `generate_dataset.py` | In progress (batch running) |
| 3. WSS computation | `compute_wss.py` | Needs rewrite (STL-mapped output) |
| 4. Graph construction | - | Not started |
| 5. GNN training | - | Not started |
| 6. Clinical output | - | Not started |
| 7. Validation | - | Not started |

---

## Research Overview

**Problem:** Cerebral aneurysm rupture is lethal. Predicting rupture risk from patient geometry without running full CFD is the goal.

**Approach:** Use HemeLB (LBM-based CFD solver) to compute steady-state WSS across 95 AneuRisk patient geometries. Train a GNN surrogate that maps vessel geometry directly to normalised WSS distribution. Correlate predicted WSS with rupture status (available in AneuRisk manifest).

**Clinical framing:** A steady-state WSS surrogate as a computationally efficient geometry-to-hemodynamics screening tool, not a replacement for full clinical CFD.

---

## What This Project Does (Plain English)

**The medical problem:**
When a blood vessel in the brain develops a bulge (an aneurysm), it can burst and cause a stroke. Doctors want to know which aneurysms are dangerous and likely to burst, and which are stable and can be monitored. Right now there is no fast, reliable way to answer that question from a scan.

**What we know from research:**
The force that blood exerts on the vessel wall as it flows past is called Wall Shear Stress (WSS). Studies show that aneurysms with abnormally low WSS at their dome tend to rupture, while aneurysms with healthier WSS patterns tend to be stable. So if you can compute WSS for a patient, you have a clue about rupture risk.

**The problem with computing WSS:**
Simulating blood flow in a patient-specific vessel is expensive. It requires specialist software (we use HemeLB), a powerful computer, and hours of computation per patient. Most hospitals cannot do this routinely.

**What we are building:**
We are training a machine learning model that learns the relationship between the shape of a vessel and the WSS pattern it produces. Once trained, this model can look at a new patient's vessel geometry (from a brain scan) and predict the WSS distribution in seconds, without running a simulation at all. That is the surrogate model.

**How we build the training data:**
We take 95 real patient geometries from a public database (AneuRisk). For each one, we run a full blood flow simulation using HemeLB and record the WSS at every point on the vessel wall. This gives us 95 pairs of (geometry, WSS pattern), which we use to train the model.

**The machine learning model:**
A vessel wall is a 3D surface made of thousands of connected points. A Graph Neural Network (GNN) is a type of neural network designed to work on connected structures like this. Each point on the surface is a node in a graph, and the edges connect neighbouring points. The GNN learns how the local and global geometry at each node relates to the WSS at that node.

**The rupture risk connection:**
The AneuRisk database also tells us which of the 95 aneurysms actually ruptured. Once the GNN can predict WSS, we check whether the predicted WSS features (such as the fraction of the dome with very low WSS) are different between the ruptured and unruptured cases. If they are, the model can potentially be used as a fast screening tool for rupture risk.

**What this is not:**
This is not a clinical diagnostic tool. It is a proof-of-concept research project showing that a fast geometry-to-WSS surrogate is achievable. The simulations use simplified assumptions (steady flow, uniform blood viscosity, rigid walls) that a clinical tool would need to address before deployment.

**In one sentence:**
We simulate blood flow in 95 real patient aneurysms, use the results to train a machine learning model that predicts blood flow patterns from vessel shape alone, and check whether those patterns can distinguish dangerous aneurysms from safe ones.

---

## Dataset

**Source:** [AneuRisk Database](https://github.com/hkjeldsberg/AneuriskDatabase)

- 95 patient cases (C0001 to C0095, not all sequential)
- All cases have: STL surface mesh, VMTK centerlines, rupture status label
- Rupture status: R (ruptured) or U (unruptured) per `manifest.csv`
- Aneurysm locations: ICA, MCA, basilar, ACA (from `manifest.csv`)
- Aggregated metadata: `data/dataset_manifest.csv`

---

## Environment Setup

```bash
# Assumes conda is installed and HemeLB tools are available
conda activate gmy-tool

# Python dependencies
pip install trimesh networkx scipy numpy pandas

# Required external binaries (must be on PATH)
# hemelb        - compiled with LaddIolet (velocity BC support confirmed)
# hlb-gmy-cli   - HemeLB geometry setup tool
# hlb-dump-extracted-properties - XTR to CSV extractor
```

**Verify LaddIolet is compiled in:**
```bash
strings $(which hemelb) | grep -i "laddiolet"
# Should return: LADDIOLET, InOutLetParabolicVelocityE, etc.
```

---

## Simulation Parameters

All parameters are fixed across all 95 cases for dataset consistency.

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Voxel size (dx) | 0.2 mm | Within published range (0.2-0.27mm), fast enough for MacBook Air M5 |
| Time step (dt) | 2e-4 s | Gives tau=0.5499, stable BGK-LBM |
| Steps | 20,000 | 4 seconds physical time, ~4 convective timescales |
| Extraction period | 2000 | 10 timestep blocks per case |
| Inlet BC | Velocity parabolic | v_max=0.05 m/s, published HemeLB standard |
| Outlet BC | Pressure cosine | 0 mmHg, equal outlets (Murray's law self-similar tree) |
| Blood viscosity | 0.0035 Pa.s | Newtonian approximation |
| Blood density | 1060 kg/m³ | Standard |
| MPI cores | 4 | Stable on 16GB RAM MacBook Air |

**LBM stability checks:**
```
tau = 0.5 + 3 * nu * dt / dx^2 = 0.5499   (stable, well above 0.5)
c_s = dx / (sqrt(3) * dt) = 0.577 m/s
Ma  = v_max / c_s = 0.087                  (below 0.1 limit)
Re  = rho * v_mean * R / mu ~ 14           (at 0.025 m/s mean velocity)
```

**Non-physiological note:** Mean inlet velocity is 0.025 m/s vs physiological 0.3-0.5 m/s for ICA. WSS magnitudes are proportionally suppressed. The GNN trains on **normalised WSS** (nWSS = WSS / case_mean_WSS) which is pressure and velocity independent. Spatial distribution patterns are fully preserved.

---

## Inlet Detection

Inlet (parent vessel) is identified using VMTK centerlines from AneuRisk:

1. For each boundary loop, find nearest centerline point
2. The loop with the largest `MaximumInscribedSphereRadius` = parent vessel = inlet
3. Falls back to largest-loop heuristic if `centerlines.csv` is absent (none in current dataset)
4. Mismatch warning logged if centerline and geometric heuristic disagree

---

## Pipeline Scripts

### Step 0: Download AneuRisk Data

```bash
python download_aneurisk.py
```

- Shallow-clones `https://github.com/hkjeldsberg/AneuriskDatabase`
- Copies all `surface/model.stl` files to `data/raw_meshes/<PATIENT_ID>.stl`
- Idempotent: skips existing files and existing clone

### Step 1: Build Dataset Manifest

```bash
python build_dataset_manifest.py
```

- Aggregates all per-patient `manifest.csv` files into `data/dataset_manifest.csv`
- Reports missing STLs, missing manifests, missing centerlines
- Output columns: `patient_id, has_stl, has_manifest, has_centerlines, ruptureStatus, aneurysmLocation, age, sex, ...`

### Step 2: Generate CFD Dataset

```bash
python generate_dataset.py
```

Interactive prompt on startup:

```
[1] Full run      - rerun everything from scratch
[2] Auto-resume   - skip completed steps (default)
[3] Quality only  - rerun quality checks on existing CSVs
[4] Single case   - process one patient ID
```

**Pipeline steps per patient:**

| Step | Action |
|------|--------|
| 1 | Generate `.pr2` boundary profile (centerline-based inlet detection) |
| 2 | Voxelize STL to `.gmy` via `hlb-gmy-cli` |
| 3 | Patch XML: absolute .gmy path, 20000 steps, extraction block, **replace inlet with velocity parabolic BC** |
| 4 | Purge stale output directory |
| 5 | Run HemeLB via `mpirun -n 4` |
| 6 | Extract `whole.xtr` to CSV via `hlb-dump-extracted-properties` |
| 7 | Convergence check (rel velocity L2 change < 5% between steps 18000 and 20000) |
| 8 | Mass conservation logged (not used for exclusion - unreliable for pressure outlet BCs) |

**Outputs:**
- `data/processed_gmy/<ID>.gmy` - lattice geometry
- `data/processed_gmy/<ID>_input.xml` - HemeLB input XML with velocity BC
- `data/outputs_csv/<ID>/<ID>_fluid_data.csv` - velocity and pressure at all fluid nodes, 10 timesteps
- `data/excluded_cases.csv` - cases failing quality checks

**Safe stopping:** Press Ctrl+C between cases (after `All checks PASSED` line). On resume, option 2 auto-detects completed steps from existing files. Partial simulation outputs (output dir without CSV) are automatically deleted and rerun. Incomplete CSVs (fewer than 10 timestep blocks) are automatically deleted and rerun.

**To force-kill HemeLB:**
```bash
pkill -f hemelb; pkill -f mpirun
```

### Step 3: Compute WSS

```bash
python compute_wss.py
```

**Current version (fluid-lattice nodes):** Identifies wall-adjacent fluid nodes via 6-neighbor lookup, estimates wall normal from missing neighbors, computes WSS = mu * v_tangential / (0.5 * dx).

**Needs rewrite before Step 4:** Output must be per-STL-vertex (not per-fluid-node) for GNN graph construction. Planned: load STL, KDTree map STL vertices to nearest fluid nodes, use STL vertex normals, finite difference between two fluid nodes.

Output: `data/wss/<ID>_wss.csv` with columns `node_id, x, y, z, wss_magnitude`.

### Step 4: Verify Simulation (diagnostic)

```bash
python verify_simulation.py C0001   # or any patient ID
```

Prints 7 diagnostic checks: header metadata, velocity statistics (including top-10 nodes), pressure statistics, iolet velocity analysis with flux direction, Reynolds number estimate, spatial velocity distribution, WSS statistics.

---

## Quality Checks

### Convergence Check
Compares velocity field between timesteps 18000 and 20000:
```
rel_change = ||V_20000 - V_18000||_F / ||V_18000||_F < 0.05
```
Failure → case logged to `excluded_cases.csv` and excluded from training.

Confirmed passing cases:
- C0001: rel_change = 2.05e-06
- C0002: rel_change = 1.90e-05
- C0003: rel_change = 1.71e-06

### Mass Conservation Check
Logged only, never used for exclusion. HemeLB pressure outlet BCs do not guarantee flux balance at measurement planes, making this check unreliable. The convergence check is the sole quality gate.

---

## File Structure

```
HemeLB_Surrogate_Pipeline/
├── generate_dataset.py         Main CFD pipeline
├── compute_wss.py              WSS computation (needs rewrite for Step 4)
├── download_aneurisk.py        AneuRisk data download
├── build_dataset_manifest.py   Metadata aggregation
├── verify_simulation.py        Per-case diagnostic tool
├── master_prompt_final.md      Full research specification
│
├── AneuriskDatabase/           Git clone of AneuRisk repo
│   └── models/
│       └── C0001/
│           ├── surface/model.stl
│           ├── morphology/centerlines.csv
│           └── manifest.csv
│
└── data/
    ├── raw_meshes/             STL files + .pr2 profiles
    │   ├── C0001.stl
    │   ├── C0001.pr2
    │   └── ...
    ├── processed_gmy/          HemeLB geometry files
    │   ├── C0001.gmy
    │   ├── C0001_input.xml
    │   └── ...
    ├── outputs_csv/            HemeLB simulation output
    │   └── C0001/
    │       └── C0001_fluid_data.csv
    ├── wss/                    WSS computation output
    │   └── C0001_wss.csv
    ├── dataset_manifest.csv    Aggregated AneuRisk metadata
    └── excluded_cases.csv      Quality exclusion log
```

---

## Timing (MacBook Air M5, 16GB RAM, MPI_CORES=4)

| Stage | Time per case |
|-------|--------------|
| PR2 generation | < 1 s |
| Voxelization (0.2mm) | ~8 s |
| HemeLB simulation (20000 steps) | ~25-35 min (geometry-dependent) |
| CSV extraction | ~1-2 min |
| Quality checks | ~15 s |
| **Total per case** | **~27-37 min** |
| **Total for 95 cases** | **~45-60 hours** |

Split across two or more overnight runs. Pipeline is safe to interrupt and resume.

---

## Known Limitations

1. **Steady-state only.** Pulsatile flow requires Womersley BCs and ~240,000 steps at dt=1e-5s, which causes LBM instability on this hardware at 0.2mm resolution.

2. **Non-physiological inlet velocity.** v_mean = 0.025 m/s vs physiological 0.3-0.5 m/s for ICA. Required to stay below Mach limit (Ma < 0.1) at dx=0.2mm. GNN trains on normalised WSS, so this does not affect spatial pattern learning.

3. **Newtonian blood.** mu = 0.0035 Pa.s constant. Non-Newtonian effects (Carreau-Yasuda) are important at low shear rates near stagnation zones.

4. **Rigid walls.** Fluid-structure interaction not modelled.

5. **0.2mm resolution.** Borderline for WSS accuracy (10-20 voxels across vessel diameter). Published automated workflows use 0.2-0.27mm range.

6. **Retrospective geometry.** Ruptured aneurysms may have changed shape post-rupture. The geometry simulated is post-event, not pre-rupture.

7. **Equal outlet pressures.** Justified under Murray's law self-similar tree assumption but not patient-specific.

All limitations must be stated explicitly in the methods section.

---

## Next Steps (Remaining Pipeline)

### Immediate
- [ ] Complete 95-case batch with velocity BCs
- [ ] Rewrite `compute_wss.py` for STL-vertex-mapped output (Step 3 depends on this)
- [ ] Check class balance (ruptured vs unruptured) in `dataset_manifest.csv`

### Step 4: Graph Construction
- STL vertices as graph nodes
- Edge connectivity from STL faces
- Node features: position, local curvature, geodesic distance to inlet, vessel radius from centerlines
- Graph label: per-vertex normalised WSS from updated `compute_wss.py`

### Step 5: GNN Training
- Architecture: Message Passing GNN (MP-GNN) or Graph Transformer
- Input: patient geometry graph
- Output: per-node normalised WSS
- Loss: MSE on nWSS
- Validation: held-out cases

### Step 6: Clinical Output
- Rupture risk score from WSS distribution features
- Correlation with `ruptureStatus` from manifest
- Stratification by `aneurysmLocation`

### Step 7: Validation and Paper
- Resolution sensitivity: compare 0.2mm vs 0.1mm WSS on 3 cases
- Surrogate accuracy: compare GNN-predicted vs HemeLB-computed WSS
- Clinical correlation: WSS features vs rupture status

---

## Citation

If using this pipeline or dataset, cite:

**AneuRisk Database:**
```
Sangalli LM, et al. (2014). Functional data analysis for haemodynamics.
AneuRisk Web Repository: http://ecm2.mathcs.emory.edu/aneuriskweb
```

**HemeLB:**
```
Mazzeo MD, Coveney PV. (2008). HemeLB: A high performance parallel lattice-Boltzmann
code for large scale fluid flow in complex geometries.
Computer Physics Communications. https://doi.org/10.1016/j.cpc.2008.02.013
```

---

## Hardware

Developed on MacBook Air M5, 16GB RAM. All simulation parameters are tuned for this hardware constraint. Running on an HPC cluster would allow 0.1mm resolution, pulsatile BCs, and full convergence in a fraction of the time.

---

## Contact

Joyanta Mondal
