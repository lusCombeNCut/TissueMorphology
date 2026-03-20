# TissueMorphology Project

A Chaste-based computational biology project for studying organoid formation, tissue morphology, and intestinal crypt budding under varying ECM stiffness.

## Project Structure

```
TissueMorphology/
├── src/                                    # Shared library classes (forces, writers, cell cycle, mutation states)
│   ├── StochasticFourTypeCellCycleModel.hpp    # Stochastic SC/TA/PC/EC cell cycle (Montes-Olivas 2023)
│   ├── UniformContactInhibitionCellCycleModel.hpp
│   ├── UniformContactInhibitionGenerationalCellCycleModel.hpp
│   ├── TACellMutationState.{hpp,cpp}           # MutationState: Transit-Amplifying
│   ├── PanethCellMutationState.{hpp,cpp}       # MutationState: Paneth Cell
│   ├── EnterocyteCellMutationState.{hpp,cpp}   # MutationState: Enterocyte
│   ├── ECMConfinementForce.hpp                 # Cell-ECM spring confinement (2D)
│   ├── DynamicECMField.hpp                     # ECM density field (2D)
│   ├── RingSmoothingForce.hpp                  # Discrete Laplacian bending (node2d)
│   ├── RingSpringForce.hpp                     # Topology-based springs (node2d)
│   ├── RingTopologyTracker.hpp                 # Ring neighbour tracking (node2d)
│   ├── FastNagaiHondaForce.hpp                 # Nagai-Honda with per-type adhesion (vertex2d)
│   ├── LumenPressureForce.hpp                  # Hydrostatic lumen pressure
│   ├── DifferentialAdhesionForce.hpp           # Cell-type-dependent adhesion
│   ├── ApicalConstrictionForce.hpp             # Apical surface contractility
│   ├── CellPolarityForce.hpp                   # Apicobasal polarity force
│   ├── BasementMembraneForce.hpp               # Radial BM constraint (3D)
│   ├── ...                                     # + 3D forces, writers, modifiers
│   └── _unused/                                # Archived unused classes
├── apps/
│   ├── src/
│   │   ├── CryptBuddingApp.cpp                 # Main entry point (CLI + INI dispatch)
│   │   ├── CryptBuddingParams.hpp              # INI parser + parameter struct + autotuning
│   │   ├── RunNode2d.hpp                       # 2D node-based runner
│   │   ├── RunVertex2d.hpp                     # 2D vertex-based runner
│   │   ├── RunNode3d.hpp                       # 3D node-based runner
│   │   ├── RunVertex3d.hpp                     # 3D vertex-based runner
│   │   ├── CryptBuddingSummaryModifier.hpp     # Per-timestep summary CSV writer
│   │   └── CryptBuddingUtils.hpp               # Shared utilities (cell assignment, killers)
│   ├── params-Node2d.ini                       # Tuned parameters for 2D node model
│   ├── params-Vertex2d.ini                     # Tuned parameters for 2D vertex model
│   ├── config/default_params.ini               # Template matching SetDefaults()
│   └── README.md                               # Configuration reference
├── test/
│   ├── CryptBudding/TestCryptBudding.hpp       # Integration test
│   ├── ECMForces/                              # ECM force unit tests
│   ├── Invasion/                               # Invasive front experiments
│   └── ContinuousTestPack.txt                  # CTest registration
├── docs/                                       # Constitutive equations + design docs
├── experiments/                                # Sweep scripts + analysis
├── hpc/                                        # HPC submission scripts
├── scripts/                                    # Utility scripts
└── visualization/                              # Python visualization tools
```

## HPC Pipeline (BluePebble)

### Prerequisites

| Item | Path on HPC |
|------|-------------|
| Container | `/user/work/$(whoami)/containers/tissuemorphology.sif` |
| Source (auto-cloned) | `/user/work/$(whoami)/TissueMorphology` |
| Build dir (auto-seeded) | `/user/work/$(whoami)/chaste_build` |
| Logs | `/user/work/$(whoami)/logs/crypt_budding/` |
| Results | `/user/work/$(whoami)/chaste_output/CryptBudding_*_results/` |

### How It Works

The sweep script uses a **two-phase Slurm pipeline** to avoid shared-library corruption:

1. **BUILD job** (single, ~1 hr) — `cmake` + `make` of the test binary
2. **RUN job** (array of 70, `--dependency=afterok`) — 7 stiffness levels × 10 replicates, simulations only

On any task failure, `scancel` cancels all remaining siblings. The last completing task creates a zip archive.

### Common Commands

```bash
# --- On login node ---

# Submit node-based stiffness sweep (build + 70 run tasks)
cd /user/work/$(whoami)/TissueMorphology/hpc
./submit_crypt_budding_sweep.sh node

# Submit vertex-based sweep
./submit_crypt_budding_sweep.sh vertex

# Submit both
./submit_crypt_budding_sweep.sh both

# Monitor jobs
squeue -u $(whoami)

# Check a specific run log
cat /user/work/$(whoami)/logs/crypt_budding/node_s0.5_r0_*.log | tail -30

# Find failed jobs
grep -l "JOB FAILED" /user/work/$(whoami)/logs/crypt_budding/*.log

# Cancel all your jobs
scancel -u $(whoami)

# Quick compile-only check (30 min, no simulation)
sbatch submit_compile_test.sh

# Nuke build dir to force clean rebuild
rm -rf /user/work/$(whoami)/chaste_build
```

### Updating Code on HPC

The sweep script auto-pulls `origin/main` at job start. Just push locally:

```bash
# On your machine
cd Chaste/projects/TissueMorphology
git add -A && git commit -m "description" && git push
# Then submit on HPC — it will git pull automatically
```

### Downloading Results

```bash
# From your local machine
scp -r sv22482@bp1-login.acrc.bris.ac.uk:/user/work/sv22482/chaste_output/CryptBudding_node_results ./
# Or download the zip archive (created when all 70 tasks finish)
scp sv22482@bp1-login.acrc.bris.ac.uk:/user/work/sv22482/chaste_output/CryptBudding_node_*.zip ./
```

### Post-Processing

```bash
python3 scripts/analyse_crypt_budding.py \
  --node-dir ./CryptBudding_node_results \
  --vertex-dir ./CryptBudding_vertex_results \
  --output-dir ./analysis_output
```

Produces: crypt count vs stiffness box plots, mean±SD plots, CSV summary, model comparison overlay.

## Crypt Budding Simulations

### Research Question

> Do changes in stiffness of in-silico hydrogel simulations have the same effect as changes in in-vivo ECM stiffness?

### Cell Cycle Model

The simulation uses the **Stochastic Four-Type Cell Cycle Model** (Montes-Olivas et al. 2023)
with contact inhibition, supporting four cell types: Stem (SC), Transit-Amplifying (TA),
Paneth (PC), and Enterocyte (EC). See `docs/ConstitutiveEquations_CellCycles.md` for full
equations and `apps/README.md` for parameter mapping.

### Key Parameters

| Parameter | Node-2D | Vertex-2D |
|-----------|---------|-----------|
| Initial cells | 18 (ring) | 12 (annulus) |
| Cell cycle | `StochasticFourTypeCellCycleModel` | same |
| Quiescent fraction | 0.7 | 0.5 |
| ECM confinement stiffness | 50 | 3.0 |
| Spring stiffness | 15 | 5 (Nagai-Honda) |
| Lumen mode | target-volume (K=30) | target-volume (K=500) |
| ECM grid | hex, spacing=1 | hex, spacing=0.5 |
| End time | 48 h (local) / sweep-dependent | 168 h |
| dt | 0.005 | 0.005 (overrides autotuning) |

See `apps/params-Node2d.ini` and `apps/params-Vertex2d.ini` for full parameter listings.

### Environment Variables

| Variable | Description | Set by |
|----------|-------------|--------|
| `ECM_STIFFNESS` | Basement membrane stiffness value | Sweep script |
| `RUN_NUMBER` | Replicate index (0–9) | Sweep script |
| `CHASTE_TEST_OUTPUT` | Output directory | Sweep script |

## Building Locally

```bash
cd /home/orlando/Thesis          # build directory
cmake Chaste -DChaste_ERROR_ON_WARNING=OFF
make -j$(nproc) Test2dCryptBuddingNodeBased
make -j$(nproc) Test2dCryptBuddingVertexBased
```

## Running Locally

```bash
# Single run with specific stiffness
ECM_STIFFNESS=5.0 RUN_NUMBER=0 ctest -R Test2dCryptBuddingNodeBased -V
```
