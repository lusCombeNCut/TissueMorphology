# TissueMorphology Project

A Chaste-based computational biology project for studying organoid formation, tissue morphology, and intestinal crypt budding under varying ECM stiffness.

## Project Structure

```
TissueMorphology/
├── src/
│   ├── BasementMembraneForce.hpp       # Radial BM constraint force
│   ├── DifferentialAdhesionForce.hpp    # Cell-type-dependent adhesion
│   └── OrganoidCellFactory.hpp         # Cell factory with BM properties
├── test/
│   ├── Test2dCryptBuddingNodeBased.hpp # 2D node-based crypt budding sweep
│   ├── Test2dCryptBuddingVertexBased.hpp # 2D vertex-based crypt budding sweep
│   ├── Test2dPainterReplication.hpp    # Painter cell migration replication
│   ├── Test3dCryptOrganoid.hpp         # 3D node-based crypt organoid
│   ├── Test3dVertexCryptOrganoid.hpp   # 3D vertex-based crypt organoid
│   ├── TestOrganoidFormation.hpp       # Basic organoid formation tests
│   └── ContinuousTestPack.txt          # Registered tests for CTest
├── hpc/
│   ├── submit_crypt_budding_sweep.sh   # Two-phase build→run sweep (main)
│   ├── submit_compile_test.sh          # Quick compile-only verification
│   ├── submit_test.sh                  # Single-test submission
│   └── submit_vertex_organoid.sh       # 3D vertex organoid submission
├── scripts/
│   ├── analyse_crypt_budding.py        # Post-processing & crypt counting
│   └── ...                             # Other analysis scripts
└── README.md
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

### Parameters

| Parameter | Node-based | Vertex-based |
|-----------|-----------|-------------|
| Mesh | 20×6 honeycomb → NodesOnlyMesh | 14×8 honeycomb vertex |
| Cell cycle | ContactInhibitionCellCycleModel | ContactInhibitionCellCycleModel |
| Quiescent fraction | 0.8 | 0.8 |
| BM stiffness | env `ECM_STIFFNESS` | 0.5 × env `ECM_STIFFNESS` |
| Spring stiffness | 30.0 | NagaiHonda (deformation = `ECM_STIFFNESS`) |
| Sloughing | PlaneBasedCellKiller at y=20 | PlaneBasedCellKiller at y=15 |
| End time | 200 hrs | 200 hrs |
| dt | 0.005 | 0.002 |
| Stiffness sweep | 0.5, 1, 2, 5, 10, 20, 50 | same |
| Replicates | 10 per stiffness | 10 per stiffness |

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
