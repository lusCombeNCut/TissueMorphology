# CryptBudding Experiments

Stiffness sweep experiments for crypt budding / organoid morphogenesis using the
**CryptBuddingApp** standalone executable (no CTest).

## Directory Structure

```
experiments/CryptBudding/
├── submit_sweep.sh              # SLURM array job script
├── analyse_crypt_budding.py     # Post-processing & crypt detection
└── README.md                    # This file
```

## Quick Start (BluePebble HPC)

### 1. Pull the container (one-time)
```bash
module load apptainer
apptainer pull /user/work/$(whoami)/containers/tissuemorphology.sif \
    docker://ghcr.io/luscombencut/tissuemorphology:latest
```

### 2. Submit a sweep
```bash
# Single model type (70 jobs: 7 stiffness × 10 replicates)
sbatch submit_sweep.sh node2d
sbatch submit_sweep.sh vertex2d
sbatch submit_sweep.sh node3d
sbatch submit_sweep.sh vertex3d

# All four models at once (4 × 70 = 280 jobs)
./submit_sweep.sh all

# Just 2D or 3D models
./submit_sweep.sh all2d
./submit_sweep.sh all3d
```

### 3. Monitor jobs
```bash
squeue -u $(whoami) -o "%.10i %.8j %.2t %.10M %.6D %R"
```

### 4. Analyse results
```bash
# All models found automatically
python analyse_crypt_budding.py \
    --base-dir /user/work/$(whoami)/sim_output \
    -o analysis_results/

# Specific model
python analyse_crypt_budding.py \
    --base-dir /user/work/$(whoami)/sim_output \
    --model node2d vertex2d \
    -o analysis_results/
```

## Sweep Parameters

| Parameter     | Values                              |
|---------------|-------------------------------------|
| **Stiffness** | 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0 |
| **Replicates** | 10 per stiffness level             |
| **Models**    | node2d, vertex2d, node3d, vertex3d  |
| **Array size** | 70 jobs per model (7 × 10)         |

## CryptBuddingApp CLI

The app replaces the old env-var-driven CTest approach. Full usage:
```
./CryptBuddingApp -model <type> [options]

Required:
  -model <node2d|vertex2d|node3d|vertex3d>

Options:
  -stiffness <double>   ECM stiffness (default: 5.0)
  -run <int>            run/replicate number (default: 0)
  -lumen <0|1>          lumen pressure (default: 1)
  -apical <0|1>         apical constriction (default: 1)
  -ecm <0|1>            ECM guidance, 3D only (default: 0)
  -relax <0|1>          relaxation phase (default: 1)
  -slough <0|1>         sloughing, 2D only (default: 1)
  -diffadh <0|1>        differential adhesion (default: 1)
  -endtime <double>     override end time
  -dt <double>          override timestep
  -help                 print this message
```

## Output Structure

The app writes output to `$CHASTE_TEST_OUTPUT/CryptBudding/<model>/stiffness_<X>/run_<N>/`:
```
sim_output/
└── CryptBudding/
    ├── node2d/
    │   ├── stiffness_0.5/
    │   │   ├── run_0/
    │   │   │   ├── results.viznodes
    │   │   │   ├── results_*.vtu
    │   │   │   └── crypt_summary.csv
    │   │   ├── run_1/
    │   │   └── ...
    │   └── stiffness_50.0/
    └── vertex3d/
        └── ...
```

## Changes from Old Approach

| Feature | Old (CTest) | New (App) |
|---|---|---|
| **Executable** | `ctest -R TestCryptBudding` | `./CryptBuddingApp -model node2d` |
| **Parameters** | Environment variables | CLI flags `-flag value` |
| **Output** | stdout buffered by CTest | Real-time stdout |
| **Framework** | CxxTest suite | Standalone `main()` |
| **Build target** | `make TestCryptBudding` | `make CryptBuddingApp` |
