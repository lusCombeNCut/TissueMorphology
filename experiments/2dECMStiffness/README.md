# Experiment: 2D Crypt Budding — ECM Stiffness Sweep

## Aim

Determine how **ECM (basement membrane) stiffness** controls the number of
**crypt buds** formed during intestinal organoid development. The hypothesis is
that softer ECM permits more outward budding events while stiffer ECM suppresses
them.

## Models

Two complementary 2D cell-based models simulate a cross-section of a circular
organoid:

| Model | Test file | Framework |
|-------|-----------|-----------|
| **Node-based** | `test/CryptBudding2d/Test2dCryptBuddingNodeBased.hpp` | Overlapping spheres, linear springs, differential adhesion |
| **Vertex-based** | `test/CryptBudding2d/Test2dCryptBuddingVertexBased.hpp` | Nagai–Honda vertex model (area + perimeter energy) |

Both share basement membrane confinement, lumen pressure, and apical
constriction forces. Full constitutive equations and modelling choices are
documented in `test/CryptBudding2d/experiment_design.md`.

## Sweep Design

| Property | Value |
|----------|-------|
| Independent variable | `ECM_STIFFNESS` |
| Levels | 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0 |
| Replicates per level | 10 |
| Total jobs per model | 70 (7 × 10) |
| Total jobs (both models) | 140 |
| Simulation duration | 200 h (node), 200 h (vertex) |
| Primary output | Number of crypt buds detected |

## Files in This Directory

| File | Purpose |
|------|---------|
| `submit_crypt_budding_sweep.sh` | Slurm array job script for BluePebble HPC |
| `analyse_crypt_budding.py` | Post-processing: crypt detection + plotting |
| `README.md` | This file |

## Usage

### On HPC (BluePebble)

```bash
# Node-based sweep (70 jobs)
sbatch experiments/2dECMStiffness/submit_crypt_budding_sweep.sh node

# Vertex-based sweep (70 jobs)
sbatch experiments/2dECMStiffness/submit_crypt_budding_sweep.sh vertex

# Both models (submits two array jobs)
./experiments/2dECMStiffness/submit_crypt_budding_sweep.sh both
```

### Local (single stiffness)

```bash
ECM_STIFFNESS=5.0 RUN_NUMBER=0 ctest -R Test2dCryptBuddingNodeBased --output-on-failure
ECM_STIFFNESS=5.0 RUN_NUMBER=0 ctest -R Test2dCryptBuddingVertexBased --output-on-failure
```

### Analysis

```bash
# Single model
python experiments/2dECMStiffness/analyse_crypt_budding.py \
    --data-dir /path/to/testoutput/CryptBudding2d_NodeBased

# Both models (overlay comparison)
python experiments/2dECMStiffness/analyse_crypt_budding.py \
    --node-dir /path/to/CryptBudding2d_NodeBased \
    --vertex-dir /path/to/CryptBudding2d_VertexBased \
    --output-dir analysis_output
```

### Output

- `crypts_vs_stiffness.png` — box plot + mean ± SD
- `comparison_crypts_vs_stiffness.png` — node vs vertex overlay
- `crypt_counts.csv` — raw counts per (stiffness, replicate)

## Crypt Detection Method

1. Compute population centroid
2. Convert cell positions to polar coordinates $(r, \theta)$
3. Bin by $\theta$, compute mean $\bar{r}(\theta)$
4. Smooth with Savitzky–Golay filter (circular padding)
5. Detect peaks in $\bar{r}(\theta)$ via `scipy.signal.find_peaks`
   (prominence ≥ 0.5, angular separation ≥ 30°)
6. Number of peaks = number of crypts

## Related Documentation

- `test/CryptBudding2d/experiment_design.md` — full experiment design
- `test/CryptBudding2d/parameters.md` — all parameter values and justification status
