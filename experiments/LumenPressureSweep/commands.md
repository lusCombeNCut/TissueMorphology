# LumenPressureSweep — Run Commands

Re-run to collect both `node2d` and `vertex2d` output data.

---

## 1. Submit the sweep (BluePebble login node)

```bash
ssh sv22482@bp1-login.acrc.bris.ac.uk
cd /user/work/sv22482/TissueMorphology/experiments/LumenPressureSweep

./submit_sweep.sh all2d
```

This submits three chained Slurm jobs:
- **Build** — git pull + cmake + make CryptBuddingApp
- **node2d array** — 8 pressures × 10 replicates = 80 tasks (depends on build)
- **vertex2d array** — 8 pressures × 10 replicates = 80 tasks (depends on build)
- **Analysis** — runs after both arrays finish

The script prints the exact `scp` command with the timestamped output path when submitted. Note it down.

---

## 2. Monitor jobs

```bash
squeue -u sv22482
```

---

## 3. Download results locally

Replace `<TIMESTAMP>` with the value printed during submission (format: `YYYY-MM-DD_HH-MM-SS`).

### node2d

```bash
mkdir -p /home/orlando/Thesis/sim_results/lumen_pressure_sweep

scp -r sv22482@bp1-login.acrc.bris.ac.uk:/user/work/sv22482/sim_output/<JOBID>_node2d_<TIMESTAMP>/archives/ \
    /home/orlando/Thesis/sim_results/lumen_pressure_sweep/node2d_<TIMESTAMP>/
```

### vertex2d

```bash
scp -r sv22482@bp1-login.acrc.bris.ac.uk:/user/work/sv22482/sim_output/<JOBID>_vertex2d_<TIMESTAMP>/archives/ \
    /home/orlando/Thesis/sim_results/lumen_pressure_sweep/vertex2d_<TIMESTAMP>/
```

### Analysis plots

```bash
scp -r sv22482@bp1-login.acrc.bris.ac.uk:/user/work/sv22482/sim_output/analysis_lumen_pressure_<TIMESTAMP>/plots/ \
    /home/orlando/Thesis/sim_results/lumen_pressure_sweep/plots_<TIMESTAMP>/
```

---

## 4. Extract archives and run local analysis

```bash
cd /home/orlando/Thesis/Chaste/projects/TissueMorphology/experiments/LumenPressureSweep

# node2d
./fetch_and_analyse.sh \
    --skip-scp \
    --output-base /home/orlando/Thesis/sim_results/lumen_pressure_sweep \
    --model node2d

# vertex2d
./fetch_and_analyse.sh \
    --skip-scp \
    --output-base /home/orlando/Thesis/sim_results/lumen_pressure_sweep \
    --model vertex2d
```

Or use the combined fetch+extract approach (skipping the separate scp step above):

```bash
./fetch_and_analyse.sh \
    --job-dir <JOBID>_node2d_<TIMESTAMP> \
    --output-base /home/orlando/Thesis/sim_results/lumen_pressure_sweep \
    --model node2d

./fetch_and_analyse.sh \
    --job-dir <JOBID>_vertex2d_<TIMESTAMP> \
    --output-base /home/orlando/Thesis/sim_results/lumen_pressure_sweep \
    --model vertex2d
```
