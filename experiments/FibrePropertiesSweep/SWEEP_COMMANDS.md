# Fibre Anisotropy Strength Sweep — Commands

Sweeps `ghostAnisotropyStrength` (0.0, 0.1, 0.25, 0.5, 0.75, 1.0, 2.0) × 10 replicates = 70 tasks.

---

## 1. Submit on BluePebble

```bash
ssh sv22482@bp1-login.acrc.bris.ac.uk
cd /user/work/sv22482/TissueMorphology/experiments/FibrePropertiesSweep
./submit_sweep.sh vertex2d
```

Note the job IDs printed. The sweep runs in three automatic phases (build → sim → analyse).

---

## 2. Check job status

```bash
squeue -u sv22482
```

---

## 3. Merge scattered output into one folder (on BluePebble)

Replace `<JOBID>` with the sim array job ID printed at submission.

```bash
MERGED="/user/work/sv22482/sim_output/anisotropy_sweep_vertex2d_merged"
mkdir -p "$MERGED"

for dir in /user/work/sv22482/sim_output/<JOBID>_vertex2d_*/; do
    for run in "$dir"a*_r*/; do
        [ -d "$run" ] || continue
        cp -r "$run" "$MERGED/$(basename "$run")"
    done
done

cd /user/work/sv22482/sim_output
tar czf anisotropy_sweep_vertex2d_merged.tar.gz anisotropy_sweep_vertex2d_merged/
echo "Done — $(ls $MERGED | wc -l) runs merged"
```

Expected: 70 directories (`a0.0_r0` … `a2.0_r9`).

---

## 4. Copy to local machine

```bash
LOCAL_OUT="/home/orlando/Thesis/sim_results/FibrePropertiesSweep"
mkdir -p "$LOCAL_OUT"
scp sv22482@bp1-login.acrc.bris.ac.uk:/user/work/sv22482/sim_output/anisotropy_sweep_vertex2d_merged.tar.gz "$LOCAL_OUT/"
```

---

## 5. Extract and analyse locally

```bash
cd "$LOCAL_OUT"
tar xzf anisotropy_sweep_vertex2d_merged.tar.gz

python3 /home/orlando/Thesis/Chaste/projects/TissueMorphology/experiments/FibrePropertiesSweep/analyse_and_plot.py \
    --data-dir "$LOCAL_OUT/anisotropy_sweep_vertex2d_merged" \
    --model vertex2d \
    -o "$LOCAL_OUT/plots" \
    2>&1 | tee "$LOCAL_OUT/analysis.log"
```

Plots saved to `$LOCAL_OUT/plots/`.
