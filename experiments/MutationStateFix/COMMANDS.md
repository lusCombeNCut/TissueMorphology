# MutationStateFix — Experiment Commands

Re-run vertex2d + node2d at **300 Pa** and **1300 Pa** (5 replicates each)
with two bug fixes:

1. **Writer bug**: `CellTrueMutationStateWriter` replaces the default
   `CellMutationStatesWriter` (which overwrote mutation state with
   CellLabel colour=5 for quiescent cells, making them appear as Paneth).
2. **Enterocyte assignment bug**: `AssignStochasticFourType` now sets
   `EnterocyteCellMutationState` (ID=6) on the enterocyte branch —
   previously enterocytes kept `WildTypeCellMutationState` (ID=0),
   making them indistinguishable from stem cells. Also fixed
   `cell_type_id` for enterocytes (was 2.0, now 3.0).

## 1. Push

```bash
cd /home/orlando/Thesis/Chaste/projects/TissueMorphology
git add -A && git commit -m "fix: assign EnterocyteCellMutationState in AssignStochasticFourType + CellTrueMutationStateWriter"
git push
```

## 2. Submit on BluePebble

```bash
ssh sv22482@bp1-login.acrc.bris.ac.uk
cd /user/work/sv22482/TissueMorphology
git fetch origin && git reset --hard origin/main
cd experiments/MutationStateFix
./submit.sh all2d
```

## 3. Monitor

```bash
squeue -u sv22482
```

## 4. Download and extract

```bash
scp sv22482@bp1-login.acrc.bris.ac.uk:/user/work/sv22482/sim_output/mutation_state_fix_<TIMESTAMP>.tar.gz \
    /home/orlando/Thesis/sim_results/
cd /home/orlando/Thesis/sim_results && tar -xzf mutation_state_fix_<TIMESTAMP>.tar.gz
```

## 5. Analyse locally

```bash
cd /home/orlando/Thesis/Chaste/projects/TissueMorphology/experiments/MutationStateFix

# Update VERTEX_ROOT / NODE_ROOT paths in each script first, then:
python3 plot_cell_type_counts.py
python3 plot_quiescent_fraction.py
python3 plot_combined_grid.py
```

Outputs saved to `experiments/MutationStateFix/plots/`.
