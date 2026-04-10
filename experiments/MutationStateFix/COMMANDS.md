# MutationStateFix — Experiment Commands

Re-run vertex2d + node2d at **300 Pa** and **1300 Pa** (5 replicates each)
after fixing the `CellMutationStatesWriter` / `LegacyCellProliferativeTypesWriter`
CellLabel colour-override bug that caused quiescent cells to be mislabelled as
Paneth cells (CellLabel default colour = 5 = PanethCellMutationState colour).

## Bug Summary

| File | Bug | Fix |
|---|---|---|
| `cell_based/src/writers/cell_writers/CellMutationStatesWriter.cpp` | `GetCellDataForVtkOutput()` and `VisitCell()` replaced mutation state with `CellLabel::GetColour()` (=5) when contact-inhibited | Removed CellLabel override — always output `GetMutationState()->GetColour()` |
| `cell_based/src/writers/cell_writers/LegacyCellProliferativeTypesWriter.cpp` | Same CellLabel override in both methods | Same fix |

## 1. Commit and push the fix

```bash
cd /home/orlando/Thesis/Chaste/projects/TissueMorphology

# Stage the Chaste core writer fixes (these live in the parent Chaste repo)
cd /home/orlando/Thesis/Chaste
git add cell_based/src/writers/cell_writers/CellMutationStatesWriter.cpp
git add cell_based/src/writers/cell_writers/LegacyCellProliferativeTypesWriter.cpp
git commit -m "fix: remove CellLabel colour override from VTU mutation state writers

ContactInhibitionCellCycleModel adds a CellLabel (default colour=5) to
quiescent cells. CellMutationStatesWriter and LegacyCellProliferativeTypesWriter
both checked for CellLabel and substituted its GetColour() value into the
output, overriding the actual mutation state. Since CellLabel colour=5
coincides with PanethCellMutationState, every contact-inhibited cell appeared
as Paneth in VTU files, causing oscillating/inflated Paneth counts."

# Now commit the experiment files to TissueMorphology
cd /home/orlando/Thesis/Chaste/projects/TissueMorphology
git add experiments/MutationStateFix/
git commit -m "exp: add MutationStateFix re-run experiment (300Pa, 1300Pa)"
git push
```

## 2. SSH to BluePebble and update the source

```bash
ssh sv22482@bp1-login.acrc.bris.ac.uk

cd /user/work/sv22482/TissueMorphology
git fetch origin && git reset --hard origin/main

# IMPORTANT: Also update the Chaste core writer fix.
# The Chaste build dir has the old source. Copy the fixed files:
cp /user/work/sv22482/TissueMorphology/../chaste_src/cell_based/src/writers/cell_writers/CellMutationStatesWriter.cpp \
   /user/work/sv22482/chaste_build/cell_based/src/writers/cell_writers/ 2>/dev/null || true

# Or, if Chaste source is inside the container, bind-mount the fix.
# The build phase will rebuild from the updated TissueMorphology source,
# but the Chaste core files need to be updated in the build tree manually
# (or in the container source). See step 2b below.
```

### 2b. Patch Chaste core writers in the build tree

The Chaste core `.cpp` files live inside the container image at
`/home/chaste/src/cell_based/src/writers/cell_writers/`. Since the build
job bind-mounts `/user/work/sv22482/chaste_build` as `/home/chaste/build`,
you need to ensure the **source** used by cmake also has the fix.

**Option A** — Patch the build tree source directly:

```bash
# Copy fixed files into the persistent build tree source cache
BUILD_SRC="/user/work/sv22482/chaste_build"

# The container cmake copies source into the build tree on first run.
# Find where the writers are compiled from:
find "${BUILD_SRC}" -name "CellMutationStatesWriter.cpp" -type f 2>/dev/null

# Then overwrite those files with the fixed versions from TissueMorphology:
# (adjust path based on find output)
```

**Option B** — Bind-mount fixed source over container source:

Add to the apptainer exec calls in submit.sh:
```bash
--bind "/user/work/sv22482/TissueMorphology/chaste_fixes/cell_based:/home/chaste/src/cell_based"
```

**Option C (simplest)** — Delete the build tree and let the build phase
re-seed from container + rebuild with the updated writers:

```bash
rm -rf /user/work/sv22482/chaste_build
```

The build phase will re-seed from the container and rebuild everything.
This is slowest (~20 min) but guarantees a clean build.

## 3. Submit the experiment

```bash
cd /user/work/sv22482/TissueMorphology/experiments/MutationStateFix

# Run both vertex2d and node2d
./submit.sh all2d

# Or run individually:
# ./submit.sh vertex2d
# ./submit.sh node2d
```

This submits:
- 1 build job (shared)
- 10 simulation tasks per model (2 stiffness × 5 replicates)
- 1 archive job (creates `.tar.gz` after all sims finish)

## 4. Monitor jobs

```bash
squeue -u sv22482
# or watch specific jobs:
watch -n 30 'squeue -u sv22482 --format="%.10i %.15j %.8T %.10M %.6D %R"'
```

## 5. Download results

The archive job prints the exact scp command, but it will be:

```bash
# From local machine:
scp sv22482@bp1-login.acrc.bris.ac.uk:/user/work/sv22482/sim_output/mutation_state_fix_<TIMESTAMP>.tar.gz \
    /home/orlando/Thesis/sim_results/

# Extract:
cd /home/orlando/Thesis/sim_results
tar -xzf mutation_state_fix_<TIMESTAMP>.tar.gz
```

## 6. Run analysis locally

```bash
cd /home/orlando/Thesis/Chaste/projects/TissueMorphology/experiments
python3 plot_cell_type_counts.py
```

(Update data paths in the script to point at the new extracted results.)
