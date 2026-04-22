# TissueMorphology

Chaste-based simulation of intestinal crypt budding under varying ECM stiffness.
Masters thesis — Engineering Mathematics, University of Bristol.

## Structure

```
apps/src/          — CryptBuddingApp.cpp, CryptBuddingParams.hpp, Run{Node,Vertex}{2d,3d}.hpp
src/               — Forces, cell cycle models, writers, mutation states
test/              — Integration + unit tests
docs/              — Equations, architecture, unit conversion, quickstart
experiments/       — Sweep scripts and analysis
hpc/               — BluePebble Slurm submission scripts
scripts/           — Utility scripts
latex_snippets/    — Thesis table/figure snippets
visualization/     — Python visualisation tools
```

## Build & Run (local)

```bash
cd /home/orlando/Thesis
cmake Chaste -DChaste_ERROR_ON_WARNING=OFF   # re-run only after adding new .cpp files
make -j$(nproc) CryptBuddingApp

./projects/TissueMorphology/apps/CryptBuddingApp \
  -config Chaste/projects/TissueMorphology/apps/params-Vertex2d.json \
  -model vertex2d -shear-modulus 1300 -run 0
```

See `docs/quickstart/` for more detail.

## HPC (BluePebble)

```bash
# On login node — submits build job + 70-task array (7 stiffness × 10 replicates)
cd /user/work/$(whoami)/TissueMorphology/hpc
./submit_crypt_budding_sweep.sh node     # node-based
./submit_crypt_budding_sweep.sh vertex   # vertex-based
./submit_crypt_budding_sweep.sh both

squeue -u $(whoami)                      # monitor
scancel -u $(whoami)                     # cancel all
```

The sweep script auto-pulls `origin/main` at job start — just push locally and submit.

**Results paths on HPC:**
- Logs: `/user/work/$(whoami)/logs/crypt_budding/`
- Output: `/user/work/$(whoami)/chaste_output/CryptBudding_*_results/`

```bash
# Download results
scp sv22482@bp1-login.acrc.bris.ac.uk:/user/work/sv22482/chaste_output/CryptBudding_node_*.zip ./
```

## API Docs

Generated with Doxygen from `docs/Doxyfile`:

```bash
cd docs && doxygen Doxyfile
# Opens docs/api/html/index.html
```

`docs/api/` is gitignored — regenerate locally as needed.

## Key docs

| Path | Contents |
|------|----------|
| `docs/ARCHITECTURE.md` | Dependency matrix, naming conventions |
| `docs/equations/` | Constitutive equations for all forces and cell cycles |
| `docs/units/` | CHASTE→SI unit conversion |
| `docs/PARAMETER_JUSTIFICATION.md` | Literature-sourced parameter validation |
| `apps/README.md` | JSON parameter → C++ class mapping |
