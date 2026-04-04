# Quick Start Guide

## Prerequisites

- **Chaste** compiled (or Docker container `chaste/release`)
- CMake ≥ 3.16, GCC/Clang with C++17
- Python 3.8+ (for analysis scripts)

## Build

All commands from the CMake build root (`/home/orlando/Thesis/`):

```bash
cmake .                                    # Detect new .cpp files
make -j$(nproc) CryptBuddingApp            # Build main executable
```

After adding new `.cpp` files to `src/`, re-run `cmake .` before `make`.
Header-only `.hpp` files do **not** require a cmake re-run.

## Run

```bash
# Default vertex2d simulation:
./projects/TissueMorphology/apps/CryptBuddingApp \
  -config Chaste/projects/TissueMorphology/apps/params-Vertex2d.json

# With shear modulus (Pa) — derives SLS parameters automatically:
./projects/TissueMorphology/apps/CryptBuddingApp \
  -config apps/params-Vertex2d.json \
  -model vertex2d \
  -shear-modulus 1300 \
  -run 0
```

### CLI flags

| Flag | Description | Default |
|------|-------------|---------|
| `-model <type>` | `node2d`, `vertex2d`, `node3d`, `vertex3d` | Required |
| `-config <path>` | JSON parameter file | Required |
| `-shear-modulus <Pa>` | ECM shear modulus; derives `ghostE0`, `ghostE1` | — |
| `-stiffness <val>` | Cell-ghost spring stiffness (sim units) | 5.0 |
| `-run <int>` | Replicate index (affects random seed) | 0 |
| `-lumen <0\|1>` | Enable lumen pressure | 1 |
| `-apical <0\|1>` | Enable apical constriction | 1 |

CLI flags override JSON values.

## Configuration (JSON)

Four model-specific configs are provided:
- `apps/params-Node2d.json`
- `apps/params-Vertex2d.json`
- `apps/params-Node3d.json`
- `apps/params-Vertex3d.json`

Structure:
```json
{
  "unitConversion": { "usePhysicalUnits": false },
  "simulation":     { "modelType": "vertex2d", "endTime": 168 },
  "features":       { "enableLumenPressure": true, "enableGhostNodeECM": true },
  "geometry":       { "numCells2dVertex": 12 },
  "forces":         { "SpringForce": {}, "GhostNodeECM": { "ViscoelasticECM": {} } },
  "cellTypes":      { "stemFraction": 0.42 },
  "cellCycle":      { "quiescentFraction": 0.7 }
}
```

Missing keys silently keep compiled-in defaults. Unknown keys are ignored.

### Physical units mode

Set `"usePhysicalUnits": true` in `unitConversion` and specify parameters in SI (N/m, Pa).
See [docs/units/](../units/) for conversion details.

### Shear modulus → SLS derivation

When `-shear-modulus G` is provided (or `ecmShearModulusPa` in JSON):
- $E_0 = G \times 2.9 \times \Delta x / r$, where $r = E_1/E_0 = 0.05$ (Fertala et al., 2025)
- $E_1 = r \times E_0$
- Conversion: $G_u \approx (E_0 + E_1) / (2.9a)$

## Output

Results go to `testoutput/CryptBudding/<git-hash>/<model>/<ecm-label>/run_<N>/`:
- `.vtu` / `.pvd` — VTK time series (ParaView)
- `summary.csv` — per-timestep metrics from `CryptBuddingSummaryModifier`
- `params_used.json` — snapshot of all parameters used

## Analysis

```bash
cd experiments/ECMStiffnessSweep/
python3 analyse_crypt_budding.py \
  --base-dir /path/to/testoutput/CryptBudding/<hash> \
  --model vertex2d
```

## Tests

```bash
make -j$(nproc) TestCryptBudding && ctest -R TestCryptBudding -V
```

Environment variables: `ECM_STIFFNESS`, `RUN_NUMBER`.

## Adding a new parameter

1. Add member + default in `CryptBuddingParams.hpp` (`SetDefaults()`)
2. Add JSON read in `LoadFromFile()`
3. Add output in `SaveToFile()`
4. Use in the relevant `Run*.hpp` file(s)
