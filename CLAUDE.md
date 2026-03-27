# TissueMorphology — Claude Code Context

## Project Overview

Chaste-based computational biology project simulating intestinal crypt budding and organoid formation under varying ECM (extracellular matrix) properties (stiffness, lumen pressure etc.). Part of a Masters thesis investigating whether in-silico hydrogel stiffness changes reproduce in-vivo ECM stiffness effects on morphogenesis.

**Project root:** `/home/orlando/Thesis/Chaste/projects/TissueMorphology`
**Build root:** `/home/orlando/Thesis` (Chaste CMake out-of-source build)

---

## Tech Stack

- **C++** — all simulation code (src/, apps/src/, test/)
- **CMake** — build system, integrates with Chaste core
- **Python** — post-processing and visualisation scripts
- **Bash** — HPC and utility scripts
- **Dependencies:** Chaste, OrganoidChaste (3D vertex), CGAL, VTK, PETSc, Boost

---

## Build System

```
/home/orlando/Thesis/                       ← CMake build root
├── Chaste/                                 ← Chaste source (read-only upstream)
│   └── projects/TissueMorphology/          ← THIS PROJECT (git-tracked)
│       ├── CMakeLists.txt                  ← chaste_do_project() — auto-discovers src/*.cpp
│       ├── src/                            ← Shared library: chaste_project_TissueMorphology
│       ├── apps/
│       │   ├── CMakeLists.txt              ← chaste_do_apps_project() — builds CryptBuddingApp
│       │   └── src/CryptBuddingApp.cpp     ← Main executable entry point
│       └── test/                           ← CxxTest-based tests
```

**Key build commands** (run from `/home/orlando/Thesis/`):
```bash
cmake .                                     # Re-detect new .cpp files
make -j$(nproc) CryptBuddingApp             # Build the main app
make -j$(nproc) chaste_project_TissueMorphology  # Build shared library only
make -j$(nproc) TestCryptBudding && ctest -R TestCryptBudding -V  # Build & run tests
```

**CRITICAL:** After adding new `.cpp` files to `src/`, you MUST re-run `cmake .` before `make` — the build system uses file globbing at configure time. Header-only `.hpp` files do NOT require a cmake re-run.

---

## Directory Structure

```
TissueMorphology/
├── src/                    # Shared library (48 headers)
│   ├── forces/             # ECM & tissue force classes
│   ├── cell_types/         # Mutation states (Paneth, TA, Enterocyte)
│   ├── cell_cycle/         # Cell cycle models
│   ├── division_rules/     # Division rules (node/vertex)
│   ├── topology/           # Ring/surface topology tracking
│   ├── modifiers/          # Simulation modifiers
│   ├── writers/            # VTK/CSV output writers
│   └── utilities/          # Helpers
├── apps/
│   ├── src/
│   │   ├── CryptBuddingApp.cpp        # Main entry point
│   │   ├── CryptBuddingParams.hpp     # JSON parser & parameter struct
│   │   ├── RunNode2d/Vertex2d/Node3d/Vertex3d.hpp
│   │   └── CryptBuddingSummaryModifier.*
│   ├── params-Node2d.json / params-Vertex2d.json / params-Node3d.json / params-Vertex3d.json
│   ├── params-Node2d-physical.json    # Physical-unit variants (SI)
│   └── experiments/        # JSON configs for parameter sweeps
├── test/
│   ├── CryptBudding/       # Unified integration tests
│   ├── ECMForces/          # Force unit tests
│   └── Invasion/           # Invasive front tests
├── experiments/            # Python analysis & sweep scripts
├── hpc/                    # BluePebble Slurm submission scripts
├── scripts/                # Utility Python scripts
└── visualization/          # Python visualisation tools
```

---

## Execution Flow

```
CryptBuddingApp.cpp
  ├── Parses CLI args: -model {node2d|vertex2d|node3d|vertex3d} -config <file.json>
  ├── CryptBuddingParams.hpp
  │     ├── SetDefaults()        ← compiled-in defaults
  │     ├── LoadFromFile()       ← JSON parser (nested objects, model type inside "simulation")
  │     └── Finalise()           ← auto-derives dt, endTime, t1/t2 thresholds, radii
  ├── Dispatches to one of:
  │     ├── RunNode2d.hpp        ← 2D node-based (ring of cells + ECM field)
  │     ├── RunVertex2d.hpp      ← 2D vertex-based (annulus + Nagai-Honda)
  │     ├── RunNode3d.hpp        ← 3D node-based (sphere surface)
  │     └── RunVertex3d.hpp      ← 3D vertex-based (OrganoidChaste monolayer)
  └── PvdFixUtils.hpp            ← Post-simulation PVD file repair
```

**Run command:**
```bash
# From the build root (/home/orlando/Thesis):
./projects/TissueMorphology/apps/CryptBuddingApp -config Chaste/projects/TissueMorphology/apps/params-Node2d.json

# CLI flags override JSON values, e.g.:
./projects/TissueMorphology/apps/CryptBuddingApp -config params-Node2d.json -model vertex2d
```

---

## Parameter System

Configuration is **JSON** (migrated from INI). The `modelType` field inside `"simulation"` selects the runner.

**JSON structure** (from `params-Node2d.json`):
```json
{
  "unitConversion": { "usePhysicalUnits": false, ... },
  "simulation":    { "modelType": "node2d", "endTime": 5.0, "dt": 0.005, ... },
  "features":      { "enableLumenPressure": true, "enableGhostNodeECM": true, ... },
  "geometry":      { "numCells2dNode": 50, "interactionCutoff2d": 1.5 },
  "forces": {
    "SpringForce":        { "springStiffness": 15.0, ... },
    "LumenPressureForce": { "lumenPressure": 5.0 },
    "RingSmoothingForce": { "bendingStiffness": 20.0, ... },
    "ECMConfinementForce":{ "ecmStiffness": 50.0, "ecmGridType": "hex", ... },
    "GhostNodeECM":       { "ghostDamping": 5.0, ..., "ViscoelasticECM": { ... } }
  },
  "cellTypes":  { "stemFraction": 0.2, "transitFraction": 0.5 },
  "cellCycle":  { "quiescentFraction": 0.6, ..., "StochasticFourTypeModel": { ... } }
}
```

**Physical units variant:** `params-*-physical.json` files set `"usePhysicalUnits": true` and define parameters in SI units (N/m, Pa, N·m). The effective damping η = referenceForcenN × 10⁻⁹ / (referenceVelocityUmPerH × 10⁻⁶/3600).

**Key parameter notes:**
- **Master stiffness knob:** `ecmStiffness` (under `ECMConfinementForce`) propagates to multiple forces
- **Auto-tuning in `Finalise()`:** `dt`, `endTime`, `t1Threshold`, `t2Threshold`, radii — setting these in JSON prevents the override
- Missing keys silently keep their `SetDefaults()` value; unknown keys are ignored
- **Units (non-physical mode):** time = hours, length = cell diameters (1 CD ≈ 10 µm), force = arbitrary (η=1)

**How to add a new parameter:**
1. Add member variable + comment in `apps/src/CryptBuddingParams.hpp`
2. Set default in `SetDefaults()`
3. Add JSON read call in `LoadFromFile()`
4. Add output line in `SaveToFile()`
5. Use in the relevant `Run*.hpp` file(s)
6. Document in `apps/params-Node2d.json` and/or other config files

**Parameter flow:**
```
JSON file ──► LoadFromFile() ──► struct members ──► Finalise() ──► Run*.hpp
```

---

## Key Source Files

| File | Purpose |
|------|---------|
| `apps/src/CryptBuddingApp.cpp` | CLI dispatcher — start here |
| `apps/src/CryptBuddingParams.hpp` | All simulation parameters, JSON parsing, auto-tuning |
| `apps/src/RunNode2d.hpp` | 2D node-based simulation setup |
| `apps/src/RunVertex2d.hpp` | 2D vertex-based simulation setup |
| `apps/src/RunNode3d.hpp` | 3D node-based simulation setup |
| `apps/src/RunVertex3d.hpp` | 3D vertex simulation setup (OrganoidChaste) |
| `src/forces/ECMConfinementForce.hpp` | Grid-based ECM confinement (primary ECM force) |
| `src/forces/GhostNodeEcmForce.hpp` | Particle-based viscoelastic ECM |
| `src/forces/ViscoelasticGhostNodeEcmForce.hpp` | Maxwell viscoelastic ECM variant |
| `src/forces/DynamicECMContactGuidanceForce.hpp` | Active migration along ECM fibres |
| `src/cell_cycle/StochasticFourTypeCellCycleModel.hpp` | Stem/TA/Paneth/Enterocyte stochastic transitions |

---

## Force Models

1. **ECMConfinementForce** (2D/3D) — Grid-based ECM with density-weighted springs; primary stiffness model
2. **DynamicECMContactGuidanceForce** (2D/3D) — Cell migration guidance along ECM fibres
3. **GhostNodeEcmForce** (2D/3D) — Particle-based viscoelastic ECM
4. **ViscoelasticGhostNodeEcmForce** — Maxwell viscoelastic variant
5. **LumenPressureForce** — Hydrostatic lumen expansion
6. **ApicalConstrictionForce** — Apical surface contractility
7. **CellPolarityForce** — Apicobasal polarity
8. **BasementMembraneForce** — Radial BM constraint
9. **DifferentialAdhesionForce** — Cell-type-dependent adhesion
10. **FastNagaiHondaForce** — Vertex model with per-type adhesion

---

## Dependency Graph (src/ → Runner)

| Class | node2d | vertex2d | node3d | vertex3d | Tests |
|-------|:------:|:--------:|:------:|:--------:|:-----:|
| **Cell Cycle** |
| StochasticFourTypeCellCycleModel | ✓ | ✓ | — | — | — |
| UniformContactInhibitionGenerationalCCM | ✓ | ✓ | ✓ | ✓ | — |
| UniformContactInhibitionCellCycleModel | ✓ | ✓ | ✓ | ✓ | — |
| **Mutation States** |
| TACellMutationState | ✓ | ✓ | — | — | — |
| PanethCellMutationState | ✓ | ✓ | — | — | — |
| EnterocyteCellMutationState | ✓ | ✓ | — | — | — |
| **Forces (2D)** |
| RingSpringForce | ✓ | — | — | — | — |
| RingSmoothingForce | ✓ | — | — | — | — |
| ECMConfinementForce | ✓ | ✓ | — | — | — |
| DynamicECMField | ✓ | ✓ | — | — | ✓ |
| LumenPressureForce | ✓ | ✓ | ✓ | — | ✓ |
| FastNagaiHondaForce | — | ✓ | — | — | — |
| DifferentialAdhesionForce | ✓ | — | ✓ | — | ✓ |
| ApicalConstrictionForce | ✓ | ✓ | ✓ | — | ✓ |
| CellPolarityForce | ✓ | — | ✓ | — | — |
| **Forces (3D)** |
| SurfaceSpringForce | — | — | ✓ | — | — |
| SurfaceSmoothingForce | — | — | ✓ | — | — |
| ECMConfinementForce3d | — | — | ✓ | ✓ | — |
| DynamicECMField3d | — | — | ✓ | ✓ | ✓ |
| DynamicECMContactGuidanceForce3d | — | — | ✓ | ✓ | — |
| BasementMembraneForce | — | — | ✓ | ✓ | ✓ |
| **Topology** |
| RingTopologyTracker | ✓ | — | — | — | — |
| SurfaceTopologyTracker | — | — | ✓ | — | — |
| **Writers** |
| ECMFieldWriter / ECMFieldWriter3d | ✓ | ✓ | ✓ | ✓ | ✓ |
| RingOutlineWriter | ✓ | — | — | — | — |
| CellPolarityWriter | ✓ | — | ✓ | — | — |
| CellContactInhibitionStatusWriter | ✓ | ✓ | ✓ | ✓ | — |
| SurfaceMeshWriter | — | — | ✓ | — | — |
| **Modifiers** |
| ContinuousPvdModifier | ✓ | ✓ | ✓ | — | — |
| CryptBuddingSummaryModifier | ✓ | ✓ | ✓ | ✓ | — |
| **Division Rules** |
| TangentialCentreBasedDivisionRule | ✓ | — | ✓ | — | — |
| LocalTangentVertexBasedDivisionRule | — | ✓ | — | — | — |
| **Utilities** |
| TimedForce | — | ✓ | ✓ | ✓ | — |
| SimProfiler | ✓ | ✓ | ✓ | ✓ | — |

---

## Naming Conventions

| Pattern | Meaning |
|---------|---------|
| `Run{Node,Vertex}{2d,3d}.hpp` | Runner for a specific model type |
| `*Force.hpp` | AbstractForce subclass |
| `*Writer.hpp` | Cell/population writer for VTK output |
| `*Modifier.hpp` | Simulation modifier (runs each timestep) |
| `*Tracker.hpp` | Topology data structure |
| `*CellCycleModel.hpp` | AbstractCellCycleModel subclass |
| `*MutationState.{hpp,cpp}` | AbstractCellMutationState subclass |
| `*DivisionRule.hpp` | Division direction rule |

---

## How to Add a New Force/Class

1. Create `.hpp` (and optionally `.cpp`) in `src/`
2. If `.cpp` was added: re-run `cmake .` from build root
3. Include in the relevant `Run*.hpp` runner
4. Wire up to a parameter toggle in `CryptBuddingParams.hpp` if needed
5. Add to the dependency table above

## How to Add a New Cell Type / Mutation State

Follow the pattern of `TACellMutationState`:
1. Create `NewTypeMutationState.hpp` + `.cpp` in `src/`
2. Inherit from `AbstractCellMutationState`, pass unique colour index
3. Add Chaste serialization boilerplate
4. Re-run `cmake .`
5. Include and `MAKE_PTR()` in the runner(s)

---

## Test Suite

| Test | Description |
|------|-------------|
| `TestCryptBudding` | Unified integration test for all 4 model types (env vars: `ECM_STIFFNESS`, `RUN_NUMBER`) |
| `Test2dCryptBuddingNodeBased` | 2D node ring (80 cells) |
| `Test2dCryptBuddingVertexBased` | 2D vertex annulus (40 elements) |
| `Test3dCryptOrganoid` | 3D node sphere (100 cells) |
| `Test3dVertexCryptOrganoid` | 3D vertex monolayer (OrganoidChaste) |
| `Test2dApicalConstriction` | Apical constriction validation |
| `Test2dDynamicECMInvasion` | Invasive migration |
| `TestHello_TissueMorphology` | Sanity check |

---

## Output Formats

- **VTK/ParaView** (`.vtu`, `.pvtu`, `.pvd`) — cell positions, forces, time series
- **CSV** — summary statistics via `CryptBuddingSummaryModifier`

---

## Gotchas & Common Pitfalls

### Force wiring order in Run*.hpp

The `if/else` chains in each runner are **mutually exclusive**. A force enabled by an earlier branch shadows all later branches. For example, `enableDifferentialAdhesion` is checked before `useTopologyBasedSprings` — if both are true, only the first branch fires. When adding new force options, be aware of the priority order and ensure combinations work correctly.

### Feature flag interactions (non-obvious)

| Flags | Behaviour |
|-------|-----------|
| `enableGhostNodeECM=true` + `enableViscoelasticECM=true` | Uses `ViscoelasticGhostNodeEcmForce/Field` |
| `enableGhostNodeECM=true` + `enableViscoelasticECM=false` | Uses elastic `GhostNodeEcmForce/Field` |
| `enableGhostNodeECM=false` + `enableEcmConfinement=true` | Uses grid-based `ECMConfinementForce` |
| `enableDifferentialAdhesion=true` + `useTopologyBasedSprings=true` | DifferentialAdhesionForce with topology filtering (inherits RingSpringForce) |
| `enableDifferentialAdhesion=true` + `useTopologyBasedSprings=false` | DifferentialAdhesionForce with distance-based springs |

### Cell data keys (string-based, no compile-time checking)

These keys are written/read via `GetCellData()->SetItem(key, value)`. Typos are silent bugs.

| Key | Set by | Used by |
|-----|--------|---------|
| `"polarity_theta"`, `"polarity_phi"` | Runner (init), CellPolarityForce | TangentialCentreBasedDivisionRule, CellPolarityWriter |
| `"is_apical"` | Runner (init) | DifferentialAdhesionForce |
| `"cell_type_id"` | Runner (init), StochasticFourTypeCellCycleModel | Writers, summary modifier |
| `"volume"` | Runner (init), population | ContactInhibition cell cycle models |

### Runner code duplication

The 4 runners (`RunNode2d/Vertex2d/Node3d/Vertex3d.hpp`) share ~700 lines of duplicated logic (cell cycle setup, cell type assignment, ECM force wiring, relaxation phase). Changes to shared logic must be replicated in all 4 files. This is a known tech-debt item.

### Serialization

- **Header-only classes** (forces, modifiers, writers, trackers): no `.cpp` needed, no `cmake .` re-run
- **Mutation states**: require `.cpp` with `CHASTE_CLASS_EXPORT()` macro + `cmake .` re-run
- Forces are **not checkpointed** — simulations cannot resume from mid-run. Fine for one-shot runs.

### Parameter auto-derivation in Finalise()

`Finalise()` overwrites `dt`, `endTime`, `t1Threshold`, `t2Threshold`, and radii **unless** the corresponding `*Overridden` flag is true (set automatically when JSON provides an explicit value). If you set `dt` in JSON, the auto-tuning is disabled for `dt`. This is silent — no warning is printed.

---

## HPC (BluePebble, University of Bristol)

- Singularity containers for environment encapsulation
- Slurm job arrays for parameter sweeps
- Scripts in `hpc/`; see `hpc/README.md`

---

## Documentation

| File | Contents |
|------|---------|
| `README.md` | Architecture, HPC pipeline, build & run |
| `apps/README.md` | JSON parameter → C++ class mapping |
| `apps/PARAMS_README.md` | Full parameter reference with units/ranges |
| `test/README.md` | Test suite documentation |
| `ECM_Forces_Comparison.md` | Mathematical summary of ECM force models |
| `src/forces/ViscoelasticECM_Methodology.md` | Viscoelastic ECM derivation |
| `docs/ConstitutiveEquations_Forces.md` | Force constitutive equations |
| `docs/ConstitutiveEquations_CellCycles.md` | Cell cycle equations |
| `experiments/CryptBudding/README.md` | Crypt budding experiment setup |
| `experiments/2dECMStiffness/README.md` | 2D stiffness sweep setup |
