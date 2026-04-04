# TissueMorphology — Architecture Guide

> **Purpose:** Quick-reference for AI coding agents and new contributors.
> Provides the dependency graph, naming conventions, and edit recipes
> needed to modify the CryptBuddingApp without reading every file.

---

## 1. Build System

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
```

**CRITICAL:** After adding new `.cpp` files to `src/`, you MUST re-run `cmake .`
before `make` — the build system uses file globbing at configure time.
Header-only files (`.hpp`) in `src/` do NOT require cmake re-run.

---

## 2. Execution Flow

```
CryptBuddingApp.cpp
  ├── Parses CLI args: -model {node2d|vertex2d|node3d|vertex3d} -config <file.ini>
  ├── CryptBuddingParams.hpp
  │     ├── SetDefaults()        ← compiled-in defaults
  │     ├── LoadFromFile()       ← INI parser (key=value, # comments, [section] headers ignored)
  │     └── Finalise()           ← auto-derives dt, endTime, t1/t2 thresholds, radii
  ├── Dispatches to one of:
  │     ├── RunNode2d.hpp        ← 2D node-based (ring of cells + ECM field)
  │     ├── RunVertex2d.hpp      ← 2D vertex-based (annulus + Nagai-Honda)
  │     ├── RunNode3d.hpp        ← 3D node-based (sphere surface)
  │     └── RunVertex3d.hpp      ← 3D vertex-based (OrganoidChaste monolayer)
  └── PvdFixUtils.hpp            ← Post-simulation PVD file repair
```

---

## 3. Dependency Graph (src/ → Runner)

### Which classes are used where:

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
| OrganoidCellFactory | — | — | — | — | ✓ |
| Hello_TissueMorphology | — | — | — | — | ✓ |

---

## 4. How to Add a New Parameter

1. **Declare** in `apps/src/CryptBuddingParams.hpp`:
   - Add member variable (with comment)
   - Set default in `SetDefaults()`
   - Add `getDouble/getBool/getUnsigned` call in `LoadFromFile()`
   - Add output line in `SaveToFile()`

2. **Use** in the relevant `Run*.hpp` file(s)

3. **Document** in `apps/params-Node2d.ini` and/or `apps/params-Vertex2d.ini`

4. **Update** `apps/config/default_params.ini` if it's a broadly-applicable parameter

**INI parser behaviour:** Missing keys silently keep their `SetDefaults()` value.
No error is thrown for unknown keys either (they're ignored).

---

## 5. How to Add a New Force/Class

1. Create `.hpp` (and optionally `.cpp`) in `src/`
2. If `.cpp` was added: re-run `cmake .` from build root
3. Include in the relevant `Run*.hpp` runner
4. Wire up to a parameter toggle in `CryptBuddingParams.hpp` if needed
5. Add to the dependency table above

---

## 6. How to Add a New Cell Type / Mutation State

Follow the pattern of `TACellMutationState`:
1. Create `NewTypeMutationState.hpp` + `.cpp` in `src/`
2. Inherit from `AbstractCellMutationState`, pass unique colour index
3. Add Chaste serialization boilerplate
4. Re-run `cmake .`
5. Include and `MAKE_PTR()` in the runner(s)

---

## 7. Naming Conventions

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

## 8. Parameter Flow

```
INI file ──► LoadFromFile() ──► struct members ──► Finalise() ──► Run*.hpp
                                                      ▲
                                             auto-derives:
                                             dt, endTime,
                                             t1/t2 thresholds,
                                             radii, ecmDomainHalf
```

Parameters that trigger autotuning in `Finalise()`:
- `dt` — derived from model type + stiffness if not set in INI
- `endTime` — 168h for most models if not set
- `t1Threshold2d` — from ecmConfinementStiffness threshold
- `t2Threshold2d` — 0.05 default
- `ecmConfinementStiffness` — falls back to `cellGhostSpringStiffness` if sentinel
- Radii (`bmRadius2d`, `ecmMaxRadius2d`, etc.) — from fraction × organoidRadius

Setting any of these explicitly in the INI file prevents the override.

---

## 9. Testing

```bash
# Build and run integration test
make -j$(nproc) TestCryptBudding && ctest -R TestCryptBudding -V

# Build app for manual testing
make -j$(nproc) CryptBuddingApp
./projects/TissueMorphology/apps/CryptBuddingApp -model node2d \
    -config Chaste/projects/TissueMorphology/apps/params-Node2d.ini
```

---

## 10. File Locations Quick Reference

| What | Where |
|------|-------|
| Main app | `apps/src/CryptBuddingApp.cpp` |
| Parameter struct | `apps/src/CryptBuddingParams.hpp` |
| 2D runners | `apps/src/RunNode2d.hpp`, `apps/src/RunVertex2d.hpp` |
| 3D runners | `apps/src/RunNode3d.hpp`, `apps/src/RunVertex3d.hpp` |
| INI configs | `apps/params-Node2d.ini`, `apps/params-Vertex2d.ini` |
| Default config | `apps/config/default_params.ini` |
| Config reference | `apps/README.md` |
| Force equations | `docs/ConstitutiveEquations_Forces.md` |
| Cell cycle equations | `docs/ConstitutiveEquations_CellCycles.md` |
| Sweep scripts | `experiments/CryptBudding/submit_sweep.sh` |
| HPC scripts | `hpc/` |
| Visualizations | `visualization/` |
