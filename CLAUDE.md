# TissueMorphology — AI Assistant Context

Chaste-based C++ project simulating intestinal crypt budding under varying ECM properties.
Masters thesis (Eng. Maths, Bristol) investigating hydrogel stiffness effects on organoid morphogenesis.

**Paths:** project=`/home/orlando/Thesis/Chaste/projects/TissueMorphology`, build=`/home/orlando/Thesis`

## Build & Run

```bash
cd /home/orlando/Thesis
cmake .                          # after adding new .cpp files only
make -j$(nproc) CryptBuddingApp # build
./projects/TissueMorphology/apps/CryptBuddingApp \
  -config Chaste/projects/TissueMorphology/apps/params-Vertex2d.json \
  -model vertex2d -shear-modulus 1300 -run 0
```

## Architecture

```
CryptBuddingApp.cpp → CryptBuddingParams.hpp (JSON→struct→Finalise()) → Run{Node,Vertex}{2d,3d}.hpp
```

Key app source files in `apps/src/`:
- `CryptBuddingParams.hpp` — all parameters, JSON parsing, unit conversion, auto-tuning
- `ECMWiringHelpers.hpp` — shared ECM force setup logic
- `CellSetupHelpers.hpp` — shared cell cycle/type setup logic
- `SimulationHelpers.hpp` — shared simulation creation logic

## Parameter naming (current)

| Parameter | Purpose |
|-----------|---------|
| `cellGhostSpringStiffness` | Cell-ghost node spring stiffness (formerly `ecmStiffness`) |
| `ghostE0` | SLS relaxed modulus E₀ (formerly `ghostRelaxedStiffness`) |
| `ghostE1` | SLS transient modulus E₁ (formerly `ghostRelaxationModulus`) |
| `ecmShearModulusPa` | Physical shear modulus (Pa); auto-derives ghostE0/E1 via $G_u \approx (E_0+E_1)/(2.9a)$ |

## Feature flag interactions

| Flags | Result |
|-------|--------|
| `enableGhostNodeECM` + `enableViscoelasticECM` | `ViscoelasticGhostNodeEcmForce/Field` (SLS model) |
| `enableGhostNodeECM` only | Elastic `GhostNodeEcmForce/Field` |
| `enableEcmConfinement` (no ghost) | Grid-based `ECMConfinementForce` |

## Gotchas

- **cmake re-run:** Required after adding `.cpp` files (file globbing). Not needed for `.hpp`.
- **Runner duplication:** 4 runners share ~700 lines; changes must be replicated. Shared logic is being extracted to `ECMWiringHelpers.hpp`, `CellSetupHelpers.hpp`, `SimulationHelpers.hpp`.
- **Finalise() auto-tuning:** Silently overwrites `dt`, `endTime`, `t1/t2Threshold`, radii unless JSON explicitly sets them.
- **Cell data keys:** String-based (`"cell_type_id"`, `"volume"`, `"polarity_theta"`). Typos are silent bugs.
- **Force wiring order:** `if/else` chains in runners are mutually exclusive. Earlier branches shadow later ones.
- **Serialization:** Forces not checkpointed; simulations cannot resume mid-run.

## Documentation map

| Path | Contents |
|------|---------|
| `docs/quickstart/` | Build, run, configure |
| `docs/equations/` | Constitutive equations for all forces and cell cycles |
| `docs/units/` | Unit conversion system (CHASTE→SI) |
| `docs/ARCHITECTURE.md` | Dependency matrix, naming conventions, recipes |
| `docs/PARAMETER_JUSTIFICATION.md` | Literature-sourced parameter validation |
| `docs/dev_history/` | Archived historical design documents |
| `apps/README.md` | JSON parameter → C++ class mapping |
| `test/README.md` | Test suite documentation |
