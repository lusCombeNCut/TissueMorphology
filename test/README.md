# Test Suite Overview: CryptBudding

## 0. Unified Test — `TestCryptBudding.hpp`

A single entry-point test (`CryptBudding/TestCryptBudding.hpp`) drives all four
model configurations via environment variables:

```bash
# Select model and run
MODEL_TYPE=node2d   ECM_STIFFNESS=5.0 RUN_NUMBER=0 ctest -R TestCryptBudding
MODEL_TYPE=vertex2d ECM_STIFFNESS=5.0 RUN_NUMBER=0 ctest -R TestCryptBudding
MODEL_TYPE=node3d   ECM_STIFFNESS=5.0 RUN_NUMBER=0 ctest -R TestCryptBudding
MODEL_TYPE=vertex3d ECM_STIFFNESS=5.0 RUN_NUMBER=0 ctest -R TestCryptBudding
```

**Feature toggles** (0 = off, 1 = on; defaults shown):

| Variable | Default | Applies to |
|----------|---------|------------|
| `ENABLE_LUMEN_PRESSURE` | 1 | All models |
| `ENABLE_APICAL_CONSTRICTION` | 1 | All models |
| `ENABLE_ECM_GUIDANCE` | 0 | 3D only |
| `ENABLE_RELAXATION` | 1 | All models |
| `ENABLE_SLOUGHING` | 1 | 2D only |
| `ENABLE_DIFFERENTIAL_ADHESION` | 1 | Node-based only |

All four models share:
- `ContactInhibitionCellCycleModel` (density-dependent proliferation)
- Three-tier cell type assignment (Stem / TA / Differentiated)
- `BasementMembraneForce` with ECM degradation
- `VolumeTrackingModifier`
- Relaxation phase (geometry settling before proliferation)
- Summary CSV writer and progress logging

The individual per-model test files (`CryptBudding2d/`, `CryptBudding3d/`) remain
for reference and backwards compatibility.

---

## 1. Initial Conditions

### CryptBudding2d

| Test | Geometry | Cell placement | Num cells | Radii | Symmetry breaking |
|------|----------|---------------|-----------|-------|-------------------|
| **Test2dCryptBuddingNodeBased** | Annular ring of point particles on a circle | Single ring at radius $R_0 = 8.0$ | 80 | $R_0 = 8.0$ | Radial noise $\delta r \sim \mathcal{U}(-0.15, 0.15)$ |
| **Test2dCryptBuddingVertexBased** | Annular ring of quadrilateral vertex elements | Two rings of nodes (inner + outer) forming wedge quads | 40 elements (80 nodes) | $R_\text{inner} = 6.0$, $R_\text{outer} = 8.0$ | Radial noise $\delta r \sim \mathcal{U}(-0.1, 0.1)$ on each node |

Both 2D tests assign cell types by angular position on the ring (bottom = stem, flanks = TA, top = differentiated).

### CryptBudding3d

| Test | Geometry | Cell placement | Num cells | Radii | Symmetry breaking |
|------|----------|---------------|-----------|-------|-------------------|
| **Test3dCryptOrganoid** (node-based) | Spherical shell of point particles | Fibonacci sphere at radius $R = 25\,\mu m$ with shell thickness $3\,\mu m$ | 100 | $R = 25$ | Random perturbation within shell thickness |
| **Test3dVertexCryptOrganoid** (vertex) | Spherical monolayer vertex mesh (OrganoidChaste finite-thickness model) | `FiniteThicknessRandomizedSphereMeshGenerator` | 100 (relax) / 200 (growth) | $R_\text{inner} = 10$ | Built-in mesh randomisation |
| **Test3dOrganoidFormation** (node-based) | Solid sphere of random points | Random $(r, \theta, \phi)$ uniformly in volume | 30–40 | $R = 18$–$20$ | Intrinsically random |

**Key difference:** The 3D tests use z-coordinate fraction ($z / R$) for cell type assignment, whereas the 2D tests use angular fraction on the ring.

---

## 2. Cell Cycle Models

All primary CryptBudding tests (and the unified `TestCryptBudding.hpp`) now use
`ContactInhibitionCellCycleModel` for consistency. The older standalone tests
(`Test3dVertexCryptOrganoid`, `Test3dOrganoidFormation`) still use their original
cell cycle models.

| Test file | Cell cycle model | Key parameters | Density feedback? |
|-----------|-----------------|----------------|-------------------|
| **TestCryptBudding** (all models) | `ContactInhibitionCellCycleModel` | `quiescentVolumeFraction = 0.7`, `equilibriumVolume` varies by model | **Yes** |
| **Test2dCryptBuddingNodeBased** | `ContactInhibitionCellCycleModel` | `quiescentVolumeFraction = 0.7`, `equilibriumVolume = 1.0` | **Yes** |
| **Test2dCryptBuddingVertexBased** | `ContactInhibitionCellCycleModel` | `quiescentVolumeFraction = 0.6`, `equilibriumVolume = target_area` | **Yes** |
| **Test3dCryptOrganoid** | `ContactInhibitionCellCycleModel` | `quiescentVolumeFraction = 0.7`, `equilibriumVolume = 1.0` | **Yes** |
| **Test3dVertexCryptOrganoid** (growth) | `FixedG1GenerationalCellCycleModel` | `stemG1 = 25h`, `transitG1 = 15h` | **No** |
| **Test3dOrganoidFormation** | `UniformG1GenerationalCellCycleModel` | Default durations | **No** |

### Where are the cell cycle models defined?

All cell cycle models live in Chaste core at:

- `ContactInhibitionCellCycleModel` → [cell_based/src/cell/cycle/ContactInhibitionCellCycleModel.hpp](../../cell_based/src/cell/cycle/ContactInhibitionCellCycleModel.hpp)
  - Inherits `AbstractSimplePhaseBasedCellCycleModel`. Extends G1 when cell is compressed below `quiescentVolumeFraction × equilibriumVolume`. Requires `VolumeTrackingModifier`.
- `FixedG1GenerationalCellCycleModel` → [cell_based/src/cell/cycle/FixedG1GenerationalCellCycleModel.hpp](../../cell_based/src/cell/cycle/FixedG1GenerationalCellCycleModel.hpp)
  - Inherits `AbstractSimpleGenerationalCellCycleModel`. Deterministic G1 duration for stem/transit cells. Uses generation counting (`maxTransitGenerations`).
- `UniformG1GenerationalCellCycleModel` → [cell_based/src/cell/cycle/UniformG1GenerationalCellCycleModel.hpp](../../cell_based/src/cell/cycle/UniformG1GenerationalCellCycleModel.hpp)
  - Like `FixedG1`, but G1 is drawn from a uniform distribution. Also generational.
- `NoCellCycleModel` — built-in sentinel; cells never divide.

---

## 3. Cell Proliferative Types

All tests draw from the same three Chaste built-in types:

| Type | Class | Proliferative? | Assignment rule |
|------|-------|----------------|-----------------|
| **Stem** | `StemCellProliferativeType` | Yes | 2D: bottom 40% of ring; 3D: $z/R < -0.5$ |
| **Transit-amplifying (TA)** | `TransitCellProliferativeType` | Yes | 2D: flanks 30%; 3D: $-0.5 \leq z/R < 0.3$ |
| **Differentiated** | `DifferentiatedCellProliferativeType` | No | 2D: top 30%; 3D: $z/R \geq 0.3$ |

---

## 4. Feature & Force Comparison Table

### Unified Test (`TestCryptBudding.hpp`)

| Feature / Force | node2d | vertex2d | node3d | vertex3d |
|---|:---:|:---:|:---:|:---:|
| **Population type** | Node-based | Vertex-based | Node-based | Monolayer vertex (OrganoidChaste) |
| **Initial geometry** | Annular ring | Annular ring (wedge quads) | Fibonacci sphere | Spherical monolayer mesh |
| **Cell–cell force** | `DifferentialAdhesionForce` | `NagaiHondaForce` | `DifferentialAdhesionForce` | `SurfaceTensionForce` |
| **Basement membrane** | ✅ | ✅ | ✅ | ✅ |
| **ECM degradation** | ✅ | ✅ | ✅ | ✅ |
| **Lumen pressure** | ✅ `LumenPressureForce` | ✅ `LumenPressureForce` | ✅ `LumenPressureForce` | ✅ `LumenPressureSubForce` (OrganoidChaste) |
| **Apical constriction** | ✅ | ✅ | ✅ | ✅ |
| **ECM contact guidance** | ❌ (2D) | ❌ (2D) | ✅ (toggle) | ✅ (toggle) |
| **Cell cycle** | `ContactInhibitionCellCycleModel` | `ContactInhibitionCellCycleModel` | `ContactInhibitionCellCycleModel` | `ContactInhibitionCellCycleModel` |
| **Relaxation phase** | ✅ | ✅ | ✅ | ✅ |
| **Cell sloughing** | ✅ (toggle) | ✅ (toggle) | ❌ | ❌ |
| **Volume tracking** | ✅ | ✅ | ✅ | ✅ |
| **Target area/volume** | ❌ | ✅ `SimpleTargetAreaModifier` | ❌ | ✅ `GeometricalTargetVolumeModifier` |
| **Simulation class** | `OffLatticeSimulation<2>` | `OffLatticeSimulation<2>` | `OffLatticeSimulation<3>` | `FiniteThicknessSimulation3d` |

### Standalone Tests (original, retained for reference)

| Feature / Force | 2D Node | 2D Vertex | 3D Node | 3D Vertex | 3D Formation |
|---|:---:|:---:|:---:|:---:|:---:|
| **Cell cycle** | `ContactInhibition` | `ContactInhibition` | `ContactInhibition` | `FixedG1Generational` | `UniformG1Generational` |
| **Lumen pressure** | ✅ | ✅ | ❌ | ❌ | ❌ |
| **Apical constriction** | ✅ | ✅ | ❌ | ❌ | ❌ |
| **Relaxation phase** | ✅ | ❌ | ❌ | ✅ | ❌ |
| **ECM contact guidance** | ❌ | ❌ | ✅ | ✅ | ❌ |
| **Differential adhesion** | ✅ | ❌ | ❌ | ❌ | ❌ |
| **Cell sloughing** | ✅ | ✅ | ❌ | ❌ | ❌ |

---

## 5. NagaiHondaForce vs SurfaceTensionForce

These are **not interchangeable** — they operate on different mesh/population types:

| | `NagaiHondaForce` | `SurfaceTensionForce` |
|---|---|---|
| **Source** | Chaste core (2D vertex) | OrganoidChaste (3D finite-thickness) |
| **Population** | `VertexBasedCellPopulation<2>` | `MonolayerVertexBasedCellPopulation<3>` |
| **Mesh** | `MutableVertexMesh<2,2>` | `MutableMonolayerVertexMesh<3,3>` |
| **Energy** | Area + perimeter + adhesion ($E = \sum_\alpha \kappa_A (A_\alpha - A_0)^2 + \kappa_P P_\alpha^2 + \sum_{\langle ij \rangle} \gamma_{ij} l_{ij}$) | Per-face surface tensions (apical $\gamma_a$, basal $\gamma_b$, lateral $\gamma_l$) + volume elasticity |
| **Reference** | Nagai & Honda (2001) | Drozdowski & Schwarz (2025) |

Both are energy-based vertex models, but `SurfaceTensionForce` is the 3D
extension for finite-thickness epithelia with distinct apical–basal polarity.

---

## 6. Standardisation Decisions

The unified `TestCryptBudding.hpp` applies the following standardisation:

1. **Cell cycle** → `ContactInhibitionCellCycleModel` everywhere (was `FixedG1GenerationalCellCycleModel` in old 3D vertex test)
2. **Relaxation phase** → enabled for all models (was missing in 2D vertex and 3D node)
3. **Lumen pressure** → available for all models (was 2D-only); 3D vertex uses `LumenPressureSubForce` from OrganoidChaste
4. **Apical constriction** → available for all models (was 2D-only)
5. **ECM stiffness sweep** → parameterised via `ECM_STIFFNESS` for all models
6. **Summary writer** → `CryptBuddingSummaryModifier` for all models
7. **Deleted** → `TestOrganoidFormation.hpp` (2D formation test, superseded)
