# Experiment Design: ECM Force Validation Tests

## 1. Overview

This test suite validates two fundamental force components used across the
TissueMorphology project: **apical constriction** and **basement membrane
confinement with ECM degradation**. These are unit-/integration-level tests
that verify individual force mechanics before they are combined in full
organoid or invasion simulations.

| Test | Dimension | Framework | Purpose |
|------|-----------|-----------|---------|
| `Test2dApicalConstriction` | 2D | Node-based (overlapping spheres) | Verify apical wedging drives invagination in a flat sheet |
| `Test2dThresholdECMForce` | 2D | Node-based | Validate BM confinement + ECM degradation in a circular organoid |
| `TestThresholdECMForce` | 3D | Node-based | Same as above in 3D spherical geometry |

---

## 2. Test 1: Apical Constriction (`Test2dApicalConstriction`)

### Geometry

A flat 10×10 grid of cells at unit spacing (100 cells total). Cells within a
2-unit radius of the grid centre are marked `is_apical = 1.0`; all others are
non-apical.

### Constitutive Equations

#### Cell–Cell Springs

Standard generalised linear spring force (Meineke et al., 2001):

$$\mathbf{F}_{ij} = k_s \cdot (|\mathbf{r}_{ij}| - s_{ij}) \cdot \hat{\mathbf{r}}_{ij}$$

with $k_s = 30.0$ (force/length, dimensionless Chaste units).

#### Apical Constriction Force

For cells with `is_apical = 1.0`, an inward radial force drives area reduction
(Martin et al., 2009):

$$\mathbf{F}_\text{apical} = -\sigma_\text{ac} \cdot (A - A_\text{target}) \cdot \hat{\mathbf{r}}$$

where:
- $\sigma_\text{ac} = 10.0$ — constriction strength
- $A_\text{target} = A_\text{ref} \cdot (1 - f_\text{reduction}),\quad f_\text{reduction} = 0.6$
- Force acts only when $A > A_\text{target}$ (cell needs to constrict)

### Modelling Choices

- **Flat sheet geometry** deliberately tests whether the constriction force
  alone generates an invagination pit from a planar starting configuration.
- **No basement membrane or lumen pressure** — isolates the apical constriction
  effect.
- Birth time set to $-10$ h so cells are mature and non-dividing during the
  10 h simulation window.

### Expected Outcome

Central apical cells should form a visible invagination pit as their area
contracts, pulling neighbouring cells inward.

**References:**
- Martin, A.C. et al. (2009) *Nature* 457:495–499
- Meineke, F.A. et al. (2001) *Cell Prolif.* 34(4):253–266

---

## 3. Test 2: 2D Threshold ECM Force (`Test2dThresholdECMForce`)

### Geometry

50 cells distributed uniformly within a disk of radius 5.0 µm (using
$r = R \sqrt{U}$ for uniform-in-area sampling). This models a 2D
cross-section of an organoid.

### Cell Heterogeneity

Cells are classified by distance from the origin:

| Region | Distance | Type | BM stiffness |
|--------|----------|------|-------------|
| Inner core | $r < 3.0$ | Core (type 0) | 1.0 |
| Outer surface | $r \geq 3.0$ | Surface (type 1) | 5.0 |

The `OrganoidCellFactory` assigns proliferative types (25% stem probability).

### Constitutive Equations

#### Cell–Cell Springs

$$\mathbf{F}_{ij} = k_s \cdot (|\mathbf{r}_{ij}| - s_{ij}) \cdot \hat{\mathbf{r}}_{ij}, \quad k_s = 35.0$$

#### Basement Membrane Force

A radial restoring force confining cells within a target radius $R_\text{BM}$
(Dunn et al., 2012):

$$\mathbf{F}_\text{BM} = -k_\text{BM} \cdot (r - R_\text{BM}) \cdot \hat{\mathbf{r}}, \quad \text{for } r > R_\text{BM}$$

with $k_\text{BM} = 2.0$ and initial $R_\text{BM} = 6.0$ µm.

#### ECM Degradation

Time-dependent radius expansion modelling MMP-mediated BM degradation
(Alcaraz et al., 2011):

$$R_\text{BM}(t) = R_\text{BM}^0 + \dot{R} \cdot t, \quad R_\text{BM}(t) \leq R_\text{max}}$$

with $\dot{R} = 0.4$ µm/time-unit and $R_\text{max} = 20.0$ µm.

### Modelling Choices

- **Heterogeneous BM stiffness** per cell (stored in `CellData`) allows surface
  cells to experience stronger confinement than core cells.
- **Long simulation** (20 time units) allows observing expansion dynamics.
- Small dt (0.005) with high movement threshold (50.0) ensures stability during
  initial settling.

**References:**
- Dunn, S.-J. et al. (2012) *J. Theor. Biol.* 298:82–91
- Alcaraz, J. et al. (2011) *EMBO J.* 30(11):2284–2295

---

## 4. Test 3: 3D Threshold ECM Force (`TestThresholdECMForce`)

### Geometry

50 cells distributed uniformly within a sphere of radius 18.0 µm (using
$r = R \cdot U^{1/3}$ for uniform-in-volume sampling, with spherical
coordinates).

### Cell Heterogeneity

| Region | Distance | Type | BM stiffness |
|--------|----------|------|-------------|
| Inner core | $r < 12.0$ | Core (type 0) | 1.0 |
| Outer surface | $r \geq 12.0$ | Surface (type 1) | 5.0 |

### Constitutive Equations

Identical to the 2D version (§3) but in 3D, with:
- $k_s = 35.0$, interaction cutoff 5.0 µm
- $k_\text{BM} = 2.0$, $R_\text{BM}^0 = 20.0$ µm
- $\dot{R} = 0.4$ µm/time-unit, $R_\text{max} = 80.0$ µm

### Modelling Choices

- **3D spherical geometry** validates that the `BasementMembraneForce` correctly
  computes radial distances and forces in three dimensions.
- Higher movement threshold (500.0) accommodates initial 3D settling dynamics.
- Uses PETSc parallelisation for performance.

**References:**
- Same as §3.
