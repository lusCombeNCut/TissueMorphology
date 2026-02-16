# Experiment Design: 3D Organoid Formation and Crypt Morphogenesis

## 1. Overview

This test suite simulates **3D intestinal organoid development** — from an
initial spherical cyst of epithelial cells to a budding structure with crypt-like
protrusions. Two complementary cell-based frameworks are used:

| Test | Framework | Cell representation | Key mechanics |
|------|-----------|-------------------|---------------|
| `Test3dOrganoidFormation` | Node-based (overlapping spheres) | Point particles + springs | `BasementMembraneForce` confinement |
| `Test3dCryptOrganoid` | Node-based + dynamic ECM | Point particles + springs | BM + 3D ECM fiber guidance + cell type heterogeneity |
| `Test3dVertexCryptOrganoid` | Vertex-based (finite-thickness monolayer) | 3D polyhedral cells | `SurfaceTensionForce` + BM + ECM + epithelial buckling |

The vertex model is the most biophysically detailed, explicitly tracking cell
shapes, apical/basal surfaces, and T1 transitions required for crypt invagination.

---

## 2. Geometry

### Initial Configuration

All tests start from a **spherical shell** (or solid sphere) of cells,
representing the cystic stage of organoid development (Sato et al., 2009).

| Test | Shape | Radius | Cells | Notes |
|------|-------|--------|-------|-------|
| `Test3dOrganoidFormation` (Basic) | Solid sphere | 20.0 µm | 30 | Random in volume |
| `Test3dOrganoidFormation` (Stiffness/LongTerm) | Solid sphere | 18.0 µm | 30–40 | Random in volume |
| `Test3dCryptOrganoid` | Spherical shell | 25.0 µm | 100 | Fibonacci-distributed on shell, 3 µm thickness |
| `Test3dVertexCryptOrganoid` | Spherical monolayer | 10.0 µm inner | 100–200 | Finite-thickness vertex mesh |

**Shell vs solid:**
- `Test3dOrganoidFormation` uses $r = R \cdot U^{1/3}$ sampling (uniform in volume)
- `Test3dCryptOrganoid` uses Fibonacci sphere lattice with random radial perturbation in $[R - h/2, R + h/2]$
- `Test3dVertexCryptOrganoid` uses `FiniteThicknessRandomizedSphereMeshGenerator` from OrganoidChaste

---

## 3. Constitutive Equations

### 3.1 Cell–Cell Springs (Node-Based Models)

Generalised linear spring force (Meineke et al., 2001):

$$\mathbf{F}_{ij} = k_s \cdot (|\mathbf{r}_{ij}| - s_{ij}) \cdot \hat{\mathbf{r}}_{ij}$$

| Test | $k_s$ | Cutoff | Division rest length | Growth duration |
|------|-------|--------|---------------------|----------------|
| Basic | 30.0 | default | 0.5 | 1.0 h |
| Stiffness | 25.0 | default | — | — |
| LongTerm | 35.0 | default | — | — |
| CryptOrganoid | 20.0 | 15.0 µm | 0.5 | 1.0 h |

### 3.2 Surface Tension Force (Vertex Model)

The finite-thickness vertex model from OrganoidChaste (Drozdowski & Schwarz,
2025) computes forces from surface tensions on apical, basal, and lateral cell
faces:

$$E = \sum_\alpha \left[ \gamma_A S_A^\alpha + \gamma_B S_B^\alpha + \gamma_L S_L^\alpha + K_V (V^\alpha - V_0^\alpha)^2 \right]$$

where:
- $\gamma_A, \gamma_B, \gamma_L$ are apical, basal, and lateral surface tensions
- $S_A^\alpha, S_B^\alpha, S_L^\alpha$ are face areas
- $K_V$ enforces volume homeostasis via `GeometricalTargetVolumeModifier`

Forces on each vertex: $\mathbf{F}_i = -\nabla_i E$.

**Cell height** is derived from equilibrium cell shape:

$$h = \frac{2}{3\sqrt{3}} \left(\frac{9}{2}\right)^{2/3} \left(\frac{\gamma_A + \gamma_B}{\gamma_L}\right)^{2/3}$$

**T1 transition threshold:**

$$l_{T1} = \frac{0.66}{\left(3^2 \cdot \frac{1 + \gamma_L}{\gamma_L}\right)^{1/3}}$$

### 3.3 Basement Membrane Force

Radial restoring force confining cells within a target radius (Dunn et al., 2012):

$$\mathbf{F}_\text{BM} = -k_\text{BM} \cdot (r - R_\text{BM}) \cdot \hat{\mathbf{r}}, \quad r > R_\text{BM}$$

With optional time-dependent ECM degradation:

$$R_\text{BM}(t) = R_\text{BM}^0 + \dot{R} \cdot t, \quad R_\text{BM}(t) \leq R_\text{max}$$

### 3.4 3D ECM Contact Guidance Force (CryptOrganoid and VertexCryptOrganoid)

Migration bias from a 3D fiber field stored on a voxel grid (Painter, 2009;
Metzcar et al., 2025):

$$\mathbf{F}_\text{ECM} = v_0 \cdot s \cdot \sqrt{\rho} \cdot \mathbf{n}_\text{migration}$$

where $\mathbf{n}_\text{migration}$ combines fiber-aligned, perpendicular, and
random components weighted by fiber density. The ECM field uses **radial**
fiber orientation (fibers pointing outward from the organoid centre), mimicking
collagen/Matrigel gel surrounding an organoid.

**ECM dynamics:**

| Process | Equation | Rate |
|---------|----------|------|
| Degradation | $\partial\rho/\partial t = -k_d \rho$ | $k_d = 0.002$ |
| Remodeling | $\partial\phi/\partial t = k_r(\theta_\text{cell} - \phi)$ | $k_r = 0.05$ |
| Deposition | $\partial\rho/\partial t = +k_p$ | $k_p = 0.0003$ |

### 3.5 Volume Homeostasis

- **Node-based**: `VolumeTrackingModifier` + `ContactInhibitionCellCycleModel`
- **Vertex-based**: `GeometricalTargetVolumeModifier` with reference volume computed from sphere shell geometry:

$$V_\text{avg} = \frac{4\pi}{3N}(R_\text{outer}^3 - R_\text{inner}^3)$$

---

## 4. Cell Type Assignment

### Test3dCryptOrganoid

Cells are assigned types based on normalised $z$-coordinate ($z_\text{frac} = z / R$):

| Region | $z_\text{frac}$ | Cell type | Proliferative | ECM degradation factor |
|--------|-----------------|-----------|---------------|----------------------|
| Crypt base | $< -0.5$ | Stem (+Paneth at 1/3) | Yes | 2.0 (high) |
| Transit zone | $[-0.5, 0.3)$ | Transit-amplifying | Yes | 1.0 (baseline) |
| Villus-like | $\geq 0.3$ | Differentiated | No | 0.5 (low) |

This mirrors the crypt–villus axis (Barker et al., 2007).

### Test3dVertexCryptOrganoid

Same spatial logic, with cell types assigned by element centroid $z$:

| Region | $z_\text{frac}$ | Cell type | G1 duration |
|--------|-----------------|-----------|-------------|
| Crypt base | $< -0.5$ | Stem | 25.0 h |
| Transit zone | $[-0.5, 0.3)$ | Transit-amplifying | 15.0 h |
| Villus-like | $\geq 0.3$ | Differentiated | $10^6$ h (no division) |

### Test3dOrganoidFormation

Uses `OrganoidCellFactory` with core/surface heterogeneity:
- Core ($r < 30$): type 0, BM stiffness 1.0
- Surface ($r \geq 30$): type 1, BM stiffness 5.0

---

## 5. Simulation Protocol

### Two-Phase Approach (Vertex Model)

Standard practice for vertex simulations (Osborne et al., 2017):

**Phase 1 — Relaxation** (5 time units):
- No cell division
- Surface tension forces active
- BM confinement without degradation
- Simulated annealing parameters suppress T1 transitions
- Fine timestep ($dt = 0.001$)

**Phase 2 — Growth** (100 time units):
- Cell division enabled (stem + TA cells)
- BM degradation activated
- ECM guidance force added
- Active T1 transitions enabled
- Coarser timestep ($dt = 0.006$)

### Single-Phase (Node-Based)

`Test3dOrganoidFormation`: no explicit relaxation (short run, $t = 2$–100 h).
`Test3dCryptOrganoid`: single continuous run for 168 h (7 days) with all
forces active from the start.

---

## 6. Modelling Choices

### Why two frameworks?

| Feature | Node-based | Vertex-based |
|---------|-----------|-------------|
| Cell shape | Implicit (spheres) | Explicit (polyhedra) |
| Epithelial buckling | Emergent from spring forces | Natural from surface tension |
| T1 transitions | N/A | Explicit neighbour exchanges |
| Computational cost | Low | High |
| Crypt detection | Surface curvature (future) | Direct from mesh geometry |

The node-based model tests whether proliferation + BM confinement + ECM guidance
alone produce crypt-like protrusions. The vertex model provides explicit cell
shape mechanics needed for realistic buckling instability (Hannezo et al., 2011).

### Why radial ECM?

Intestinal organoids are cultured in Matrigel gel where collagen fibers extend
radially from the organoid surface. Radial ECM bias guides budding protrusions
outward, consistent with experimental observations (Gjorevski et al., 2016).

### Cell number choices

- 30–40 cells (OrganoidFormation): minimal tests for force validation
- 100 cells (CryptOrganoid): sufficient for meaningful heterogeneity and statistics
- 100–200 cells (VertexCryptOrganoid): realistic for a small organoid cyst

---

## 7. Output and Post-Processing

### Output files

- `results_*.vtu`: cell positions and properties (ParaView)
- `ecm3d_*.vti`: 3D ECM field voxel data (ParaView)
- `ecm3d_results.pvd`: ECM time series index

### Dependent variables

- **Cell count** over time (growth curve)
- **Cell type composition** (stem / TA / differentiated ratios)
- **Maximum radial distance** (organoid extent / expansion ratio)
- **Crypt detection** — planned via surface curvature analysis (not yet implemented)

---

## 8. References

- Sato, T. et al. (2009) *Nature* 459:262–265
- Barker, N. et al. (2007) *Nature* 449:1003–1007
- Painter, K.J. (2009) *J. Math. Biol.* 58:511–543
- Metzcar, J. et al. (2025) — ECM contact guidance model
- Drozdowski, O. & Schwarz, U.S. (2025) — OrganoidChaste vertex model
- Dunn, S.-J. et al. (2012) *J. Theor. Biol.* 298:82–91
- Hannezo, E. et al. (2011) *Phys. Rev. Lett.* 107:078104
- Gjorevski, N. et al. (2016) *Nature* 539:560–564
- Osborne, J.M. et al. (2017) *PLoS Comput. Biol.* 13(7):e1005387
- Meineke, F.A. et al. (2001) *Cell Prolif.* 34(4):253–266
- Nagai, T. & Honda, H. (2001) *Phil. Mag. B* 81(7):699–719
