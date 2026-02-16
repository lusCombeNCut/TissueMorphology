# Experiment Design: 2D Invasive Front into Structured ECM

## 1. Overview

This test suite reproduces and extends the **invasive cellular front** scenario
from **Metzcar et al. (2025)** and the original computational study by
**Painter (2009)**. Cells migrate from a boundary into a domain filled with
structured ECM fibers, and the fiber orientation pattern controls the invasion
morphology.

Three test classes implement progressively more realistic models:

| Test | ECM model | Cell influx | Proliferation | Key feature |
|------|-----------|-------------|---------------|-------------|
| `Test2dInvasiveFront` | Static (analytical) | No | No | Minimal ECM guidance validation |
| `Test2dDynamicECMInvasion` | Dynamic grid field | No | Contact-inhibited | Cells remodel and degrade ECM |
| `Test2dPainterReplication` | Dynamic grid field | Yes (30 cells / 180 min) | No | Full Painter (2009) replication with cell influx |

All three test four ECM orientation patterns:
1. **Random** — baseline diffusive invasion
2. **Parallel** — fibers horizontal (perpendicular to invasion direction → slow)
3. **Perpendicular** — fibers vertical (aligned with invasion → fast)
4. **Mixed** — perpendicular in centre, parallel on flanks

---

## 2. Geometry

### Domain

Rectangular domain: 600 µm (width) × 1000 µm (height). Cells enter at the
bottom edge ($y \approx 0$) and migrate upward into the ECM.

### Initial Cell Placement

30 cells placed randomly along the bottom edge:
- $x \sim \mathcal{U}(50, 550)$ µm
- $y \sim \mathcal{U}(0, 50)$ µm

### ECM Grid

A discrete field of fiber orientations on a regular grid:
- Grid spacing: 12.5 µm (`Test2dDynamicECMInvasion`) or 25.0 µm (`Test2dPainterReplication`)
- Each grid voxel stores fiber angle $\phi \in [0, 2\pi)$ and density $\rho \in [0, 1]$

---

## 3. Constitutive Equations

### 3.1 ECM Contact Guidance Force

The central force model from Painter (2009). Cell migration velocity is biased
by the local ECM fiber orientation:

$$\mathbf{v}_\text{cell} = v_0 \cdot s \cdot \sqrt{\rho} \cdot \mathbf{n}_\text{migration}$$

where:
- $v_0$ is the base migration speed (µm/min)
- $s$ is the ECM sensitivity
- $\rho$ is the local ECM fiber density
- $\mathbf{n}_\text{migration}$ is the migration direction, composed of:

$$\mathbf{n}_\text{migration} = w_\parallel \hat{\mathbf{f}} + w_\perp \hat{\mathbf{f}}_\perp + w_r \hat{\mathbf{r}}_\text{random}$$

where $\hat{\mathbf{f}}$ is the local fiber direction, $\hat{\mathbf{f}}_\perp$
is perpendicular to it, and the weights depend on the **anisotropy** parameter
$a$:
- $a = 0$: isotropic (random walk)
- $a = 1$: highly anisotropic (strong fiber following)

#### Static vs Dynamic ECM

| Property | `ECMContactGuidanceForce` | `DynamicECMContactGuidanceForce` |
|----------|--------------------------|----------------------------------|
| Fiber orientation | Fixed analytical function | Grid-stored, updated each step |
| Degradation | No | Cells reduce local $\rho$ |
| Remodeling | No | Cells align fibers toward traction direction |
| Deposition | No | Optional; cells deposit new ECM |

**ECM remodeling equation** (Metzcar et al., 2025):

$$\frac{\partial \phi}{\partial t} = k_\text{remodel} \cdot (\theta_\text{cell} - \phi)$$

where $\theta_\text{cell}$ is the cell's migration direction and $k_\text{remodel}$
is the remodeling rate.

**ECM degradation:**

$$\frac{\partial \rho}{\partial t} = -k_\text{degrade} \cdot \rho$$

when a cell occupies the voxel.

### 3.2 Cell–Cell Adhesion

In `Test2dInvasiveFront`, a generalised linear spring force provides:

$$\mathbf{F}_{ij} = k_s \cdot (|\mathbf{r}_{ij}| - s_{ij}) \cdot \hat{\mathbf{r}}_{ij}, \quad k_s = 0.4$$

**Note:** Cell–cell adhesion is **commented out** in `Test2dDynamicECMInvasion`
and `Test2dPainterReplication` to isolate ECM guidance effects. The force
infrastructure remains for future reactivation.

### 3.3 Contact Inhibition Cell Cycle (Dynamic Invasion Only)

Proliferation in `Test2dDynamicECMInvasion` uses `ContactInhibitionCellCycleModel`:

$$\text{Quiescent if } \frac{V_\text{cell}}{V_\text{eq}} < \phi_q$$

with $\phi_q = 0.7$ and $V_\text{eq} = \pi \cdot 5^2 \approx 78.5$ µm².

Cell cycle durations:
- G1 (stem): 8 h, G1 (transit): 12 h
- S: 5 h, G2: 4 h, M: 1 h
- Total: 18–22 h

### 3.4 Cell Influx (Painter Replication Only)

A custom `CellInfluxModifier` adds 30 new cells at the bottom edge every 180
minutes, replicating the constant cell supply in Painter (2009). New cells are
placed at random positions along the bottom 50 µm and assigned non-dividing
cell cycle models ($T_\text{cycle} = 10^6$ min).

---

## 4. ECM Orientation Patterns

| Pattern | Description | Expected invasion morphology |
|---------|-------------|------------------------------|
| **Random** | $\phi \sim \mathcal{U}(0, 2\pi)$ per voxel | Diffusive front, roughly isotropic |
| **Parallel** | $\phi = 0$ (horizontal) everywhere | Invasion suppressed; cells spread laterally |
| **Perpendicular** | $\phi = \pi/2$ (vertical) everywhere | Deep, narrow invasion front |
| **Mixed** | $\phi = \pi/2$ in centre third, $\phi = 0$ on flanks | Finger-like protrusion in centre |

---

## 5. Modelling Choices

### Why node-based (not vertex)?

The invasion front scenario involves cells migrating through open ECM rather
than forming a densely packed epithelium. Node-based models are appropriate
because:
1. Cells are loosely connected (no tight junctions)
2. Individual cell migration is the primary mechansim
3. No need for explicit cell shape tracking

### Why dynamic ECM?

Static ECM (as in `Test2dInvasiveFront`) produces qualitatively correct
patterns but misses the feedback loop where cells remodel their environment.
Dynamic ECM enables:
- **Path clearing**: leading cells degrade ECM, facilitating followers
- **Track formation**: cell traction aligns fibers, creating migration highways
- These are hallmark features of collective cell invasion (Friedl & Alexander, 2011)

### Simulation duration

- `Test2dInvasiveFront`: 1 day (1440 min) — quick validation
- `Test2dDynamicECMInvasion`: 5 days (7200 min) — full invasion dynamics
- `Test2dPainterReplication`: 5 days (7200 min) — matches Painter (2009)

### Multiple replicates

All tests support `RUN_NUMBER` environment variable for parallel batch
execution with different random seeds (seed = base + run × 1000).

---

## 6. Output and Post-Processing

### Output files

- VTU files: cell positions, types, migration directions
- ECM grid snapshots: fiber angles and densities (via `ECMFieldWriter`)

### Dependent variables

- **95th percentile invasion depth** (Painter metric): $y_{95}$ — the $y$-coordinate
  below which 95% of cells lie
- **Invasion front shape**: spatial distribution of cells at final time
- **ECM remodeling pattern**: final fiber alignment map

---

## 7. References

- Painter, K.J. (2009) *J. Math. Biol.* 58:511–543
- Metzcar, J. et al. (2025) — ECM contact guidance and dynamic remodeling
- Friedl, P. & Alexander, S. (2011) *Cell* 147(5):992–1009
- Meineke, F.A. et al. (2001) *Cell Prolif.* 34(4):253–266
