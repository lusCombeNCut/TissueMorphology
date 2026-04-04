# ECM Force Comparison — CryptBuddingApp

A mathematical summary of the three ECM force models available in the CryptBudding simulation framework.

---

## 1. ECMConfinementForce (Grid-Based Confinement)

**Purpose:** Confine cells within a tissue boundary by applying density-weighted spring forces between each cell and nearby fixed ECM grid elements.

**ECM representation:** Fixed Eulerian grid (`DynamicECMField` in 2D, `DynamicECMField3d` in 3D). Each grid voxel stores density $\rho \in [0,1]$, fibre angle/direction, and anisotropy.

### Force law

For each cell $i$ interacting with every ECM grid element $j$ within cutoff $r_c$:

$$\mathbf{F}_{ij} = k \, \rho_j \, (d_{ij} - s_0) \, \hat{\mathbf{d}}_{ij}$$

| Symbol | Meaning |
|--------|---------|
| $k$ | Spring stiffness (`mConfinementStiffness`, default 5.0) |
| $\rho_j$ | ECM density at grid element $j$ ∈ [0,1] |
| $d_{ij}$ | Euclidean distance between cell $i$ and ECM element $j$ |
| $s_0$ | Rest length (`mEcmSpringRestLength`, default 0) |
| $\hat{\mathbf{d}}_{ij}$ | Unit vector from cell to ECM element |

**Regime behaviour:**
- $d_{ij} > s_0$: adhesive pull toward ECM (confinement)
- $d_{ij} < s_0$: repulsive push away from ECM

The total force on cell $i$ is $\mathbf{F}_i = \sum_j \mathbf{F}_{ij}$.

### ECM dynamics
- **Degradation:** cells reduce $\rho$ at their position (MMP secretion), rate $\propto \Delta t$
- **Fibre remodeling:** cell traction (radially outward from centroid) rotates local fibre direction
- **Diffusion:** density field is smoothed each timestep ($D > 0$)

### 3D variant (`ECMConfinementForce3d`)
Uses a simpler radial-only formulation:

$$\mathbf{F}_i = -k \, \rho(\mathbf{x}_i) \, \hat{\mathbf{r}}_i$$

where $\hat{\mathbf{r}}_i$ points radially outward from the tissue centroid. No per-element summation — just a single density sample at the cell position.

### Key characteristics
- **Grid is fixed** — ECM elements do not move
- **Purely linear** spring law (no log/exp corrections)
- Forces are **reactive/passive** — the ECM resists cell displacement but does not guide migration direction
- No Newton's 3rd law reaction on the grid

---

## 2. DynamicECMContactGuidanceForce (Grid-Based Migration Guidance)

**Purpose:** Drive active cell migration along ECM fibre directions, modulated by local density and anisotropy.

**ECM representation:** Same fixed Eulerian grid as above (`DynamicECMField` / `DynamicECMField3d`).

### Force law

The force on cell $i$ at position $\mathbf{r}_i$ is:

$$\mathbf{F}_i = v_0 \, s \, \sqrt{\rho(\mathbf{r}_i)} \, \hat{\mathbf{n}}_{\text{mig}}$$

| Symbol | Meaning |
|--------|---------|
| $v_0$ | Base migration speed (`mBaseSpeed`, default 1.25) |
| $s$ | ECM sensitivity multiplier (`mECMSensitivity`, default 1.0) |
| $\rho$ | Local ECM density ∈ [0,1] |
| $\hat{\mathbf{n}}_{\text{mig}}$ | Composite migration direction (unit vector, see below) |

### Migration direction decomposition

The migration direction is a weighted superposition of three components:

$$\hat{\mathbf{n}}_{\text{mig}} = \text{normalise}\!\Big(\,
  \alpha_{\text{eff}} \, \mathbf{m}_{\text{fibre}}
  + (1-\alpha_{\text{eff}}) \cdot 0.3 \; \mathbf{m}_{\perp}
  + (1-\alpha_{\text{eff}}) \cdot 0.5 \; \mathbf{m}_{\text{rand}}
\Big)$$

where $\alpha_{\text{eff}} = \alpha \cdot \rho$ is the effective anisotropy (product of fibre alignment $\alpha$ and density).

| Component | Formula | Description |
|-----------|---------|-------------|
| $\mathbf{m}_{\text{fibre}}$ | $\pm \hat{\mathbf{f}}$ | Bidirectional motion along/against fibre (sign chosen randomly each step) |
| $\mathbf{m}_{\perp}$ | $\xi \, \hat{\mathbf{f}}_\perp$ | Diffusion perpendicular to fibre ($\xi \sim U(-1,1)$) |
| $\mathbf{m}_{\text{rand}}$ | Random unit vector | Isotropic random walk |

**Limiting cases:**
- Dense, aligned ECM ($\alpha_{\text{eff}} \to 1$): motion is almost purely along fibres
- Sparse/isotropic ECM ($\alpha_{\text{eff}} \to 0$): dominated by random walk

### ECM dynamics
- **Degradation:** cells reduce $\rho$ at their position
- **Fibre remodeling:** cell traction rotates local fibre toward the pulling direction
- **Deposition (optional):** cells deposit new ECM aligned with their migration direction
- **Diffusion:** density field is smoothed each timestep

### Key characteristics
- **Active force** — drives cell motility (not a passive resistance)
- Speed scales as $\sqrt{\rho}$ — cells slow down in degraded ECM
- **Stochastic** — random components introduce Brownian-like migration noise
- No pairwise cell–ECM spring interactions; force depends only on local field values

---

## 3. GhostNodeEcmForce (Particle-Based Viscoelastic ECM)

**Purpose:** Model the ECM as a network of movable particles ("ghost nodes") connected by viscoelastic springs, providing mechanically realistic cell–ECM and ECM–ECM interactions.

**ECM representation:** Off-lattice Lagrangian particles (`GhostNodeEcmField`). Each ghost node carries position, velocity, density $\rho \in [0,1]$, fibre direction $\hat{\mathbf{f}}$, anisotropy $\alpha$, and a neighbour connectivity list.

*Inspired by Carrasco-Mantis et al. (2024), doi:10.1016/j.compbiomed.2024.109559*

### Cell–ghost node force law

For each cell $i$ interacting with nearby ghost node $j$ (within cutoff $r_c$):

**Compression** ($d_{ij} \leq s_0$) — log-repulsion, density-independent:

$$F_{ij} = k_{cg} \, s_0 \ln\!\left(1 + \frac{d_{ij} - s_0}{s_0}\right) \cdot A_{ij}$$

**Extension** ($d_{ij} > s_0$) — exponential-decay attraction, density-weighted:

$$F_{ij} = k_{cg} \, \rho_j \, (d_{ij} - s_0) \, e^{-\alpha_{\text{decay}} (d_{ij} - s_0)/s_0} \cdot A_{ij}$$

**Zero rest length** ($s_0 = 0$):

$$F_{ij} = k_{cg} \, \rho_j \, d_{ij} \, e^{-\alpha_{\text{decay}} \, d_{ij}/r_c} \cdot A_{ij}$$

The directed force vector is $\mathbf{F}_{ij} = F_{ij} \, \hat{\mathbf{d}}_{ij}$ (cell → ghost direction).

| Symbol | Meaning |
|--------|---------|
| $k_{cg}$ | Cell–ghost stiffness (`mCellGhostStiffness`, default 10.0) |
| $\rho_j$ | Ghost node density ∈ [0,1] |
| $s_0$ | Rest length (`mCellGhostRestLength`) |
| $\alpha_{\text{decay}}$ | Exponential decay rate (hardcoded = 5.0) |
| $A_{ij}$ | Anisotropy factor (see below) |

### Anisotropy factor

$$A_{ij} = 1 - a \, \alpha_j \, \bigl(1 - |\cos\theta_{ij}|\bigr)$$

where $\theta_{ij}$ is the angle between the cell–ghost displacement vector and the ghost node's fibre direction, $\alpha_j$ is the node's anisotropy, and $a = 0.5$ is the anisotropy strength. This means:
- Along fibres ($\cos\theta = 1$): $A = 1$ — full stiffness
- Across fibres ($\cos\theta = 0$): $A = 1 - a\alpha_j$ — reduced stiffness (easier to separate across fibres)

### Ghost–ghost (ECM–ECM) interactions

Ghost nodes are connected to their topological neighbours by viscoelastic springs. The ghost–ghost dynamics are integrated as overdamped motion:

$$\eta \, \dot{\mathbf{x}}_j = \sum_{k \in \mathcal{N}(j)} k_{gg} \, (d_{jk} - s_0^{gg}) \, A_{jk} \, \hat{\mathbf{d}}_{jk} + \mathbf{F}_j^{\text{ext}}$$

where $\eta$ is the damping coefficient (`mGhostDamping`) and $\mathbf{F}_j^{\text{ext}}$ includes Newton's 3rd law reactions from cell–ghost springs.

### Newton's 3rd law

Unlike the other two forces, the GhostNodeEcmForce applies **equal-and-opposite reactions** to ghost nodes. When a cell pushes against ECM, the ghost node is pushed outward, causing the ECM mesh to deform.

### ECM dynamics
- **Degradation:** cells reduce $\rho_j$ of nearby ghost nodes (MMP secretion)
- **Fibre remodeling:** cell traction (radially outward from centroid) rotates ghost node fibre directions
- **Node removal:** ghost nodes with $\rho < \rho_{\text{threshold}}$ are periodically deleted
- **No diffusion** — density is carried by discrete particles, not a continuum field

### Key characteristics
- **Lagrangian** — ECM particles move, deform, and can be removed
- **Nonlinear** force law (log compression + exponential extension), matching Chaste's `GeneralisedLinearSpringForce` conventions
- **Bidirectional mechanical coupling** — cells deform the ECM mesh and ECM resists cell motion
- **Anisotropic** — fibre direction modulates force magnitude
- **No continuum diffusion** — information propagates through the spring network

---

## Side-by-Side Summary

| Property | ECMConfinementForce | DynamicECMContactGuidanceForce | GhostNodeEcmForce |
|----------|--------------------|---------------------------------|-------------------|
| **ECM representation** | Fixed grid (Eulerian) | Fixed grid (Eulerian) | Movable particles (Lagrangian) |
| **Force type** | Passive spring (confinement) | Active migration (motility) | Passive spring (viscoelastic) |
| **Force law** | Linear: $k \rho (d - s_0)$ | $v_0 s \sqrt{\rho} \, \hat{n}$ | Log/exp nonlinear + anisotropy |
| **Interaction** | Cell ↔ grid elements | Cell ← local field value | Cell ↔ ghost node (Newton III) |
| **ECM deformation** | No | No | Yes (ghost nodes move) |
| **Fibre anisotropy in force** | No | Yes (migration direction) | Yes (force magnitude) |
| **Stochastic component** | No | Yes (random walk + perp diffusion) | No |
| **Node removal** | No | No | Yes (depleted nodes deleted) |
| **Density dependence** | Linear ($\rho$) | Square root ($\sqrt{\rho}$) | Linear ($\rho$) in extension; none in compression |
| **Density diffusion** | Yes (grid PDE) | Yes (grid PDE) | No (discrete particles) |
| **Typical use** | Organoid boundary confinement | Cell migration in stroma | Mechanically realistic ECM shell |
| **Dimension** | 2D templated (+3D variant) | 2D templated + 3D class | 2D/3D templated |

---

## When to Use Which

- **ECMConfinementForce** — When you need a simple confining boundary that cells can locally degrade to enable budding. Cheapest computationally.
- **DynamicECMContactGuidanceForce** — When cell migration direction should be guided by ECM fibre alignment (e.g. invasion into stroma). Does not confine; it *drives* motion.
- **GhostNodeEcmForce** — When the ECM itself must deform mechanically under cell forces (e.g. ECM remodeling, matrix stretching, realistic confinement with fibre anisotropy). Most expensive but most physically faithful.
