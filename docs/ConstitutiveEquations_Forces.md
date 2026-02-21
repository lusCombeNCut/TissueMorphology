# Constitutive Equations — Forces

All force models used in the **CryptBudding** organoid simulation app, with their
mathematical formulations, default parameters, and literature references.

---

## Table of Contents

| # | Force | Models |
|---|-------|--------|
| 1 | [Generalised Linear Spring Force](#1-generalised-linear-spring-force) | node2d, node3d |
| 2 | [Differential Adhesion Force](#2-differential-adhesion-force) | vertex2d, vertex3d † |
| 3 | [Nagai–Honda Force](#3-nagaihonda-force) | vertex2d |
| 4 | [Surface Tension Force](#4-surface-tension-force-3d-monolayer-vertex) | vertex3d |
| 5 | [Basement Membrane Force](#5-basement-membrane-force) | all |
| 6 | [Lumen Pressure Force](#6-lumen-pressure-force) | node2d, vertex2d, node3d |
| 7 | [Lumen Pressure SubForce](#7-lumen-pressure-subforce-3d-monolayer-vertex) | vertex3d |
| 8 | [Apical Constriction Force](#8-apical-constriction-force) | vertex2d |
| 9 | [Dynamic ECM Contact Guidance Force](#9-dynamic-ecm-contact-guidance-force-3d) | node3d, vertex3d |
| 10 | [Overdamped Equation of Motion](#10-overdamped-equation-of-motion) | all |
| 11 | [Sloughing (Plane-Based Cell Killer)](#11-sloughing-plane-based-cell-killer) | node2d, vertex2d, node3d |

---

## 1. Generalised Linear Spring Force

**Source:** `Chaste/cell_based/src/population/forces/GeneralisedLinearSpringForce.{hpp,cpp}`

A pairwise cell–cell interaction force derived from the Meineke spring model.
The force law differs for mesh-based and node-based populations.

### Constitutive Equations

**Rest-length evolution after division**
([GeneralisedLinearSpringForce.cpp, lines 155–165](../../../cell_based/src/population/forces/GeneralisedLinearSpringForce.cpp)):

$$
s_{ij}(t) = \lambda + (s_{ij}^{\infty} - \lambda)\,\frac{t}{\tau_{\text{grow}}}
\qquad \text{for } t < \tau_{\text{grow}}
$$

where $\lambda$ is the division resting spring length (`mMeinekeDivisionRestingSpringLength`, default 0.5),
$s_{ij}^{\infty}$ is the equilibrium rest length, and $\tau_{\text{grow}}$ is the
spring growth duration (`mMeinekeSpringGrowthDuration`, default 1.0 h).

**Apoptosis rest-length reduction**
([GeneralisedLinearSpringForce.cpp, lines 176–204](../../../cell_based/src/population/forces/GeneralisedLinearSpringForce.cpp)):

The total rest length is split into two per-cell contributions
$s_A$ and $s_B$, representing each cell's share of the spring:

*Mesh-based populations:*

$$
s_A = \tfrac{1}{2}\,s, \qquad s_B = \tfrac{1}{2}\,s
$$

*Node-based populations* (proportional to cell radii $r_A$, $r_B$):

$$
s_A = \frac{r_A}{r_A + r_B}\,s, \qquad s_B = \frac{r_B}{r_A + r_B}\,s
$$

If cell $A$ has begun apoptosis, its contribution shrinks linearly to zero over
the apoptosis time $T_{\text{apop}}$:

$$
s_A \leftarrow s_A \times \frac{t_{\text{until\ death}}^{A}}{T_{\text{apop}}}
$$

(and similarly for $s_B$ if cell $B$ is apoptotic). The final rest length is:

$$
s_{ij} = s_A + s_B
$$

**Mesh-based populations (linear spring)**
([GeneralisedLinearSpringForce.cpp, line 213](../../../cell_based/src/population/forces/GeneralisedLinearSpringForce.cpp)):

$$
\mathbf{F}_{ij} = \mu\, m\, (d_{ij} - s_{ij})\, \hat{\mathbf{r}}_{ij}
$$

**Node-based populations**
([GeneralisedLinearSpringForce.cpp, lines 217–229](../../../cell_based/src/population/forces/GeneralisedLinearSpringForce.cpp)):

*Repulsion* ($d_{ij} < s_{ij}$):

$$
\mathbf{F}_{ij} = \mu\, m\, s_0 \ln\!\left(1 + \frac{d_{ij} - s_{ij}}{s_0}\right) \hat{\mathbf{r}}_{ij}
$$

*Attraction* ($d_{ij} \ge s_{ij}$):

$$
\mathbf{F}_{ij} = \mu\, m\, (d_{ij} - s_{ij})\, \exp\!\left(-\frac{5\,(d_{ij} - s_{ij})}{s_0}\right) \hat{\mathbf{r}}_{ij}
$$

| Symbol | Meaning | Default |
|--------|---------|---------|
| $\mu$ | Spring stiffness (`mMeinekeSpringStiffness`) | 15.0 (Chaste default); set to 30.0 in CryptBudding |
| $m$ | `VariableSpringConstantMultiplicationFactor` | 1.0 |
| $d_{ij}$ | Distance between nodes $i$ and $j$ | — |
| $s_{ij}$ | Rest length | 1.0 (or sum of radii) |
| $s_0$ | Final rest length (before apoptosis correction) | — |
| $\alpha$ | Exponential decay constant | 5.0 |

### References

- Meineke, F.A., Potten, C.S. & Loeffler, M. (2001). Cell migration and organization in the intestinal crypt using a lattice-free model. *Cell Proliferation*, 34(4), 253–266. doi:[10.1046/j.0960-7722.2001.00216.x](https://doi.org/10.1046/j.0960-7722.2001.00216.x)
- Chaste documentation: [GeneralisedLinearSpringForce](https://chaste.github.io/chaste-docs/)

---

## 2. Differential Adhesion Force

**Source:** `Chaste/projects/TissueMorphology/src/DifferentialAdhesionForce.hpp` (header-only, extends `GeneralisedLinearSpringForce`)

Modifies the spring stiffness by a cell-type-dependent adhesion multiplier.
Cells are classified as *apical* or *basal* via the `is_apical` cell data item.

### Constitutive Equations

The force between nodes $i$ and $j$ is identical to the Generalised Linear Spring
Force (§1), but the multiplication factor $m$ is replaced by
([DifferentialAdhesionForce.hpp, `VariableSpringConstantMultiplicationFactor`](../../../projects/TissueMorphology/src/DifferentialAdhesionForce.hpp)):

$$
m_{\text{type}}(i,j) = \begin{cases}
\alpha_{aa} & \text{if both apical} \\
\alpha_{bb} & \text{if both basal} \\
\alpha_{ab} & \text{if mixed (apical–basal)}
\end{cases}
$$

| Symbol | Parameter | Class Default | CryptBudding |
|--------|-----------|---------------|---------------|
| $\alpha_{aa}$ | `apicalApicalAdhesion` | 1.0 | 1.2 |
| $\alpha_{bb}$ | `basalBasalAdhesion` | 1.0 | 1.0 |
| $\alpha_{ab}$ | `apicalBasalAdhesion` | 0.5 | 0.5 |

Hence apical–apical interactions are 20% stronger than baseline,
while cross-type interactions are 50% weaker — driving cell-type sorting.

### References

- Steinberg, M.S. (1963). Reconstruction of tissues by dissociated cells. *Science*, 141(3579), 401–408. (Differential adhesion hypothesis)
- Graner, F. & Glazier, J.A. (1992). Simulation of biological cell sorting using a two-dimensional extended Potts model. *Physical Review Letters*, 69(13), 2013.

---

## 3. Nagai–Honda Force

**Source:** `Chaste/cell_based/src/population/forces/NagaiHondaForce.{hpp,cpp}`

A vertex-model force law for 2D vertex-based cell populations, derived from the
Nagai–Honda energy functional.

### Constitutive Equations

**Energy functional:**

$$
U = \sum_{\alpha} \left[
  \lambda_d \left(A_\alpha - A_\alpha^0\right)^2
+ \lambda_m \left(L_\alpha - L_\alpha^0\right)^2
+ \sum_{\langle ij \rangle \in \alpha} \gamma_{ij}\, \ell_{ij}
\right]
$$

**Force on vertex $i$:**

$$
\mathbf{F}_i = -\nabla_i U
$$

This decomposes into three contributions:

1. **Deformation (area elasticity):**

$$
\mathbf{F}_i^{\text{area}} = -2\,\lambda_d \sum_{\alpha \ni i}
\left(A_\alpha - A_\alpha^0\right) \nabla_i A_\alpha
$$

2. **Membrane surface tension (perimeter elasticity):**

$$
\mathbf{F}_i^{\text{perim}} = -2\,\lambda_m \sum_{\alpha \ni i}
\left(L_\alpha - L_\alpha^0\right) \nabla_i L_\alpha
$$

where $L_\alpha^0 = 2\sqrt{\pi\, A_\alpha^0}$ (isoperimetric target perimeter).

3. **Adhesion (line tension):**

$$
\mathbf{F}_i^{\text{adh}} = -\sum_{\alpha \ni i}\sum_{j \sim i}
\gamma_{ij}\, \nabla_i \ell_{ij}
$$

where

$$
\gamma_{ij} = \begin{cases}
\gamma_{\text{cell-cell}} & \text{if edge shared by two cells} \\
\gamma_{\text{boundary}} & \text{if edge on free boundary}
\end{cases}
$$

| Symbol | Parameter | Default | CryptBudding |
|--------|-----------|---------|--------------|
| $\lambda_d$ | `NagaiHondaDeformationEnergyParameter` | 100.0 | `ecmStiffness` |
| $\lambda_m$ | `NagaiHondaMembraneSurfaceEnergyParameter` | 10.0 | 10.0 |
| $\gamma_{\text{cell-cell}}$ | `NagaiHondaCellCellAdhesionEnergyParameter` | 0.5 | 1.0 |
| $\gamma_{\text{boundary}}$ | `NagaiHondaCellBoundaryAdhesionEnergyParameter` | 1.0 | 2.0 |
| $A_\alpha^0$ | Target area (from `SimpleTargetAreaModifier`) | 1.0 | computed from mesh geometry |

### References

- Nagai, T. & Honda, H. (2001). A dynamic cell model for the formation of epithelial tissues. *Philosophical Magazine B*, 81(7), 699–719. doi:[10.1080/13642810108205772](https://doi.org/10.1080/13642810108205772)
- Farhadifar, R., Röper, J.C., Aigouy, B., Eaton, S. & Jülicher, F. (2007). The influence of cell mechanics, cell-cell interactions, and proliferation on epithelial packing. *Current Biology*, 17(24), 2095–2104. doi:[10.1016/j.cub.2007.11.049](https://doi.org/10.1016/j.cub.2007.11.049)

---

## 4. Surface Tension Force (3D Monolayer Vertex)

**Source:** `Chaste/projects/OrganoidChaste/src/monolayer_population/forces/SurfaceTensionForce.{hpp,cpp}`

A 3D monolayer vertex model force for the OrganoidChaste finite-thickness vertex
framework. Each cell is a 3D polyhedron with distinguishable apical, basal, and
lateral faces, each carrying an independent surface tension.

### Constitutive Equations

**Energy functional**
([SurfaceTensionForce.cpp, `AddForceContribution`, lines 70–540](../../../projects/OrganoidChaste/src/monolayer_population/forces/SurfaceTensionForce.cpp)):

where $\gamma_f \in \{\gamma_{\text{apical}},\, \gamma_{\text{basal}},\, \gamma_{\text{lateral}}\}$
depending on face type.

**Per-cell-type surface tensions** are assigned via mutation-state scaling:

$$
\gamma_f^{(\text{type})} = \gamma_f \times s_{\text{type}}
$$

| Cell type | Mutation state | Scale factor $s$ | Default |
|-----------|---------------|-------------------|---------|
| Stem | WildType | `gammaStemScale` | 0.7 |
| Transit-amplifying | ApcOneHit | `gammaTransitScale` | 1.0 |
| Differentiated / Paneth | ApcTwoHit | `gammaDiffScale` | 1.3 |

**Force on vertex $i$ (unconstrained):**

$$
\mathbf{F}_i = -\sum_{f \ni i} \gamma_f \, \nabla_i A_f
$$

**Volume-conserving projection:**

The force is projected onto the manifold of constant cell volumes using Lagrange
multipliers. Define $M_{\alpha\beta} = \sum_i (\nabla_i V_\alpha) \cdot (\nabla_i V_\beta)$,
then solve:

$$
\sum_\beta M_{\alpha\beta}\, \lambda_\beta^{\text{proj}} = \sum_i \mathbf{F}_i \cdot \nabla_i V_\alpha
$$

The projected force is:

$$
\mathbf{F}_i^{\text{proj}} = \mathbf{F}_i - \sum_\alpha \lambda_\alpha^{\text{proj}}\, \nabla_i V_\alpha
$$

**Volume correction:**

An additional displacement enforces target volumes:

$$
\sum_\beta M_{\alpha\beta}\, \mu_\beta = -(V_\alpha - V_\alpha^0)
$$

$$
\Delta\mathbf{x}_i = \sum_\alpha \mu_\alpha\, \nabla_i V_\alpha
$$

**Simulated annealing (optional):**

Metropolis random perturbations with Boltzmann acceptance
$P(\Delta E) = \exp(-\Delta E / T)$, where $T = T_0\,(1 - t/\tau)$.

| Symbol | Parameter | Default |
|--------|-----------|---------|
| $\gamma_{\text{apical}}$ | `gammaApical` | 0.85 |
| $\gamma_{\text{basal}}$ | `gammaBasal` | 0.85 |
| $\gamma_{\text{lateral}}$ | `gammaLateral` | 0.70 |

### References

- Hannezo, E., Prost, J. & Joanny, J.F. (2014). Theory of epithelial sheet morphology in three dimensions. *Proceedings of the National Academy of Sciences*, 111(1), 27–32. doi:[10.1073/pnas.1312076111](https://doi.org/10.1073/pnas.1312076111)
- Rauzi, M. et al. (2008). Nature and anisotropy of cortical forces orienting *Drosophila* tissue morphogenesis. *Nature Cell Biology*, 10(12), 1401–1410.
- OrganoidChaste framework (Oliver et al.)

---

## 5. Basement Membrane Force

**Source:** `Chaste/projects/TissueMorphology/src/BasementMembraneForce.hpp` (header-only, templated on DIM)

A radial confinement force that models the extracellular matrix (ECM) / basement
membrane as an elastic shell surrounding the organoid. Only active when a cell
exceeds the membrane radius.

### Constitutive Equations

([BasementMembraneForce.hpp, `AddForceContribution`](../../../projects/TissueMorphology/src/BasementMembraneForce.hpp))

$$
\mathbf{F}_i = \begin{cases}
-k_i \left(r_i - R_{\text{BM}}\right) \hat{\mathbf{r}}_i
  & \text{if } r_i > R_{\text{BM}} \\[4pt]
\mathbf{0} & \text{otherwise}
\end{cases}
$$

where:

- $r_i = \|\mathbf{x}_i - \mathbf{x}_c\|$ — distance from organoid centre
- $\hat{\mathbf{r}}_i = (\mathbf{x}_i - \mathbf{x}_c) / r_i$ — outward unit radial vector
- $k_i$ — per-cell stiffness (from cell data `basement_membrane_stiffness`, or default $k$)
- $R_{\text{BM}}$ — current basement membrane radius

**ECM degradation (time-dependent softening)**
([BasementMembraneForce.hpp, `UpdateBasementMembraneRadius`](../../../projects/TissueMorphology/src/BasementMembraneForce.hpp)):

$$
R_{\text{BM}}(t) = \min\!\big(R_0 + \dot{R}\, t,\; R_{\max}\big)
$$

| Symbol | Parameter | Default |
|--------|-----------|---------|
| $k$ | `bmStiffnessNode` or `bmStiffnessVertex` | $= $ `ecmStiffness` (node) or `ecmStiffness × 0.5` (vertex) |
| $R_0$ | Initial BM radius (2D: $R_{\text{org}} + 2$; 3D: $R_{\text{org}} + 2$) | model-specific |
| $\dot{R}$ | `ecmDegradationRate` | 0.02 |
| $R_{\max}$ | Maximum degraded radius | $4 \times R_{\text{org}}$ |

For vertex populations (2D and 3D), the force is distributed equally across all
nodes of the element:

$$
\mathbf{f}_{\text{node}} = \frac{\mathbf{F}_i}{N_{\text{nodes}}}
$$

**Organoid centre tracking:**

$\mathbf{x}_c$ is dynamically updated as the population centroid at each time step.

### References

- Dunn, S.J., Näthke, I.S. & Osborne, J.M. (2013). Computational models reveal a passive mechanism for cell migration in the crypt. *PLoS ONE*, 8(11), e80516. doi:[10.1371/journal.pone.0080516](https://doi.org/10.1371/journal.pone.0080516)
- Almet, A.A., Maini, P.K., Moulton, D.E. & Byrne, H.M. (2021). Modelling perspectives on the intestinal crypt, a canonical wnt signalling-driven system. *Journal of the Royal Society Interface*, 18, 20210013.

---

## 6. Lumen Pressure Force

**Source:** `Chaste/projects/TissueMorphology/src/LumenPressureForce.hpp` (header-only, templated on DIM)

Models hydrostatic pressure from the lumen interior. Pushes cells **outward** when
they are closer to the lumen centre than the equilibrium radius. Used in node-based
and 2D vertex models.

### Constitutive Equations

([LumenPressureForce.hpp, `AddForceContribution`](../../../projects/TissueMorphology/src/LumenPressureForce.hpp))

$$
\mathbf{F}_i = \begin{cases}
P \left(R_{\text{eq}} - r_i\right) \hat{\mathbf{r}}_i
  & \text{if } r_i < R_{\text{eq}} \\[4pt]
\mathbf{0} & \text{if } r_i \ge R_{\text{eq}}
\end{cases}
$$

where:

- $P$ — pressure strength (`mPressureStrength`, default 2.0)
- $R_{\text{eq}}$ — lumen equilibrium radius (`mLumenEquilibriumRadius`, default $R_{\text{org}} + 1$)
- $r_i = \|\mathbf{x}_i - \mathbf{x}_L\|$ — distance from lumen centre
- $\hat{\mathbf{r}}_i$ — outward unit radial vector from lumen centre

For vertex populations the force is again distributed equally across all element
nodes.

**Centre tracking:**

$\mathbf{x}_L$ is dynamically updated as the population centroid when `mTrackCenter = true`.

| Symbol | Parameter | Default |
|--------|-----------|---------|
| $P$ | `lumenPressure` | 2.0 |
| $R_{\text{eq}}$ | `lumenEqRadius2d` / `lumenEqRadius3d` | $R_{\text{org}} + 1$ |

### References

- Yan, K.S. et al. (2017). Intestinal enteroendocrine lineage cells possess homeostatic and injury-inducible stem cell activity. *Cell Stem Cell*, 21, 78–90.
- Dahl-Jensen, S. & Bhatt, D.K. (2017). Lumen formation in the developing pancreas: mechanics, signaling, and disease. *Mechanisms of Development*, 144, 55–63.

---

## 7. Lumen Pressure SubForce (3D Monolayer Vertex)

**Source:** `Chaste/projects/OrganoidChaste/src/monolayer_population/forces/sub_forces/LumenPressureSubForce.{hpp,cpp}`

A sub-force within the OrganoidChaste `SurfaceTensionForce` framework. Instead of
the simplified radial formulation (§6), it uses the **volumetric lumen gradient**
computed directly from the 3D mesh geometry.

### Constitutive Equations

([LumenPressureSubForce.cpp, `GetForceContributions`](../../../projects/OrganoidChaste/src/monolayer_population/forces/sub_forces/LumenPressureSubForce.cpp))

$$
\mathbf{F}_i = P \; \nabla_i V_{\text{lumen}}
$$

where:

- $P$ — pressure constant (`mPressureConstant`, default = `lumenPressure` = 2.0)
- $\nabla_i V_{\text{lumen}}$ — gradient of the enclosed lumen volume with respect to vertex $i$'s position, computed by `MutableMonolayerVertexMesh::CalculateLumenVolGradient`

This is the thermodynamic conjugate force to the lumen volume — i.e. the same force
that would arise from a uniform internal pressure $P$ acting on the apical surface.

| Symbol | Parameter | Default |
|--------|-----------|---------|
| $P$ | `lumenPressure` | 2.0 |

### References

- Hannezo, E., Prost, J. & Joanny, J.F. (2014). Theory of epithelial sheet morphology in three dimensions. *PNAS*, 111(1), 27–32.
- OrganoidChaste framework

---

## 8. Apical Constriction Force

**Source:** `Chaste/projects/TissueMorphology/src/ApicalConstrictionForce.hpp` (header-only, templated on DIM)

Models actomyosin-driven apical constriction that narrows the apical surface of
epithelial cells. Only acts on cells labelled as *apical* (`is_apical > 0.5`).
The force drives cells **inward** toward the population centroid.

### Constitutive Equations

For each apical cell where $A_i > A_{\text{target}}$
([ApicalConstrictionForce.hpp, `AddForceContribution`](../../../projects/TissueMorphology/src/ApicalConstrictionForce.hpp)):

$$
\mathbf{F}_i = -\kappa \left(A_i - A_{\text{target}}\right) \hat{\mathbf{r}}_i
$$

where:

- $\kappa$ — constriction strength (`mConstrictionStrength`, default 3.0 in CryptBudding)
- $A_i$ — current cell area (estimated from Voronoi / element geometry)
- $A_{\text{target}} = A_{\text{ref}} \times (1 - f_{\text{red}})$ — target area after constriction
- $A_{\text{ref}}$ — reference area (captured at initialisation)
- $f_{\text{red}}$ — target area reduction fraction (`mTargetReduction`, default 0.5)
- $\hat{\mathbf{r}}_i = (\mathbf{x}_i - \bar{\mathbf{x}}) / \|\mathbf{x}_i - \bar{\mathbf{x}}\|$ — unit radial vector from centroid

The force is directed **inward** (toward the centroid), simulating the contractile
apical actomyosin ring observed in epithelial morphogenesis.

| Symbol | Parameter | Class Default | CryptBudding |
|--------|-----------|---------------|---------------|
| $\kappa$ | `apicalConstrictionStrength` | 5.0 | 3.0 |
| $f_{\text{red}}$ | `mTargetReduction` | 0.5 | 0.5 |

### References

- Martin, A.C. & Goldstein, B. (2014). Apical constriction: themes and variations on a cellular mechanism driving morphogenesis. *Development*, 141(10), 1987–1998. doi:[10.1242/dev.108043](https://doi.org/10.1242/dev.108043)
- Rauzi, M. & Lenne, P.F. (2011). Cortical forces in cell shape changes and tissue morphogenesis. *Current Topics in Developmental Biology*, 95, 93–144.

---

## 9. Dynamic ECM Contact Guidance Force (3D)

**Source:** `Chaste/projects/TissueMorphology/src/DynamicECMContactGuidanceForce3d.hpp` (header-only)

Models cell migration guided by an extracellular matrix (ECM) fibre field. The ECM
is represented by a `DynamicECMField3d` defined on a regular 3D grid storing local
fibre direction, density, and anisotropy. Only used in 3D models (`node3d`, `vertex3d`).

### Constitutive Equations

**Migration direction (stochastic)**
([DynamicECMContactGuidanceForce3d.hpp, `AddForceContribution`](../../../projects/TissueMorphology/src/DynamicECMContactGuidanceForce3d.hpp)):

$$
\hat{\mathbf{n}}_{\text{mig}} = \text{normalise}\!\Big[
  \alpha_{\text{eff}}\,(\pm\hat{\mathbf{f}})
+ (1 - \alpha_{\text{eff}})\cdot 0.3\;\xi\,\hat{\mathbf{f}}_\perp
+ (1 - \alpha_{\text{eff}})\cdot 0.5\;\hat{\mathbf{r}}_{\text{rand}}
\Big]
$$

where:

- $\hat{\mathbf{f}}$ — local fibre direction from ECM field (sign chosen randomly — bidirectional)
- $\hat{\mathbf{f}}_\perp$ — random vector perpendicular to $\hat{\mathbf{f}}$
- $\hat{\mathbf{r}}_{\text{rand}}$ — uniform random direction on the unit sphere (Marsaglia method)
- $\xi \in [-1, 1]$ — uniform random scalar
- $\alpha_{\text{eff}} = \alpha \cdot \rho$ — effective anisotropy (fibre anisotropy $\times$ ECM density)

**Force magnitude**
([DynamicECMContactGuidanceForce3d.hpp, lines ~70–90](../../../projects/TissueMorphology/src/DynamicECMContactGuidanceForce3d.hpp)):

$$
\mathbf{F}_i = v_0 \, s \, \sqrt{\rho(\mathbf{x}_i)} \; \hat{\mathbf{n}}_{\text{mig}}
$$

| Symbol | Parameter | Default |
|--------|-----------|---------|
| $v_0$ | `mBaseSpeed` | 0.3 |
| $s$ | `mECMSensitivity` | 1.0 |
| $\rho$ | ECM density at cell position | field-dependent |
| $\alpha$ | ECM anisotropy at cell position | field-dependent |

**ECM remodelling (coupled feedback):**

| Process | Rate | Effect |
|---------|------|--------|
| Degradation | 0.002 | Cells degrade ECM at their location |
| Remodelling | 0.05 | Cell traction forces reorient fibre direction |
| Deposition | 0.0003 | Cells deposit new ECM aligned with migration direction |
| Diffusion | per timestep | Spatial smoothing of ECM field |

### References

- Fraley, S.I. et al. (2015). Three-dimensional matrix fiber alignment modulates cell migration and MT1-MMP utility by spatially and temporally directing protrusions. *Scientific Reports*, 5, 14580. doi:[10.1038/srep14580](https://doi.org/10.1038/srep14580)
- Sapudom, J. et al. (2015). The phenotype of cancer cell invasion controlled by fibril diameter and pore size of 3D collagen networks. *Biomaterials*, 52, 367–375.
- Metzcar, J., Wang, Y., Heiland, R. & Macklin, P. (2019). A review of cell-based computational modeling in cancer biology. *JCO Clinical Cancer Informatics*, 3, 1–13. doi:[10.1200/CCI.18.00069](https://doi.org/10.1200/CCI.18.00069)

---

## 10. Overdamped Equation of Motion

**Source:** `Chaste/cell_based/src/simulation/numerical_methods/ForwardEulerNumericalMethod.{hpp,cpp}`
and `Chaste/cell_based/src/simulation/numerical_methods/AbstractNumericalMethod.{hpp,cpp}`

All models in CryptBudding use **overdamped (first-order) dynamics** — inertia is
neglected and cells move in the direction of the net force, at a speed proportional
to that force. This is the fundamental equation that translates force contributions
into cell displacements.

### Constitutive Equations

**Continuous form**
([AbstractNumericalMethod.cpp, `ComputeForcesIncludingDamping`](../../../cell_based/src/simulation/numerical_methods/AbstractNumericalMethod.cpp)):

$$
\eta_i \, \dot{\mathbf{x}}_i = \mathbf{F}_i
$$

**Discrete form (Forward Euler)**
([ForwardEulerNumericalMethod.cpp, `UpdateAllNodePositions`](../../../cell_based/src/simulation/numerical_methods/ForwardEulerNumericalMethod.cpp)):

$$
\mathbf{x}_i^{n+1} = \mathbf{x}_i^{n} + \frac{\Delta t}{\eta_i} \, \mathbf{F}_i
$$

where:

- $\eta_i$ — viscous drag coefficient for node $i$ (`GetDampingConstant`)
- $\mathbf{F}_i = \sum_k \mathbf{F}_i^{(k)}$ — total force on node $i$, summed over all registered force objects
- $\Delta t$ — simulation time step (`dt`)

The ratio $\mathbf{F}_i / \eta_i$ has units of velocity; multiplying by $\Delta t$
gives the displacement applied at each step.

| Symbol | Parameter | Default |
|--------|-----------|---------|
| $\eta$ | Damping constant (normal cells) | 1.0 |
| $\eta_{\text{mutant}}$ | Damping constant (mutant cells) | 1.0 |
| $\Delta t$ | Time step | model-dependent (0.005 – 0.01 for node; 0.0001 – 0.0005 for vertex) |

### References

- Pathmanathan, P. et al. (2009). A computational study of discrete mechanical tissue models. *Physical Biology*, 6(3), 036001. doi:[10.1088/1478-3975/6/3/036001](https://doi.org/10.1088/1478-3975/6/3/036001)
- Fletcher, A.G., Osterfield, M., Baker, R.E. & Shvartsman, S.Y. (2014). Vertex models of epithelial morphogenesis. *Biophysical Journal*, 106(11), 2291–2304. doi:[10.1016/j.bpj.2013.11.4498](https://doi.org/10.1016/j.bpj.2013.11.4498)

---

## 11. Sloughing (Plane-Based Cell Killer)

**Source:** `Chaste/cell_based/src/population/killers/PlaneBasedCellKiller.{hpp,cpp}`

Used in node2d, vertex2d, and node3d to remove cells that escape the simulation
domain. CryptBudding places bounding-box killers — pairs of planes on each axis —
forming a cube (3D) or square (2D) centred at the origin.

### Killing Criterion

([PlaneBasedCellKiller.cpp, `CheckAndLabelCellsForApoptosisOrDeath`](../../../cell_based/src/population/killers/PlaneBasedCellKiller.cpp))

A cell is killed (immediately removed) if its centre lies on the outward side of
any bounding plane:

$$
(\mathbf{x}_{\text{cell}} - \mathbf{x}_{\text{plane}}) \cdot \hat{\mathbf{n}} > 0
$$

In CryptBudding, for each spatial dimension $d$:
- Positive plane at $+H$ with outward normal $+\hat{\mathbf{e}}_d$
- Negative plane at $-H$ with outward normal $-\hat{\mathbf{e}}_d$

where $H = R_{\text{org}} \times f_{\text{slough}}$ and $f_{\text{slough}}$ =
`sloughRadiusFactor` (default 5.0).

| Symbol | Parameter | Default |
|--------|-----------|---------|
| $H$ | Half-width of bounding box | $5 \times R_{\text{org}}$ |
| $f_{\text{slough}}$ | `sloughRadiusFactor` | 5.0 |

> **Note:** vertex3d does not use sloughing because `MutableMonolayerVertexMesh`
> does not support `DeleteElementPriorToReMesh` in 3D.

### References

- Chaste documentation: [PlaneBasedCellKiller](https://chaste.github.io/chaste-docs/)

---

## Summary of Forces per Model Configuration

| Force | node2d | vertex2d | node3d | vertex3d |
|-------|:------:|:--------:|:------:|:--------:|
| Generalised Linear Spring | ✓ | — | ✓ | — |
| Differential Adhesion | — | ✓ (opt) | — | ✓ † |
| Nagai–Honda | — | ✓ | — | — |
| Surface Tension (3D) | — | — | — | ✓ |
| Basement Membrane | ✓ | ✓ | ✓ | ✓ |
| Lumen Pressure | ✓ (opt) | ✓ (opt) | ✓ (opt) | — |
| Lumen Pressure SubForce | — | — | — | ✓ (opt) |
| Apical Constriction | — | ✓ (opt) | — | — |
| Dynamic ECM Guidance | — | — | ✓ (opt) | ✓ (opt) |
| Overdamped EoM ($\eta\dot{\mathbf{x}}=\mathbf{F}$) | ✓ | ✓ | ✓ | ✓ |
| Sloughing (bounding box) | ✓ (opt) | ✓ (opt) | ✓ (opt) | — |

> † In **vertex3d**, differential adhesion is not implemented via the
> `DifferentialAdhesionForce` class (which is spring-based). Instead it is
> achieved through per-cell-type **surface tension scaling** in the
> `SurfaceTensionForce`: each cell type (stem / transit / differentiated) has its
> own scale factor (`gammaStemScale`, `gammaTransitScale`, `gammaDiffScale`)
> applied to the apical, basal, and lateral surface tensions.
>
> In **vertex2d**, differential adhesion could be implemented through the
> Nagai–Honda adhesion energy parameters, using the `is_apical` cell data item
> to distinguish cell types on the inner (lumen-facing) vs outer (ECM-facing)
> surface of the annular mesh.
>
> **Differential adhesion is not used in node-based models** because cells are
> point particles with no geometric faces — there is no meaningful apical/basal
> distinction.
>
> **Apical Constriction** is only meaningful for vertex-based models where cells
> have geometric extent. In vertex2d the annular mesh provides a clear inner
> (lumen-facing) and outer (ECM-facing) surface. In vertex3d, apical constriction
> is handled implicitly through differential apical/basal surface tensions in the
> `SurfaceTensionForce`. Node-based models lack the geometric structure required
> for a meaningful apical surface.

### Division Rules

| Model | Division rule | Monolayer preservation |
|-------|--------------|----------------------|
| node2d | N/A (centre-based: daughter placed along random direction) | Monolayer maintained by radial confinement forces |
| node3d | N/A (centre-based: daughter placed along random direction) | Monolayer maintained by radial confinement forces |
| vertex2d | `RadialVertexBasedDivisionRule` — divides along radial direction | ✓ daughters placed side-by-side tangentially |
| vertex3d | `MeanLongAxisMonolayerVertexBasedDivisionRule` — divides along mean apical/basal long axis | ✓ daughters placed side-by-side within monolayer plane |

### Neighbour Detection

| Model | Method | Cutoff |
|-------|--------|--------|
| node2d, node3d | **Overlapping spheres** — cells are neighbours when $d_{ij} \le r_i + r_j$, with candidate pairs found via a spatial box collection bounded by `interactionCutoff` | `interactionCutoff2d` / `interactionCutoff3d` |
| vertex2d | Topological — cells sharing an edge are neighbours | — |
| vertex3d | Topological — cells sharing a lateral face are neighbours | — |
