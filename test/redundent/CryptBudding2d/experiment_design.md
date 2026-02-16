# Experiment Design: 2D Crypt Budding from Organoids

## 1. Overview

This experiment investigates **intestinal crypt budding** from a 2D organoid
cross-section. The primary independent variable is **ECM stiffness**, which
controls the rigidity of the basement membrane surrounding the organoid. The
hypothesis is that lower ECM stiffness permits more outward budding events
(crypt formation) while higher stiffness suppresses them.

Two complementary cell-based models are used:

| Model | Framework | Cell representation | Key mechanics |
|-------|-----------|-------------------|---------------|
| **Node-based** | Overlapping spheres | Point particles + pairwise springs | DifferentialAdhesionForce (linear springs) |
| **Vertex-based** | Nagai–Honda vertex model | Polygonal cells (wedge quads) | NagaiHondaForce (area + perimeter energy) |

Both share the same custom forces for basement membrane confinement, lumen
pressure, and apical constriction.

---

## 2. Geometry: Annular Ring (Circular Organoid Cross-Section)

### Rationale

An earlier design placed cells on a **flat rectangular grid**, which produced
trivial results (exactly 1 "crypt" regardless of stiffness) because the flat
boundary conditions did not permit the mechanical buckling instability required
for budding. Real organoids are **spherical/tubular** with a central lumen
(Serra et al., 2019; Sato et al., 2009).

The redesign uses an **annular ring** of cells:
- **Apical surface** faces inward (toward the lumen)
- **Basal surface** faces outward (toward the ECM/basement membrane)
- Cells are arranged at radius $R_0$ (node-based) or between $R_\text{inner}$
  and $R_\text{outer}$ (vertex-based)

This captures the essential topology of an organoid cross-section in 2D and
allows azimuthal buckling instabilities (Hannezo et al., 2011).

### Node-based geometry

Cells are placed as point particles on a circle of radius $R_0 = 8.0$ with
small radial noise $\delta r \sim \mathcal{U}(-0.15, 0.15)$ to break perfect
symmetry.

### Vertex-based geometry

Cells are quadrilateral vertex elements forming a closed annular ring between
$R_\text{inner} = 6.0$ (apical) and $R_\text{outer} = 8.0$ (basal). Each
element subtends an angle $\Delta\theta = 2\pi / N$ where $N$ is the number
of cells. The target area of each element is:

$$A_\text{target} = \frac{\Delta\theta}{2} \left( R_\text{outer}^2 - R_\text{inner}^2 \right)$$

---

## 3. Constitutive Equations

### 3.1 Cell–Cell Interactions

#### Node-based: Generalised Linear Spring Force with Differential Adhesion

The pairwise interaction force between nodes $i$ and $j$ follows the
generalised linear spring model (Meineke et al., 2001; Dunn et al., 2012):

$$\mathbf{F}_{ij} = \mu_{ij} \cdot k_s \cdot \text{overlap}(r_{ij}, s_{ij}) \cdot \hat{\mathbf{r}}_{ij}$$

where $k_s$ is the base spring stiffness, $s_{ij}$ is the natural rest length,
$r_{ij} = \|\mathbf{x}_i - \mathbf{x}_j\|$ is the inter-node distance, and
$\mu_{ij}$ is a **differential adhesion multiplier** that depends on cell types:

| Pair type | Multiplier $\mu$ | Value |
|-----------|-------------------|-------|
| Apical–Apical | $\mu_{AA}$ | 1.2 |
| Basal–Basal | $\mu_{BB}$ | 1.0 |
| Apical–Basal | $\mu_{AB}$ | 0.5 |

This implements the **differential adhesion hypothesis** (Steinberg, 1963;
Foty & Steinberg, 2005), where stronger homotypic adhesion drives cell sorting
and layer maintenance.

**References:**
- Meineke, F.A. et al. (2001) *Cell Prolif.* 34(4):253–266
- Dunn, S.-J. et al. (2012) *J. Theor. Biol.* 298:82–91
- Steinberg, M.S. (1963) *Science* 141:401–408
- Foty, R.A. & Steinberg, M.S. (2005) *Dev. Biol.* 278(1):255–263

#### Vertex-based: Nagai–Honda Force

The total energy of the vertex model is (Nagai & Honda, 2001):

$$E = \sum_\alpha \left[ \frac{\lambda_A}{2}(A_\alpha - A_0)^2 + \frac{\lambda_P}{2} L_\alpha^2 + \sum_{\beta \sim \alpha} \gamma_{\alpha\beta} l_{\alpha\beta} \right]$$

where for each cell $\alpha$:
- $A_\alpha$ is the cell area, $A_0$ is the target area
- $L_\alpha$ is the cell perimeter
- $l_{\alpha\beta}$ is the shared edge length with neighbour $\beta$
- $\lambda_A$ is the **area elasticity** (deformation energy parameter) — **set equal to ECM_STIFFNESS**
- $\lambda_P$ is the **membrane surface energy** parameter
- $\gamma_{\alpha\beta}$ is the **cell–cell adhesion** energy (negative = adhesive)

Forces on each vertex are $\mathbf{F}_i = -\nabla_i E$.

**References:**
- Nagai, T. & Honda, H. (2001) *Phil. Mag. B* 81(7):699–719
- Farhadifar, R. et al. (2007) *Curr. Biol.* 17(24):2095–2104
- Fletcher, A.G. et al. (2014) *Biophys. J.* 106(11):2291–2304

### 3.2 Basement Membrane Force (ECM Confinement)

A radial restoring force confining cells within the basement membrane
(Dunn et al., 2012). Cells beyond the BM radius $R_\text{BM}$ experience:

$$\mathbf{F}_\text{BM} = -k_\text{BM} \cdot (r - R_\text{BM}) \cdot \hat{\mathbf{r}}$$

where $r = \|\mathbf{x} - \mathbf{x}_c\|$ is the cell's distance from the
(dynamically tracked) organoid centroid, and $k_\text{BM}$ is the **basement
membrane stiffness** (the key independent variable).

**ECM degradation** is modelled as a time-dependent increase in the BM radius
(Alcaraz et al., 2011):

$$R_\text{BM}(t) = R_\text{BM}^0 + \dot{R} \cdot t, \qquad R_\text{BM}(t) \leq R_\text{max}$$

where $\dot{R}$ is the degradation rate and $R_\text{max}$ is a biological
upper limit.

For vertex populations, the body force is distributed equally across the $n$
vertices of the element: $\mathbf{F}_{\text{node}} = \mathbf{F}_\text{BM} / n$.

**References:**
- Dunn, S.-J. et al. (2012) *J. Theor. Biol.* 298:82–91
- Alcaraz, J. et al. (2011) *EMBO J.* 30(11):2284–2295

### 3.3 Lumen Pressure Force

The organoid lumen is maintained by apical fluid secretion and
tight-junction-mediated osmotic pressure (Madan et al., 2018; Yang et al.,
2021). This is modelled as an outward radial force for cells within the
equilibrium lumen radius $R_\text{eq}$:

$$\mathbf{F}_\text{lumen} = \begin{cases} P \cdot (R_\text{eq} - r) \cdot \hat{\mathbf{r}} & \text{if } r < R_\text{eq} \\ 0 & \text{if } r \geq R_\text{eq} \end{cases}$$

where $P$ is the pressure strength. This creates a "pressure vessel" effect:
cells are pushed outward from the centre, maintaining lumen openness. Cells
that have budded outward beyond $R_\text{eq}$ feel no additional lumen pressure.

**Design choice:** $R_\text{eq}$ is set **above** the initial organoid radius
($R_\text{eq} = R_0 + 1.0$) to ensure non-zero outward pressure at $t=0$. An
earlier design with $R_\text{eq} < R_0$ resulted in zero outward force at
startup, causing the organoid to collapse inward under apical constriction.

**References:**
- Madan, E. et al. (2018) *Nature* 563:579–583
- Yang, Q. et al. (2021) *Nat. Cell Biol.* 23:733–746

### 3.4 Apical Constriction Force

Actomyosin-driven constriction of the apical cell surface, which creates
cell wedging and drives epithelial invagination (Rivas & Bhatt, 1999;
Martin et al., 2009). For cells marked as `is_apical`:

$$\mathbf{F}_\text{apical} = -\sigma_\text{ac} \cdot (A - A_\text{target}) \cdot \hat{\mathbf{r}}$$

where:
- $\sigma_\text{ac}$ is the constriction strength
- $A$ is the estimated current cell area
- $A_\text{target} = A_\text{ref} \cdot (1 - f_\text{reduction})$
- $\hat{\mathbf{r}}$ is the outward radial unit vector (force acts inward)

The force only acts when $A > A_\text{target}$ (cell needs to constrict).

**References:**
- Martin, A.C. et al. (2009) *Nature* 457:495–499
- Rivas, R.J. & Bhatt, P. (1999) in *Cell Mechanics* (eds. U. Bhatt)
- Hannezo, E. et al. (2011) *Phys. Rev. Lett.* 107:078104

---

## 4. Cell Cycle: Contact Inhibition

Proliferation is governed by the `ContactInhibitionCellCycleModel` (Chaste
built-in), which arrests the cell cycle when local cell volume fraction exceeds
a **quiescent fraction** threshold $\phi_q$:

- If $V_\text{cell} / V_\text{eq} < \phi_q$, the cell enters quiescence (G0)
- Otherwise, the cell progresses through the cycle

This implements density-dependent growth inhibition, a well-established
property of epithelial cells (Puliafito et al., 2012; Shraiman, 2005).

**References:**
- Puliafito, A. et al. (2012) *Proc. Natl. Acad. Sci.* 109(3):739–744
- Shraiman, B.I. (2005) *Proc. Natl. Acad. Sci.* 102:3318–3323

---

## 5. Cell Type Assignment

Cells are assigned proliferative types based on their angular position $\theta$
on the ring, measured from the bottom ($\theta = -\pi/2$):

| Angular region | Fraction of ring | Cell type | Proliferative |
|---------------|-----------------|-----------|---------------|
| Bottom 40% ($f < 0.2$ or $f > 0.8$) | ~32 cells (N) / ~16 cells (V) | Stem | Yes |
| Flanks 30% ($0.2 \leq f < 0.35$ or $0.65 < f \leq 0.8$) | ~24 cells (N) / ~12 cells (V) | Transit-amplifying | Yes |
| Top 30% ($0.35 \leq f \leq 0.65$) | ~24 cells (N) / ~12 cells (V) | Differentiated | No |

This mimics the crypt–villus axis:
- **Bottom** = crypt base (stem cell niche)
- **Flanks** = transit-amplifying zone
- **Top** = villus-like differentiated region

**References:**
- Barker, N. et al. (2007) *Nature* 449:1003–1007
- Snippert, H.J. et al. (2010) *Cell* 143:134–144

---

## 6. Relaxation Period (Phase 1 / Phase 2)

### Rationale

The initial cell positions are generated algorithmically (evenly spaced on a
ring with small noise). This configuration is **not mechanically equilibrated**
— forces are unbalanced and the system needs time to relax to a
quasi-equilibrium before biologically meaningful dynamics (proliferation,
budding) can occur. Without a relaxation phase, transient force artefacts from
the initial condition contaminate the early dynamics.

This two-phase approach is standard practice in Chaste simulations (Osborne
et al., 2017; see also Test3dVertexCryptOrganoid in this project).

### Implementation

**Phase 1 — Geometry relaxation** (currently 10 hours):
1. All cells are temporarily set to `DifferentiatedCellProliferativeType`
   (non-dividing)
2. All forces are active (BM, lumen pressure, apical constriction, springs)
3. `simulator.Solve()` runs for `relaxation_time` hours
4. Cells reach mechanical quasi-equilibrium without any proliferation

**Phase 2 — Growth and budding** (currently 200 hours):
1. Original proliferative types are restored (stem, TA, differentiated)
2. `simulator.SetEndTime(relaxation_time + end_time)`
3. `simulator.Solve()` continues with full dynamics

### Design note

The vertex-based test does **not** currently use a relaxation phase because
the Nagai–Honda force provides intrinsic area/perimeter restoring forces
that stabilise the geometry. This may be added if vertex simulations show
initial transient artefacts.

**References:**
- Osborne, J.M. et al. (2017) *PLoS Comput. Biol.* 13(7):e1005387

---

## 7. Experimental Sweep

### Independent variable

**ECM_STIFFNESS** — controls:
- Node-based: `basement_membrane_stiffness` (BM restoring force $k_\text{BM}$)
- Vertex-based: Nagai–Honda `deformation_energy_parameter` ($\lambda_A$)
  **and** `bm_stiffness = 0.5 × ECM_STIFFNESS`

### Levels

| ECM_STIFFNESS | Physical interpretation |
|---------------|------------------------|
| 0.5 | Very soft ECM (e.g. Matrigel) |
| 1.0 | Soft ECM |
| 2.0 | Moderate |
| 5.0 | Default / moderate–stiff |
| 10.0 | Stiff |
| 20.0 | Very stiff |
| 50.0 | Extremely stiff (rigid confinement) |

### Replicates

10 replicates per stiffness level (different random seeds), giving
**70 simulations per model** (140 total for node + vertex).

### Dependent variables (output)

- **Number of crypts**: detected from polar coordinate $r(\theta)$ profile
  using `scipy.signal.find_peaks` with prominence threshold
- **Radial statistics**: mean radius, variance, max/min, range
- **Cell counts**: total, by proliferative type

---

## 8. Post-Processing: Polar Coordinate Crypt Detection

Crypt buds are detected as **outward radial protrusions** from the organoid
surface. The algorithm (`analyse_crypt_budding.py`):

1. Compute population centroid $\mathbf{x}_c$
2. Convert all cell positions to polar coordinates $(r_i, \theta_i)$
3. Bin cells by $\theta$ and compute mean $\bar{r}(\theta)$ per bin
4. Smooth with circular padding (wraps around $2\pi$)
5. Detect peaks in $\bar{r}(\theta)$ using `scipy.signal.find_peaks`
   with minimum prominence (0.5) and minimum angular separation (30°)
6. Number of peaks = number of crypts

---

## 9. Sloughing / Boundary Conditions

Four plane-based cell killers form a bounding box at $\pm R_\text{max}$ in
$x$ and $y$ ($R_\text{max} = 5 \times R_0$). Cells crossing this boundary
are removed (sloughed), mimicking anoikis / cell shedding at the villus tip.
