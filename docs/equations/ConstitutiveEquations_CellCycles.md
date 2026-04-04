# Constitutive Equations — Cell Cycle Models

All cell cycle models used in the **CryptBudding** organoid simulation app, with
their mathematical formulations, phase structure, and literature references.

---

## Uniform Contact Inhibition Cell Cycle Model

**Source:**
- `projects/TissueMorphology/src/UniformContactInhibitionGenerationalCellCycleModel.hpp` **(default)**
- `projects/TissueMorphology/src/UniformContactInhibitionCellCycleModel.hpp` (simple variant)

**Used in:** All four model configurations (node2d, vertex2d, node3d, vertex3d)

**Selection:** Controlled by config toggles:
- `enableStochasticFourType = true` → **StochasticFourTypeCellCycleModel** (Montes-Olivas 2023) — **current default**
- `enableGenerationalCascade = true` → **UniformContactInhibitionGenerationalCellCycleModel** with Meineke cascade — legacy
- Both false → **UniformContactInhibitionCellCycleModel** (no generation tracking)

A custom cell cycle model that extends Chaste's `ContactInhibitionCellCycleModel` with:
1. **Stochastic total cycle duration** drawn from a uniform distribution
2. **Generational cascade** (Meineke et al. 2001) for Stem → TA → Differentiated progression (when enabled)

Cells become **quiescent** (G1-arrested) when mechanically compressed below a threshold
volume, implementing *contact inhibition of proliferation*.

### Cell Cycle Phases

([UniformContactInhibitionGenerationalCellCycleModel.hpp, `SetG1Duration`](../src/UniformContactInhibitionGenerationalCellCycleModel.hpp))

The S and G2 phases are **collapsed to zero**, leaving only M and G1:

| Phase | Condition | Duration |
|-------|-----------|----------|
| M (mitosis) | $t < T_M$ | 1.0 h (fixed) |
| G1 (gap 1) | $T_M \le t < T_{\text{cycle}}$ | $T_{\text{total}} - T_M$ (stochastic) |
| S (synthesis) | — | 0 h (collapsed) |
| G2 (gap 2) | — | 0 h (collapsed) |
| G0 (quiescence) | Differentiated cells | permanent arrest |

> **Note:** Unlike the base `ContactInhibitionCellCycleModel`, this custom model draws
> the **total cycle time** from a uniform distribution, then computes G1 as the remainder
> after subtracting the M-phase duration.

**Differentiated cells** are permanently arrested in G0 and never divide.

### Stochastic Cycle Duration

([UniformContactInhibitionGenerationalCellCycleModel.hpp, `SetG1Duration`](../src/UniformContactInhibitionGenerationalCellCycleModel.hpp))

Total cell cycle duration is drawn from a uniform distribution:

**Stem cells:**
$$
T_{\text{total}} \sim U(T_{\min}, T_{\max})
$$

**Transit-amplifying cells:**
$$
T_{\text{total}} \sim U(T_{\min}, T_{\max}) \times r_{\text{TA}}
$$

where $r_{\text{TA}}$ is the transit cycle ratio (default 1.0 = same as stem).

**G1 duration** is then computed as the remainder:
$$
T_{G_1} = \max\bigl(T_{\text{gap,min}},\; T_{\text{total}} - T_M\bigr)
$$

**Differentiated cells:**
$$
T_{G_1} = \infty \quad \text{(permanent G0 arrest)}
$$

### Contact Inhibition Rule

([ContactInhibitionCellCycleModel.cpp, `UpdateCellCyclePhase`](../../../cell_based/src/cell/cycle/ContactInhibitionCellCycleModel.cpp))

The quiescent volume threshold is:

$$
V_q = V_{\text{eq}} \times f_q
$$

At each time step $\Delta t$, if the cell is in G1 phase:

$$
\text{If } V_{\text{cell}} < V_q : \quad T_{G_1} \leftarrow T_{G_1} + \Delta t
\qquad \text{(G1 arrest — cell is quiescent)}
$$

$$
\text{If } V_{\text{cell}} \ge V_q : \quad \text{quiescent timer reset}
\qquad \text{(cell resumes cycling)}
$$

When a cell is quiescent, it is labelled with a `CellLabel` property (removed when
cycling resumes). This extends G1 indefinitely as long as the cell remains
compressed, providing negative feedback between cell density and proliferation.

### Division Rule

A cell divides when:

$$
t_{\text{age}} \ge T_M + T_{G_1}
$$

Since S and G2 phases are collapsed ($T_S = T_{G_2} = 0$), the total cell cycle time simplifies to:

$$
T_{\text{cycle}} = T_M + T_{G_1}(\text{extended by contact inhibition})
$$

With the stochastic duration drawn from $U(T_{\min}, T_{\max})$, a non-quiescent cell
will typically divide in 12–14 hours (default parameters).

### Parameters

| Symbol | Parameter | Default | CryptBudding Config |
|--------|-----------|---------|---------------------|
| $V_{\text{eq}}$ | `mEquilibriumVolume` | must be set | 1.0 (node) or target area (vertex) |
| $f_q$ | `mQuiescentVolumeFraction` | must be set | **0.7** (`quiescentFraction`) |
| $T_M$ | M-phase duration | 1.0 h | 1.0 h (fixed) |
| $T_S$ | S-phase duration | 5.0 h | **0 h** (collapsed) |
| $T_{G_2}$ | G2-phase duration | 4.0 h | **0 h** (collapsed) |
| $T_{\min}$ | `mTotalCycleMin` | 12.0 h | **12.0 h** (`stemCycleMin`) |
| $T_{\max}$ | `mTotalCycleMax` | 14.0 h | **14.0 h** (`stemCycleMax`) |
| $r_{\text{TA}}$ | `mTransitCycleRatio` | 1.0 | **1.0** (`taCycleRatio`) |
| — | `maxTransitGenerations` | 3 | **3** (generational model only) |

### Cell Proliferative Types

Three cell types are used in CryptBudding, assigned at initialisation with uniform
random probabilities:

| Type | Fraction | Behaviour |
|------|----------|-----------|
| **Stem** (`StemCellProliferativeType`) | **42%** | Full cell cycle, self-renewal (gen=0) |
| **Transit-amplifying** (`TransitCellProliferativeType`) | **49%** | Full cell cycle, limited divisions (gen 1–3) |
| **Differentiated** (`DifferentiatedCellProliferativeType`) | **9%** | Permanently in G0, no division |

The type fractions are configurable via `stemFraction` and `transitFraction`
parameters in the INI config file.

### Phase Diagram

```
                            ┌─────────┐
                            │  BIRTH  │
                            └────┬────┘
                                 │
                            ┌────▼────┐
                            │    M    │  (1h fixed)
                            └────┬────┘
                                 │
                 ┌───────────────▼───────────────┐
                 │              G1               │
                 │    (total − 1h, stochastic)   │
                 │                               │
                 │   V_cell < V_q ? ──► ARREST   │  (extended)
                 │   V_cell ≥ V_q ? ──► PROCEED  │
                 └───────────────┬───────────────┘
                                 │
                        (S and G2 collapsed)
                                 │
                            ┌────▼────┐
                            │ DIVIDE  │
                            └─────────┘
```

---

## Generational Cascade Model (Meineke et al. 2001)

**Source:** `projects/TissueMorphology/src/UniformContactInhibitionGenerationalCellCycleModel.hpp`

**Enabled by:** `enableGenerationalCascade = true` in config

When enabled, stem cell divisions follow the **Meineke generational cascade**:

### Division Rules

([UniformContactInhibitionGenerationalCellCycleModel.hpp](../src/UniformContactInhibitionGenerationalCellCycleModel.hpp))

**Stem cell division** (generation = 0):
- Parent cell: remains Stem, generation reset to 0
- Daughter cell: becomes TA, generation set to 1

**Transit-amplifying cell division** (generation 1 to `maxTransitGenerations`):
- Both parent and daughter increment their generation
- If generation > `maxTransitGenerations`: cell becomes Differentiated

**Differentiated cells**: Permanently arrested in G0, never divide.

### Generation Cascade Diagram

```
                    ┌──────────────────┐
                    │   STEM (gen=0)   │
                    │  self-renewing   │
                    └────────┬─────────┘
                             │ division
              ┌──────────────┴──────────────┐
              │                             │
       ┌──────▼──────┐              ┌───────▼──────┐
       │ STEM (gen=0)│              │  TA (gen=1)  │
       │   (parent)  │              │  (daughter)  │
       └─────────────┘              └───────┬──────┘
                                            │ division
                                     ┌──────▼──────┐
                                     │  TA (gen=2) │
                                     └──────┬──────┘
                                            │ division
                                     ┌──────▼──────┐
                                     │  TA (gen=3) │
                                     └──────┬──────┘
                                            │ division (gen > max)
                                   ┌────────▼────────┐
                                   │ DIFFERENTIATED  │
                                   │   (no division) │
                                   └─────────────────┘
```

### Parameters (Generational Model)

| Parameter | Config Key | Default | Description |
|-----------|------------|---------|-------------|
| Max TA generations | `maxTransitGenerations` | 3 | TA divisions before differentiation |
| Enable cascade | `enableGenerationalCascade` | true | Use generational vs simple model |

### References (Generational Cascade)

- Meineke, F.A., Potten, C.S. & Loeffler, M. (2001). Cell migration and organization
  in the intestinal crypt using a lattice-free model. *Cell Proliferation*, 34(4),
  253–266. doi:[10.1046/j.0960-7722.2001.00216.x](https://doi.org/10.1046/j.0960-7722.2001.00216.x)

---

## Stochastic Four-Type Cell Cycle Model (Montes-Olivas et al. 2023)

**Source:** `projects/TissueMorphology/src/StochasticFourTypeCellCycleModel.hpp` (header-only)

**Enabled by:** `enableStochasticFourType = true` in config (**current default in INI files**)

**Inherits from:** `ContactInhibitionCellCycleModel` — contact inhibition is preserved.

This model replaces the deterministic generational cascade with **stochastic, time-dependent
cell fate transitions** following the approach of Montes-Olivas et al. (2023). Four distinct
intestinal cell types are modelled via `AbstractCellMutationState` subclasses.

### Cell Types

| Cell Type | MutationState | ProliferativeType | Proliferates? |
|-----------|---------------|-------------------|---------------|
| Stem Cell (SC) | `WildTypeCellMutationState` | `TransitCellProliferativeType` | Yes |
| Transit-Amplifying (TA) | `TACellMutationState` | `TransitCellProliferativeType` | Yes |
| Paneth Cell (PC) | `PanethCellMutationState` | `DifferentiatedCellProliferativeType` | No ($T_{G_1} = \infty$) |
| Enterocyte (EC) | `EnterocyteCellMutationState` | `DifferentiatedCellProliferativeType` | No ($T_{G_1} = \infty$) |

### Division Rules

**Stem cell division** — the daughter cell's fate is drawn from:

$$
\text{daughter} = \begin{cases}
\text{SC} & \text{with probability } p_{\text{sc}\to\text{sc}} \\
\text{PC} & \text{with probability } p_{\text{sc}\to\text{pc}} \\
\text{TA} & \text{with probability } 1 - p_{\text{sc}\to\text{sc}} - p_{\text{sc}\to\text{pc}}
\end{cases}
$$

**Transit-amplifying cell division** — time-dependent probability:

$$
\text{daughter} = \begin{cases}
\text{TA} & \text{with probability } p_{\text{ta}\to\text{ta}}^{\text{early}},
  & t < t_{\text{switch}} \\
\text{EC} & \text{with probability } 1 - p_{\text{ta}\to\text{ta}}^{\text{early}},
  & t < t_{\text{switch}} \\
\text{TA} & \text{with probability } p_{\text{ta}\to\text{ta}}^{\text{late}},
  & t \ge t_{\text{switch}} \\
\text{EC} & \text{with probability } 1 - p_{\text{ta}\to\text{ta}}^{\text{late}},
  & t \ge t_{\text{switch}}
\end{cases}
$$

**Paneth cells and Enterocytes** do not divide.

### Parameters

| Symbol | Config Key | Default | Description |
|--------|------------|---------|-------------|
| $p_{\text{sc}\to\text{sc}}$ | `probStemToStem` | 0.89 | P(SC daughter = SC) |
| $p_{\text{sc}\to\text{pc}}$ | `probStemToPaneth` | 0.09 | P(SC daughter = PC) |
| $p_{\text{ta}\to\text{ta}}^{\text{early}}$ | `probTaToTaEarly` | 0.9 | P(TA daughter = TA) for days 1–4 |
| $p_{\text{ta}\to\text{ta}}^{\text{late}}$ | `probTaToTaLate` | 0.7 | P(TA daughter = TA) for days 5+ |
| $t_{\text{switch}}$ | `transitionTime` | 120.0 h | Time at which TA probabilities switch |
| — | `panethFraction` | 0.09 | Initial fraction of Paneth cells |

### Stochastic Cycle Duration

Same as the base contact inhibition model:

$$
T_{\text{total}} \sim U(T_{\min}, T_{\max}) \quad \text{(SC and TA)}
$$
$$
T_{\text{total}} = \infty \quad \text{(PC and EC)}
$$

### References (Stochastic Four-Type Model)

- Montes-Olivas, S., Marucci, L. & Homer, M. (2023). In-silico intestinal
  crypt–organoid model. *Cell Proliferation*, 56(1), e13370.
  doi:[10.1111/cpr.13370](https://doi.org/10.1111/cpr.13370)

---

## Interaction with Volume Tracking

The cell cycle model reads the `volume` item from `CellData`, which is updated at
every time step by the `VolumeTrackingModifier`. In vertex models, the
`SimpleTargetAreaModifier` (2D) or `GeometricalTargetVolumeModifier` (3D) also
updates a `target area` / target volume used for growth.

---

## Two-Phase Simulation Protocol

During the **relaxation phase** (Phase 1), all cells are temporarily set to
`DifferentiatedCellProliferativeType` → G0 arrest, ensuring no proliferation while
the initial geometry equilibrates. Original cell types are restored before the
**growth phase** (Phase 2).

---

## General References

- Shraiman, B.I. (2005). Mechanical feedback as a possible regulator of tissue growth. *Proceedings of the National Academy of Sciences*, 102(9), 3318–3323. doi:[10.1073/pnas.0404782102](https://doi.org/10.1073/pnas.0404782102)
- Nelson, C.M. & Chen, C.S. (2002). Cell-cell signaling by direct contact increases cell proliferation via a PI3K-dependent signal. *FEBS Letters*, 514(2–3), 238–242.
- Puliafito, A. et al. (2012). Collective and single cell behavior in epithelial contact inhibition. *Proceedings of the National Academy of Sciences*, 109(3), 739–744. doi:[10.1073/pnas.1007809109](https://doi.org/10.1073/pnas.1007809109)
- Chaste documentation: [ContactInhibitionCellCycleModel](https://chaste.github.io/chaste-docs/)
- Pitt-Francis, J. et al. (2009). Chaste: A test-driven approach to software development for biological modelling. *Computer Physics Communications*, 180(12), 2452–2471. doi:[10.1016/j.cpc.2009.07.019](https://doi.org/10.1016/j.cpc.2009.07.019)

---

## Target Area Modifier (2D Vertex — SimpleTargetAreaModifier)

**Source:** `Chaste/cell_based/src/simulation/modifiers/SimpleTargetAreaModifier.{hpp,cpp}`

**Used in:** vertex2d

Governs how target area grows after cell division and shrinks during apoptosis.
Works in tandem with the Nagai–Honda force (§3 in Forces doc) which penalises
deviations from the target area.

### Constitutive Equations

([SimpleTargetAreaModifier.cpp, `UpdateTargetAreaOfCell`](../../../cell_based/src/simulation/modifiers/SimpleTargetAreaModifier.cpp))

**Post-division growth** (cell age $< T_{\text{growth}}$):

$$
A_{\text{target}} = \frac{A_{\text{ref}}}{2}\left(1 + \frac{\text{age}}{T_{\text{growth}}}\right)
$$

Daughter cells begin at half the reference area and grow linearly to the full
reference area over $T_{\text{growth}}$ (defaults to G1 duration).

**Mature cell:**

$$
A_{\text{target}} = A_{\text{ref}}
$$

**At division (ready to divide):**

$$
A_{\text{target}} \to \frac{A_{\text{ref}}}{2}
$$

**Apoptosis shrinkage:**

$$
A_{\text{target}} \leftarrow A_{\text{target}} \times \left(1 - \frac{t_{\text{apoptotic}}}{2\,T_{\text{apoptosis}}}\right),
\qquad A_{\text{target}} \ge 0
$$

| Symbol | Parameter | Default | CryptBudding |
|--------|-----------|---------|--------------|
| $A_{\text{ref}}$ | `mReferenceTargetArea` | 1.0 | computed from mesh geometry |
| $T_{\text{growth}}$ | `mGrowthDuration` | G1 duration | G1 duration |

### References

- Farhadifar, R. et al. (2007). The influence of cell mechanics, cell-cell interactions, and proliferation on epithelial packing. *Current Biology*, 17(24), 2095–2104.

---

## Target Volume Modifier (3D Vertex — GeometricalTargetVolumeModifier)

**Source:** `Chaste/projects/OrganoidChaste/src/monolayer_population/forces/GeometricalTargetVolumeModifier.{hpp,cpp}`

**Used in:** vertex3d

The 3D analogue of `SimpleTargetAreaModifier`, with an additional **T1-transition
adaptation** mechanism that smoothly restores target volume after topological
rearrangements.

### Constitutive Equations

([GeometricalTargetVolumeModifier.cpp, `UpdateTargetVolumeOfCell`](../../../projects/OrganoidChaste/src/monolayer_population/forces/GeometricalTargetVolumeModifier.cpp))

**Post-division growth** (same linear law as 2D):

$$
V_{\text{target}}^{\text{growth}} = \frac{V_{\text{ref}}}{2}\left(1 + \frac{\text{age}}{T_{\text{growth}}}\right)
$$

**T1-transition volume adaptation** (time since last T1 $< T_{\text{adapt}}$):

$$
V_{\text{target}}^{T1} = V_{T1} + (V_{\text{ref}} - V_{T1})\,\frac{t - t_{T1}}{T_{\text{adapt}}}
$$

where $V_{T1}$ is the cell volume recorded at the moment of the T1 swap, and
$T_{\text{adapt}}$ is the adaptation duration.

**Resolution:** the final target volume is whichever adaptation is **slower**
(farther from equilibrium):

$$
V_{\text{target}} = \arg\max_{V \in \{V^{\text{growth}},\, V^{T1}\}} \left|V - V_{\text{ref}}\right|
$$

| Symbol | Parameter | Default | CryptBudding |
|--------|-----------|---------|--------------|
| $V_{\text{ref}}$ | `mReferenceTargetVolume` | must be set | `avgVol` (from mesh) |
| $T_{\text{growth}}$ | `mGrowthDuration` | G1 duration | 0.0 (instant) |
| $T_{\text{adapt}}$ | `mT1AdaptationDuration` | 0.0 | 0.1 h |

### References

- Hannezo, E., Prost, J. & Joanny, J.F. (2014). Theory of epithelial sheet morphology in three dimensions. *PNAS*, 111(1), 27–32.
- OrganoidChaste framework

---

## ECM Field Evolution Equations (DynamicECMField3d)

**Source:** `Chaste/projects/TissueMorphology/src/DynamicECMField3d.hpp` (header-only)

**Used in:** node3d (optional), vertex3d (optional)

The ECM is stored on a regular 3D Cartesian grid. Each voxel carries three
quantities: fibre direction $\hat{\mathbf{f}}$, density $\rho \in [0,1]$, and
anisotropy $\alpha \in [0,1]$. These evolve via cell-mediated feedback.

### Field Evolution Equations

([DynamicECMField3d.hpp](../../../projects/TissueMorphology/src/DynamicECMField3d.hpp))

**1. Degradation** (at cell position, each time step):

$$
\rho^{n+1} = \max\!\bigl(0,\; \rho^n - k_{\text{deg}}\,\Delta t\bigr)
$$

If $\rho < 0.5$, fibre coherence decays: $\alpha \leftarrow 0.99\,\alpha$.

**2. Fibre remodelling** (traction-driven, at cell position):

$$
\hat{\mathbf{f}}^{\,n+1} = \frac{\hat{\mathbf{f}}^{\,n}
+ k_{\text{remodel}}\,\min\!\left(1,\,\tfrac{|\mathbf{F}|}{F_{\text{ref}}}\right)
\hat{\mathbf{p}}\,\Delta t}
{\left|\hat{\mathbf{f}}^{\,n}
+ k_{\text{remodel}}\,\min\!\left(1,\,\tfrac{|\mathbf{F}|}{F_{\text{ref}}}\right)
\hat{\mathbf{p}}\,\Delta t\right|}
$$

where $\hat{\mathbf{p}}$ is the component of the cell traction direction
$\hat{\mathbf{d}}$ perpendicular to the current fibre:

$$
\hat{\mathbf{p}} = \frac{\hat{\mathbf{d}} - (\hat{\mathbf{d}} \cdot \hat{\mathbf{f}}^{\,n})\,\hat{\mathbf{f}}^{\,n}}{|\hat{\mathbf{d}} - (\hat{\mathbf{d}} \cdot \hat{\mathbf{f}}^{\,n})\,\hat{\mathbf{f}}^{\,n}|}
$$

Anisotropy increases: $\alpha \leftarrow \min(1,\; \alpha + 0.01\,\Delta t)$.

**3. Deposition** (at cell position, along migration direction):

$$
\rho^{n+1} = \min\!\bigl(1,\; \rho^n + k_{\text{dep}}\,\Delta t\bigr)
$$

Fibre direction blended toward deposition direction with weight:

$$
w = \frac{k_{\text{dep}}\,\Delta t}{\rho^n + k_{\text{dep}}\,\Delta t},
\qquad
\hat{\mathbf{f}}^{\,n+1} \propto (1-w)\,\hat{\mathbf{f}}^{\,n} + w\,\hat{\mathbf{d}}_{\text{dep}}
$$

**4. Diffusion** (6-neighbour discrete Laplacian on grid):

$$
\rho^{n+1}_{ijk} = \rho^n_{ijk} + D\,\Delta t\,\left(\bar{\rho}_{ijk} - \rho^n_{ijk}\right)
$$

where $\bar{\rho}_{ijk}$ is the mean density over the voxel and its 6
face-adjacent neighbours. The fibre direction is smoothed similarly.

| Symbol | Parameter | Default |
|--------|-----------|---------|
| $k_{\text{deg}}$ | Degradation rate | 0.002 |
| $k_{\text{remodel}}$ | Remodelling rate | 0.05 |
| $k_{\text{dep}}$ | Deposition rate | 0.0003 |
| $D$ | Diffusion coefficient | per-timestep smoothing |
| Grid spacing | `ecmGridSpacing` | 10.0 |
| Domain half-width | `ecmDomainHalf` | 80.0 |

### References

- Fraley, S.I. et al. (2015). Three-dimensional matrix fiber alignment modulates cell migration and MT1-MMP utility. *Scientific Reports*, 5, 14580.
- Metzcar, J. et al. (2019). A review of cell-based computational modeling in cancer biology. *JCO Clinical Cancer Informatics*, 3, 1–13.
