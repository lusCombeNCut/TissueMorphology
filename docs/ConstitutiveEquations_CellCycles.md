# Constitutive Equations — Cell Cycle Models

All cell cycle models used in the **CryptBudding** organoid simulation app, with
their mathematical formulations, phase structure, and literature references.

---

## Contact Inhibition Cell Cycle Model

**Source:** `Chaste/cell_based/src/cell/cycle/ContactInhibitionCellCycleModel.{hpp,cpp}`

**Used in:** All four model configurations (node2d, vertex2d, node3d, vertex3d)

A volume-dependent cell cycle model in which cells become **quiescent** (G1-arrested)
when they are mechanically compressed below a threshold volume. This implements
*contact inhibition of proliferation* — a cornerstone mechanism of epithelial
homeostasis.

### Cell Cycle Phases

([AbstractSimplePhaseBasedCellCycleModel.cpp, `SetG1Duration`](../../../cell_based/src/cell/cycle/AbstractSimplePhaseBasedCellCycleModel.cpp);
[AbstractPhaseBasedCellCycleModel.cpp, constructor](../../../cell_based/src/cell/cycle/AbstractPhaseBasedCellCycleModel.cpp))

The cell cycle follows the standard $M \to G_1 \to S \to G_2 \to M$ sequence.
Phase assignment is based on cell age $t$ since birth:

| Phase | Condition | Default Duration |
|-------|-----------|------------------|
| M (mitosis) | $t < T_M$ | 1.0 h |
| G1 (gap 1) | $T_M \le t < T_M + T_{G_1}$ | 14.0 h (stem) / 2.0 h (transit) — deterministic |
| S (synthesis) | $T_M + T_{G_1} \le t < T_M + T_{G_1} + T_S$ | 5.0 h |
| G2 (gap 2) | $T_M + T_{G_1} + T_S \le t < T_{\text{cycle}}$ | 4.0 h |
| G0 (quiescence) | Differentiated cells | permanent arrest |

> **Note:** `ContactInhibitionCellCycleModel` uses `AbstractSimplePhaseBasedCellCycleModel`,
> which assigns G1 duration **deterministically** (not stochastically). The base G1
> is exactly 14.0 h for stem cells and 2.0 h for transit cells. These are then
> extended by contact inhibition when the cell is compressed.

**Differentiated cells** are permanently arrested in G0 and never divide.

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
t_{\text{age}} \ge T_M + T_{G_1} + T_S + T_{G_2}
$$

The total cell cycle time is therefore:

$$
T_{\text{cycle}} = T_M + T_{G_1}(\text{extended by contact inhibition}) + T_S + T_{G_2}
$$

### Parameters

| Symbol | Parameter | Default | CryptBudding |
|--------|-----------|---------|--------------|
| $V_{\text{eq}}$ | `mEquilibriumVolume` | must be set | 1.0 (node) or target area (vertex) |
| $f_q$ | `mQuiescentVolumeFraction` | must be set | 0.7 |
| $T_M$ | M-phase duration | 1.0 h | 1.0 h |
| $T_S$ | S-phase duration | 5.0 h | 5.0 h |
| $T_{G_2}$ | G2-phase duration | 4.0 h | 4.0 h |
| $T_{G_1}$ | G1 duration (base, deterministic) | 14.0 h (stem) / 2.0 h (transit) | same |

### Cell Proliferative Types

Three cell types are used in CryptBudding, assigned at initialisation with uniform
random probabilities:

| Type | Fraction | Behaviour |
|------|----------|-----------|
| **Stem** (`StemCellProliferativeType`) | 20% | Full cell cycle ($T_{\text{cycle}}$ = 24.0 h base), self-renewal |
| **Transit-amplifying** (`TransitCellProliferativeType`) | 50% | Full cell cycle ($T_{\text{cycle}}$ = 12.0 h base), limited divisions |
| **Differentiated** (`DifferentiatedCellProliferativeType`) | 30% | Permanently in G0, no division |

The type fractions are configurable via `stemFraction` and `transitFraction`
parameters.

### Phase Diagram

```
                            ┌─────────┐
                            │  BIRTH  │
                            └────┬────┘
                                 │
                            ┌────▼────┐
                            │    M    │  (1h)
                            └────┬────┘
                                 │
                 ┌───────────────▼───────────────┐
                 │              G1               │
                 │                               │
                 │   V_cell < V_q ? ──► ARREST   │  (extended)
                 │   V_cell ≥ V_q ? ──► PROCEED  │
                 └───────────────┬───────────────┘
                                 │
                            ┌────▼────┐
                            │    S    │  (5h)
                            └────┬────┘
                                 │
                            ┌────▼────┐
                            │   G2    │  (4h)
                            └────┬────┘
                                 │
                            ┌────▼────┐
                            │ DIVIDE  │
                            └─────────┘
```

### Interaction with Volume Tracking

The cell cycle model reads the `volume` item from `CellData`, which is updated at
every time step by the `VolumeTrackingModifier`. In vertex models, the
`SimpleTargetAreaModifier` (2D) or `GeometricalTargetVolumeModifier` (3D) also
updates a `target area` / target volume used for growth.

### Two-Phase Simulation Protocol

During the **relaxation phase** (Phase 1), all cells are temporarily set to
`DifferentiatedCellProliferativeType` → G0 arrest, ensuring no proliferation while
the initial geometry equilibrates. Original cell types are restored before the
**growth phase** (Phase 2).

### References

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
