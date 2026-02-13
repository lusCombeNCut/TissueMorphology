# Crypt Budding 2D — Simulation Parameters

All parameters for the 2D crypt budding stiffness-sweep simulations.
Both models share the same experimental design: sweep over ECM stiffness
values (env var `ECM_STIFFNESS`) with multiple replicates (`RUN_NUMBER`).

---

## 1. Node-Based Model (`Test2dCryptBuddingNodeBased.hpp`)

### 1.1 Domain Geometry

| Parameter | Value | Justification |
|-----------|-------|---------------|
| `cells_across` | 40 | **(NEEDS JUSTIFICATION)** — chosen to produce a wide enough monolayer for multiple crypts to form side-by-side |
| `cells_up` | 12 | **(NEEDS JUSTIFICATION)** — provides sufficient vertical extent for crypt invaginations |
| `interaction_cutoff` | 3.0 cell diameters | Standard Chaste default for node-based models; ~3× rest length ensures neighbour connectivity (Meineke et al., 2001) |

### 1.2 Time

| Parameter | Value | Justification |
|-----------|-------|---------------|
| `dt` | 0.005 h | Standard for node-based Chaste simulations; small enough for numerical stability with spring stiffness 30 (CFL-like constraint) |
| `end_time` | 200.0 h (~8.3 days) | Matches typical intestinal organoid culture durations of 7–10 days (Sato et al., 2009) |
| `sampling_multiple` | 200 (= 1 h output interval) | Hourly snapshots; balances disk usage with temporal resolution |

### 1.3 Basement Membrane / Confinement

| Parameter | Value | Justification |
|-----------|-------|---------------|
| `bm_radius` | 30.0 | **(NEEDS JUSTIFICATION)** — initial confinement radius scaled to the enlarged mesh |
| `bmStiffness` | Swept: {0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0} | Key independent variable; range spans soft Matrigel (~0.5 kPa) to stiff collagen (~50 kPa) in dimensionless units (Gjorevski et al., 2016) |
| `ecm_degradation_rate` | 0.05 radius units/h | **(NEEDS JUSTIFICATION)** — allows BM to loosen gradually, preventing growth arrest from fixed confinement |
| `ecm_max_radius` | 60.0 | **(NEEDS JUSTIFICATION)** — biological upper bound to prevent runaway expansion |

### 1.4 Spring Force (GeneralisedLinearSpringForce)

| Parameter | Value | Justification |
|-----------|-------|---------------|
| `spring_stiffness` (Meineke) | 30.0 | Standard Chaste value for epithelial monolayers (Meineke et al., 2001; Van Leeuwen et al., 2009) |
| `spring_cutoff` | 1.5 | Natural rest length; cells beyond this distance do not interact. Standard for Chaste node-based models |
| `division_resting_spring_length` | 0.5 | Daughter cells start at half rest length and grow apart. Chaste default (Meineke et al., 2001) |
| `spring_growth_duration` | 1.0 h | Time for daughter-cell spring to reach full rest length. Chaste default |

### 1.5 Cell Cycle (ContactInhibitionCellCycleModel)

| Parameter | Value | Justification |
|-----------|-------|---------------|
| `quiescent_volume_fraction` | 0.6 | **(NEEDS JUSTIFICATION)** — lowered from default 0.9 to allow proliferation under moderate compression; prevents premature growth arrest |
| `equilibrium_volume` | 1.0 | Normalised cell volume; standard Chaste convention |
| `stem_g1_min` | 12.0 h | **(NEEDS JUSTIFICATION)** — approximate intestinal stem cell cycle timing |
| `stem_g1_max` | 14.0 h | **(NEEDS JUSTIFICATION)** — gives ~24 h total cell cycle with S+G2+M phases, matching murine intestinal stem cell division rate (Snippert et al., 2010) |
| `ta_g1_min` | 4.0 h | **(NEEDS JUSTIFICATION)** — transit-amplifying cells divide faster than stem cells |
| `ta_g1_max` | 6.0 h | **(NEEDS JUSTIFICATION)** — gives ~12–16 h total TA cycle, consistent with rapid crypt base turnover |
| Birth time desynchronisation | Uniform on [−12, 0] h | Prevents artificial synchronous division waves at simulation start. Standard practice |

### 1.6 Cell Type Assignment

| Parameter | Value | Justification |
|-----------|-------|---------------|
| Stem cell region | Bottom 20% of y-range | Reflects crypt base stem cell niche; Lgr5+ cells reside at crypt bottom (Barker et al., 2007) |
| Transit-amplifying region | 20%–60% of y-range | **(NEEDS JUSTIFICATION)** — TA compartment between stem niche and differentiated zone |
| Differentiated region | Top 40% of y-range | **(NEEDS JUSTIFICATION)** — post-mitotic cells that migrate upward and slough |

### 1.7 Sloughing (PlaneBasedCellKiller)

| Parameter | Value | Justification |
|-----------|-------|---------------|
| `slough_height` | 40.0 | **(NEEDS JUSTIFICATION)** — y-coordinate above which cells are removed; scaled to the enlarged domain |

### 1.8 Boundary Conditions

| Parameter | Value | Justification |
|-----------|-------|---------------|
| Bottom plane BC | y ≥ 0 enforced (normal = −y) | Prevents cells from falling below the substrate; models rigid growth surface |

### 1.9 Cell Population

| Parameter | Value | Justification |
|-----------|-------|---------------|
| `absolute_movement_threshold` | 50.0 | Generous threshold to avoid spurious "cell moved too far" exceptions in large-domain simulations. **(NEEDS JUSTIFICATION)** |

### 1.10 Sweep / Reproducibility

| Parameter | Value | Justification |
|-----------|-------|---------------|
| Default `ECM_STIFFNESS` | 5.0 | Mid-range value for local testing |
| Stiffness sweep values | {0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0} | Defined in HPC submit script; logarithmically spaced to cover soft–stiff substrate range |
| Replicates per stiffness | 10 | Provides statistical power for mean ± SD of crypt counts |
| RNG seed formula | `stiffness × 10000 + run × 137` | Deterministic and unique per (stiffness, run) pair |

---

## 2. Vertex-Based Model (`Test2dCryptBuddingVertexBased.hpp`)

### 2.1 Domain Geometry

| Parameter | Value | Justification |
|-----------|-------|---------------|
| `cells_across` | 28 | **(NEEDS JUSTIFICATION)** — width of honeycomb vertex mesh; produces an epithelial sheet wide enough for multiple invaginations |
| `cells_up` | 16 | **(NEEDS JUSTIFICATION)** — sufficient vertical rows for buckling instability |
| `flat_bottom` | true | Flat initial geometry; buckling emerges from proliferative instability rather than initial curvature |

### 2.2 Time

| Parameter | Value | Justification |
|-----------|-------|---------------|
| `dt` | 0.002 h | Smaller than node-based model; vertex models require tighter timesteps for stability due to T1/T2 transitions (Nagai & Honda, 2001) |
| `end_time` | 200.0 h (~8.3 days) | Matches node-based model and organoid culture duration |
| `sampling_multiple` | 500 (= 1 h output interval) | Hourly snapshots; 500 × 0.002 = 1.0 h |

### 2.3 Nagai-Honda Force

| Parameter | Value | Justification |
|-----------|-------|---------------|
| `nagai_honda_deformation` | = `ecmStiffness` (swept) | **Key parameter**: area elasticity models substrate stiffness. Higher values resist cell deformation → suppress buckling (Hannezo et al., 2011) |
| `nagai_honda_membrane` | 10.0 | **(NEEDS JUSTIFICATION)** — perimeter/membrane contractility energy |
| `nagai_honda_adhesion` | 1.0 | **(NEEDS JUSTIFICATION)** — cell-cell adhesion energy |
| `nagai_honda_boundary_adhesion` | 2.0 | **(NEEDS JUSTIFICATION)** — boundary adhesion, slightly higher than cell-cell to maintain monolayer integrity |

### 2.4 Basement Membrane / Confinement

| Parameter | Value | Justification |
|-----------|-------|---------------|
| `bm_stiffness` | `ecmStiffness × 0.5` | **(NEEDS JUSTIFICATION)** — BM force co-varies with ECM stiffness at half strength to avoid double-counting with Nagai-Honda area term |
| `bm_radius` | 24.0 | **(NEEDS JUSTIFICATION)** — initial confinement radius scaled to vertex mesh dimensions |
| `ecm_degradation_rate` | 0.05 radius units/h | Matches node-based model |
| `ecm_max_radius` | 50.0 | **(NEEDS JUSTIFICATION)** — slightly smaller ceiling than node model due to smaller mesh |

### 2.5 Vertex Mesh Topology

| Parameter | Value | Justification |
|-----------|-------|---------------|
| `cell_rearrangement_threshold` (T1) | 0.1 | Standard Chaste default; edge length below which T1 swap occurs (Nagai & Honda, 2001) |
| `T2_threshold` | 0.01 | Standard Chaste default; element area below which T2 (cell removal) transition occurs |
| `cell_rearrangement_ratio` | 1.5 | Standard Chaste default; new edge length after T1 swap = ratio × threshold |

### 2.6 Cell Cycle (ContactInhibitionCellCycleModel)

| Parameter | Value | Justification |
|-----------|-------|---------------|
| `quiescent_volume_fraction` | 0.6 | Matches node-based model; **(NEEDS JUSTIFICATION)** |
| `equilibrium_volume` | 1.0 (= `target_area`) | Normalised; matches target area modifier |

### 2.7 Target Area

| Parameter | Value | Justification |
|-----------|-------|---------------|
| `target_area` | 1.0 | Normalised cell area; standard for Chaste vertex models |
| `SimpleTargetAreaModifier` reference area | 1.0 | Maintains homeostatic cell area; standard Chaste modifier |

### 2.8 Cell Type Assignment

| Parameter | Value | Justification |
|-----------|-------|---------------|
| Stem cell region | Bottom 25% of y-range | Slightly wider stem niche than node model; **(NEEDS JUSTIFICATION)** — may reflect vertex model needing more proliferative cells to drive buckling |
| Transit-amplifying region | 25%–60% of y-range | **(NEEDS JUSTIFICATION)** |
| Differentiated region | Top 40% of y-range | **(NEEDS JUSTIFICATION)** |

### 2.9 Sloughing (PlaneBasedCellKiller)

| Parameter | Value | Justification |
|-----------|-------|---------------|
| `slough_height` | 28.0 | **(NEEDS JUSTIFICATION)** — matches vertex mesh vertical extent |

### 2.10 Boundary Conditions

| Parameter | Value | Justification |
|-----------|-------|---------------|
| Bottom plane BC | y ≥ 0 enforced (normal = −y) | Same as node model |
| `UseJiggledNodesOnPlane` | true | Adds small random perturbation to nodes on the boundary plane. Prevents degenerate vertex configurations at the pinned boundary |

### 2.11 Sweep / Reproducibility

| Parameter | Value | Justification |
|-----------|-------|---------------|
| Default `ECM_STIFFNESS` | 5.0 | Mid-range value for local testing |
| Stiffness sweep values | {0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0} | Same as node model |
| Replicates per stiffness | 10 | Same as node model |
| RNG seed formula | `stiffness × 10000 + run × 137` | Same as node model |

---

## 3. Parameter Comparison (Node vs Vertex)

| Parameter | Node-based | Vertex-based | Notes |
|-----------|-----------|-------------|-------|
| Initial cells | 40 × 12 = 480 | 28 × 16 = 448 | Similar cell count |
| `dt` | 0.005 h | 0.002 h | Vertex needs smaller dt |
| Cell-cell mechanics | Linear spring (k=30) | Nagai-Honda (area + perimeter) | Different force models |
| BM stiffness | = `ECM_STIFFNESS` | = `ECM_STIFFNESS × 0.5` | Halved in vertex to avoid double-counting with area elasticity |
| BM radius | 30.0 | 24.0 | Scaled to respective meshes |
| ECM max radius | 60.0 | 50.0 | Scaled to respective meshes |
| Slough height | 40.0 | 28.0 | Scaled to respective meshes |
| Quiescent fraction | 0.6 | 0.6 | Matched |
| Stem region | Bottom 20% | Bottom 25% | Slight difference |
| Bottom BC jiggle | No | Yes | Vertex needs perturbation to avoid degenerate edges |
