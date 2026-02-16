# Simulation Parameters: 2D Crypt Budding

All parameters from both the **node-based** and **vertex-based** test scripts.
Parameters marked `NEEDS_JUSTIFICATION` have not yet been rigorously derived
from experimental data or theoretical analysis and require literature support
or sensitivity analysis.

---
## 1. Geometry Parameters

| Parameter | Node-based | Vertex-based | Units | Justified? |
|-----------|-----------|-------------|-------|------------|
| `organoid_radius` | 8.0 | 7.0 (midshell) | cell diameters | NEEDS_JUSTIFICATION |
| `inner_radius` | — | 6.0 | cell diameters | NEEDS_JUSTIFICATION |
| `outer_radius` | — | 8.0 | cell diameters | NEEDS_JUSTIFICATION |
| `num_cells_in_ring` | 80 | 40 | — | NEEDS_JUSTIFICATION |
| `center_x`, `center_y` | 0.0 | 0.0 | cell diameters | Arbitrary (origin) |
| Radial noise amplitude | 0.3 (uniform) | 0.2 (uniform) | cell diameters | NEEDS_JUSTIFICATION |

**Notes:**
- Node-based uses a single-radius ring; vertex uses a finite-thickness annulus.
- The vertex `organoid_radius` is the midshell: $0.5 \times (R_\text{inner} + R_\text{outer})$.
- Cell counts differ because vertex elements are larger (finite-area quads).

---

## 2. Time Parameters

| Parameter | Node-based | Vertex-based | Units | Justified? |
|-----------|-----------|-------------|-------|------------|
| `dt` | 0.005 | Adaptive (see below) | hours | Stability-driven (CFL-like); NEEDS_JUSTIFICATION for absolute scale |
| `relaxation_time` | 10.0 | — (not used) | hours | NEEDS_JUSTIFICATION |
| `end_time` | 200.0 | 200.0 | hours | ~8 days; see Sato et al. (2009) for organoid maturation timescale |
| `sampling_multiple` | 200 | $1/\text{dt}$ | timesteps | Outputs every ~1 hour |

### Vertex dt (adaptive by stiffness)

| ECM_STIFFNESS | dt | Rationale |
|---------------|-----|-----------|
| < 1.0 | 0.0002 | Softer mesh deforms more → smaller dt for stability |
| 1.0–2.0 | 0.0005 | — |
| 2.0–5.0 | 0.001 | — |
| ≥ 5.0 | 0.001 | Stiffer elements are more numerically stable |

**Note:** These dt values were determined empirically to avoid vertex mesh
crashes (T1/T2 swap failures). NEEDS_JUSTIFICATION — a formal CFL analysis
would be preferable.

---

## 3. Basement Membrane / ECM Parameters

| Parameter | Node-based | Vertex-based | Units | Justified? |
|-----------|-----------|-------------|-------|------------|
| `bmStiffness` / `ecmStiffness` | **SWEEP VARIABLE** | **SWEEP VARIABLE** | force/length | Key independent variable |
| `bm_radius` | $R_0 + 2.0 = 10.0$ | $R_\text{outer} + 2.0 = 10.0$ | cell diameters | NEEDS_JUSTIFICATION |
| `bm_stiffness` (vertex) | — | $0.5 \times \text{ECM\_STIFFNESS}$ | force/length | NEEDS_JUSTIFICATION (scaling factor 0.5 is ad hoc) |
| `ecm_degradation_rate` | 0.02 | 0.02 | cell diameters/hour | NEEDS_JUSTIFICATION |
| `ecm_max_radius` | $4 \times R_0 = 32.0$ | $4 \times R_\text{outer} = 32.0$ | cell diameters | NEEDS_JUSTIFICATION |

**Stiffness sweep levels:**
`[0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]` — 7 levels × 10 replicates = 70 jobs per model.

**Notes:**
- `bm_radius` gap of 2.0 above the cell ring allows slight initial expansion before BM force activates.
- The 0.5 scaling factor for vertex `bm_stiffness` accounts for the fact that Nagai–Honda deformation energy already provides confinement; the BM force is supplementary.

---

## 4. Lumen Pressure Parameters

| Parameter | Node-based | Vertex-based | Units | Justified? |
|-----------|-----------|-------------|-------|------------|
| `lumen_pressure` | 2.0 | 2.0 | force/length | NEEDS_JUSTIFICATION |
| `lumen_eq_radius` | $R_0 + 1.0 = 9.0$ | $R_\text{outer} + 1.0 = 9.0$ | cell diameters | Design choice (see below) |
| `track_center` | true | true | — | Auto-tracks centroid |

**Design choice for `lumen_eq_radius`:**
Set **above** the initial cell radius so outward pressure is non-zero at $t = 0$.
An earlier choice of $R_\text{eq} = R_0 - 1.5$ resulted in zero outward pressure
at startup (all cells beyond $R_\text{eq}$), causing inward collapse under
apical constriction. The current value ensures:

$$F_\text{lumen}(t=0) = P \times (R_\text{eq} - R_0) = 2.0 \times 1.0 = 2.0 \text{ (outward)}$$

---

## 5. Apical Constriction Parameters

| Parameter | Node-based | Vertex-based | Units | Justified? |
|-----------|-----------|-------------|-------|------------|
| `apical_constriction_strength` ($\sigma_\text{ac}$) | 3.0 | 3.0 | force/area | NEEDS_JUSTIFICATION |
| `target_reduction` ($f_\text{reduction}$) | 0.5 (default) | 0.5 (default) | — | NEEDS_JUSTIFICATION; Martin et al. (2009) report 50–70% apical area reduction in *Drosophila* |
| `is_apical` marker | 1.0 (all cells) | 1.0 (all cells) | — | All cells face the lumen in this geometry |

---

## 6. Spring / Adhesion Parameters (Node-Based Only)

| Parameter | Value | Units | Justified? |
|-----------|-------|-------|------------|
| `spring_stiffness` ($k_s$) | 30.0 | force/length | NEEDS_JUSTIFICATION |
| `spring_cutoff` | 1.5 | cell diameters | Standard Chaste default (Meineke et al., 2001) |
| `division_resting_spring_length` | 0.5 | cell diameters | Standard Chaste default |
| `spring_growth_duration` | 1.0 | hours | Standard Chaste default |
| `apical_apical_adhesion` ($\mu_{AA}$) | 1.2 | multiplier | NEEDS_JUSTIFICATION |
| `basal_basal_adhesion` ($\mu_{BB}$) | 1.0 | multiplier | Baseline (unity) |
| `apical_basal_adhesion` ($\mu_{AB}$) | 0.5 | multiplier | NEEDS_JUSTIFICATION; should be weaker than homotypic per DAH |
| `interaction_cutoff` | 2.5 | cell diameters | Standard Chaste practice for NodesOnlyMesh |

---

## 7. Nagai–Honda Parameters (Vertex-Based Only)

| Parameter | Value / Formula | Units | Justified? |
|-----------|----------------|-------|------------|
| `nagai_honda_deformation` ($\lambda_A$) | = ECM_STIFFNESS | force/area | KEY VARIABLE; this couples ECM stiffness to vertex mechanics |
| `nagai_honda_membrane` ($\lambda_P$) | 10.0 | force/length | NEEDS_JUSTIFICATION |
| `nagai_honda_adhesion` ($\gamma$) | 1.0 | force/length | NEEDS_JUSTIFICATION |
| `nagai_honda_boundary_adhesion` | 2.0 | force/length | NEEDS_JUSTIFICATION |
| `target_area` ($A_0$) | $\frac{\Delta\theta}{2}(R_\text{outer}^2 - R_\text{inner}^2)$ | cell diameters² | Computed from geometry |
| `t1_threshold` | 0.15 (soft) / 0.1 (stiff) | cell diameters | NEEDS_JUSTIFICATION; empirically chosen for stability |
| `t2_threshold` | 0.01 | cell diameters² | Standard Chaste default |
| Rosette formation probability | 0.0 | — | Disabled (rosettes cause crashes in annular geometry) |

---

## 8. Cell Cycle Parameters

| Parameter | Node-based | Vertex-based | Units | Justified? |
|-----------|-----------|-------------|-------|------------|
| `quiescent_fraction` ($\phi_q$) | 0.7 | 0.6 | — | NEEDS_JUSTIFICATION; relates to contact inhibition threshold |
| `equilibrium_volume` ($V_\text{eq}$) | 1.0 | = `target_area` | cell diameters^DIM | Standard normalisation |
| Birth time (desync) | $-\mathcal{U}(0, 12)$ | $-\mathcal{U}(0, 12)$ | hours | Prevents synchronised division; 12h ≈ half cell cycle |

**Notes:**
- `quiescent_fraction` differs between models (0.7 vs 0.6). The vertex model
  uses a lower threshold because vertex cells have more precise area tracking.
  Both values NEED_JUSTIFICATION.

---

## 9. Cell Type Assignment

| Parameter | Value | Justified? |
|-----------|-------|------------|
| Stem cell arc | Bottom 40% of ring | NEEDS_JUSTIFICATION; based on crypt-villus axis anatomy |
| Transit-amplifying arc | Flanks 30% | NEEDS_JUSTIFICATION |
| Differentiated arc | Top 30% | NEEDS_JUSTIFICATION |

---

## 10. Boundary Conditions

| Parameter | Node-based | Vertex-based | Units | Justified? |
|-----------|-----------|-------------|-------|------------|
| `max_radius_for_slough` | $5 \times R_0 = 40.0$ | $5 \times R_\text{outer} = 40.0$ | cell diameters | NEEDS_JUSTIFICATION; far enough to not interfere with budding |
| Movement threshold | 50.0 | — | cell diameters/step | Generous; prevents crash on large displacement |

---

## 11. Simulation Infrastructure

| Parameter | Node-based | Vertex-based | Units | Justified? |
|-----------|-----------|-------------|-------|------------|
| Random seed formula | $\text{stiffness} \times 10000 + \text{run} \times 137$ | Same | — | Deterministic, non-overlapping seeds |
| Number of replicates | 10 | 10 | — | NEEDS_JUSTIFICATION; standard for stochastic simulations |

---

## 12. Post-Processing Parameters (`analyse_crypt_budding.py`)

| Parameter | Value | Justified? |
|-----------|-------|------------|
| Peak prominence threshold | 0.5 cell diameters | NEEDS_JUSTIFICATION |
| Minimum angular separation | 30° (~0.52 rad) | NEEDS_JUSTIFICATION; prevents counting noise as separate crypts |
| Number of angular bins | 72 (5° per bin) | NEEDS_JUSTIFICATION |

---

## Summary: Parameters Requiring Justification

The following parameters are set to reasonable but unvalidated values and
require either literature-based calibration or sensitivity analysis:

1. **Geometry**: organoid_radius, inner/outer radius, num_cells
2. **Forces**: spring_stiffness, lumen_pressure, apical_constriction_strength
3. **Adhesion**: apical-apical, apical-basal multipliers
4. **Nagai–Honda**: membrane energy, adhesion parameters, boundary adhesion
5. **Cell cycle**: quiescent_fraction (both models)
6. **ECM**: degradation rate, max radius, BM radius gap
7. **Time**: dt (absolute scale), relaxation_time
8. **Detection**: peak prominence, angular separation thresholds
9. **Sweep**: number of replicates, stiffness levels
