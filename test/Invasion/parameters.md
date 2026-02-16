# Simulation Parameters: 2D Invasive Front

Parameters from `Test2dInvasiveFront`, `Test2dDynamicECMInvasion`, and
`Test2dPainterReplication`.

---

## 1. Domain Parameters

| Parameter | Value | Units | Source | Justified? |
|-----------|-------|-------|--------|------------|
| `domain_width` | 600.0 | µm | Metzcar et al. (2025) | Yes |
| `domain_height` | 1000.0 | µm | Metzcar et al. (2025) | Yes |
| Interaction cutoff | 50.0 | µm | Metzcar et al. (2025) | Yes; ~5 cell diameters |

---

## 2. Initial Cell Placement

| Parameter | Value | Units | Justified? |
|-----------|-------|-------|------------|
| `num_initial_cells` | 30 | — | Metzcar et al. (2025); Yes |
| $x$ range | [50, 550] | µm | Leaves 50 µm margins; design choice |
| $y$ range | [0, 50] | µm | Bottom edge; design choice |

---

## 3. Time Parameters

| Parameter | InvasiveFront | DynamicECM | PainterReplication | Units | Justified? |
|-----------|--------------|------------|-------------------|-------|------------|
| `dt` | 1.0 | 2.0 | 5.0 | minutes | NEEDS_JUSTIFICATION; stability-driven |
| `end_time` | 1440 | 7200 | 7200 | minutes | 1 day / 5 days; Painter (2009) uses 5 days |
| `sampling_multiple` | 60 | 1440 | 7200 | timesteps | Output every 1 h / 1 day |

---

## 4. ECM Contact Guidance Force Parameters

| Parameter | InvasiveFront | DynamicECM | PainterReplication | Units | Source | Justified? |
|-----------|--------------|------------|-------------------|-------|--------|------------|
| `base_speed` ($v_0$) | 1.25 | 1.25 | 2.5 | µm/min | Painter (2009) | Yes |
| `ecm_sensitivity` ($s$) | 1.0 | 1.0 | 1.0 | dimensionless | Metzcar et al. (2025) | Yes |
| `anisotropy` ($a$) | 1.0 | — | — | dimensionless | Metzcar et al. (2025) | Yes |

**Note:** PainterReplication uses 2× base speed (2.5 µm/min). NEEDS_JUSTIFICATION
for this deviation from Painter's original value.

---

## 5. ECM Grid Parameters

| Parameter | DynamicECM | PainterReplication | Units | Justified? |
|-----------|------------|-------------------|-------|------------|
| Grid spacing | 12.5 | 25.0 | µm | NEEDS_JUSTIFICATION; trades resolution vs cost |

---

## 6. ECM Dynamic Remodeling Parameters

| Parameter | Value | Units | Justified? |
|-----------|-------|-------|------------|
| `enable_remodeling` | true | — | Core feature from Metzcar et al. (2025) |
| `enable_degradation` | true | — | Core feature from Metzcar et al. (2025) |
| `enable_deposition` | false | — | Disabled for simplicity |

**Note:** Remodeling and degradation rates are internal to `DynamicECMField`
and `DynamicECMContactGuidanceForce`; defaults from those classes are used.
NEEDS_JUSTIFICATION if rates differ from Metzcar et al. (2025).

---

## 7. Cell–Cell Adhesion Parameters

| Parameter | InvasiveFront | DynamicECM | PainterReplication | Units | Justified? |
|-----------|--------------|------------|-------------------|-------|------------|
| `spring_stiffness` ($k_s$) | 0.4 | — (disabled) | — (disabled) | force/length | Metzcar et al. (2025); Yes |
| `cutoff_length` | 50.0 | — | — | µm | Metzcar et al. (2025); Yes |

**Design note:** Adhesion is commented out in `DynamicECM` and
`PainterReplication` to isolate ECM-driven migration. Re-enabling adhesion with
$k_s = 0.4$ is planned for collective invasion studies.

---

## 8. Cell Cycle Parameters (DynamicECM Only)

| Parameter | Value | Units | Justified? |
|-----------|-------|-------|------------|
| Cell cycle model | ContactInhibitionCellCycleModel | — | Standard Chaste model |
| G1 duration (stem) | 8.0 | hours | Typical epithelial; NEEDS_JUSTIFICATION |
| G1 duration (transit) | 12.0 | hours | Typical epithelial; NEEDS_JUSTIFICATION |
| S duration | 5.0 | hours | Standard (Alberts et al.) |
| G2 duration | 4.0 | hours | Standard (Alberts et al.) |
| M duration | 1.0 | hours | Standard (Alberts et al.) |
| Quiescent volume fraction ($\phi_q$) | 0.7 | dimensionless | NEEDS_JUSTIFICATION |
| Equilibrium volume ($V_\text{eq}$) | $\pi \times 5^2 \approx 78.5$ | µm² | Assumes cell radius ≈ 5 µm |
| Birth time desync | $-\mathcal{U}(0, 22)$ | hours | Prevents synchronized division |

---

## 9. Cell Influx Parameters (PainterReplication Only)

| Parameter | Value | Units | Source | Justified? |
|-----------|-------|-------|--------|------------|
| Influx interval | 180 | minutes | Painter (2009) | Yes |
| Cells per influx | 30 | — | Metzcar et al. (2025) | Yes |
| New cell $x$ range | [50, 550] | µm | Matches initial placement | Design choice |
| New cell $y$ range | [0, 50] | µm | Bottom edge | Design choice |
| New cell cycle duration | $10^6$ | minutes | Effectively infinite (no division) | Design choice |

---

## 10. Numerical Parameters

| Parameter | InvasiveFront | DynamicECM | PainterReplication | Justified? |
|-----------|--------------|------------|-------------------|------------|
| Movement threshold | 50.0 | 50.0 | 500.0 | Generous; prevents crashes |
| ECM writer interval | — | 120 timesteps | 120 timesteps | Matches VTU output |

---

## 11. Random Seed / Replicate Parameters

| Parameter | DynamicECM | PainterReplication | Justified? |
|-----------|------------|-------------------|------------|
| Base seed | 12345 | 50000 | Arbitrary; non-overlapping |
| Seed formula | base + run × 1000 | base + run × 1000 | Deterministic, reproducible |
| Replicates supported | via `RUN_NUMBER` env var | via `RUN_NUMBER` env var | — |

---

## Summary: Parameters Requiring Justification

1. **Base speed** (2.5 in PainterReplication vs 1.25 in others): inconsistent
   with Painter (2009); needs rationale or correction.
2. **ECM grid spacing** (12.5 vs 25.0 µm): no convergence study performed.
3. **dt values** (1–5 min): empirically chosen; no formal stability analysis.
4. **Cell cycle durations**: standard textbook values but not calibrated to
   specific invasive cell type.
5. **Quiescent volume fraction** (0.7): controls density-dependent arrest
   threshold; not experimentally calibrated.
6. **ECM remodeling/degradation rates**: inherited from `DynamicECMField`
   defaults; need comparison with Metzcar et al. (2025) values.
