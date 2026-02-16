# Simulation Parameters: ECM Force Validation Tests

---

## 1. Geometry Parameters

### Test2dApicalConstriction

| Parameter | Value | Units | Justified? |
|-----------|-------|-------|------------|
| `grid_width` | 10 | cells | Arbitrary; sufficient for invagination |
| `grid_height` | 10 | cells | Arbitrary |
| `cell_spacing` | 1.0 | cell diameters | Standard Chaste unit spacing |
| Apical region radius | 2.0 | cell diameters | NEEDS_JUSTIFICATION |
| Interaction cutoff | 1.5 | cell diameters | Standard Chaste default |

### Test2dThresholdECMForce (2D)

| Parameter | Value | Units | Justified? |
|-----------|-------|-------|------------|
| `num_cells` | 50 | — | Moderate count for 2D |
| `organoid_radius` | 5.0 | µm | NEEDS_JUSTIFICATION |
| Inner/outer cell boundary | 3.0 | µm | NEEDS_JUSTIFICATION |
| Interaction cutoff | 1.5 | cell diameters | Standard Chaste default |

### TestThresholdECMForce (3D)

| Parameter | Value | Units | Justified? |
|-----------|-------|-------|------------|
| `num_cells` | 50 | — | Moderate count for 3D |
| `organoid_radius` | 18.0 | µm | NEEDS_JUSTIFICATION |
| Inner/outer cell boundary | 12.0 | µm | NEEDS_JUSTIFICATION |
| Interaction cutoff | 5.0 | µm | Larger for 3D; NEEDS_JUSTIFICATION |

---

## 2. Time Parameters

| Parameter | Apical Constriction | 2D Threshold | 3D Threshold | Units |
|-----------|-------------------|-------------|-------------|-------|
| `dt` | 0.01 | 0.005 | 0.005 | hours (Chaste time units) |
| `end_time` | 10.0 | 20.0 | 50.0 | hours |
| `sampling_multiple` | 10 | 40 | 40 | timesteps |
| Output interval | 0.1 | 0.2 | 0.2 | hours |

**Justification:** dt values chosen empirically for numerical stability. Smaller
dt for threshold tests due to stiff BM restoring forces. All NEED_JUSTIFICATION
via formal CFL analysis.

---

## 3. Cell–Cell Spring Parameters

| Parameter | Apical Constriction | 2D/3D Threshold | Units | Justified? |
|-----------|-------------------|----------------|-------|------------|
| `spring_stiffness` ($k_s$) | 30.0 | 35.0 | force/length | NEEDS_JUSTIFICATION |
| `division_resting_spring_length` | — | — | cell diameters | Default (0.5) |
| `spring_growth_duration` | — | — | hours | Default (1.0) |

---

## 4. Apical Constriction Parameters

| Parameter | Value | Units | Justified? |
|-----------|-------|-------|------------|
| `constriction_strength` ($\sigma_\text{ac}$) | 10.0 | force/area | NEEDS_JUSTIFICATION; Martin et al. (2009) report actomyosin-driven constriction but no direct force calibration |
| `target_reduction` ($f_\text{reduction}$) | 0.6 | dimensionless | Partially justified; Martin et al. (2009) report 50–70% apical area reduction in *Drosophila* ventral furrow |
| `is_apical` marker | 1.0 for $r < 2.0$ | — | Geometric assignment |

---

## 5. Basement Membrane / ECM Parameters

| Parameter | 2D Threshold | 3D Threshold | Units | Justified? |
|-----------|-------------|-------------|-------|------------|
| `bm_stiffness` ($k_\text{BM}$) | 2.0 | 2.0 | force/length | NEEDS_JUSTIFICATION |
| `target_radius` ($R_\text{BM}^0$) | 6.0 | 20.0 | µm | Set slightly above organoid radius to allow initial settling |
| `ecm_degradation_rate` ($\dot{R}$) | 0.4 | 0.4 | µm/time-unit | NEEDS_JUSTIFICATION |
| `ecm_max_radius` ($R_\text{max}$) | 20.0 | 80.0 | µm | NEEDS_JUSTIFICATION; allows 3–4× expansion |

---

## 6. Cell Type / Heterogeneity Parameters

| Parameter | Value | Justified? |
|-----------|-------|------------|
| Core cell BM stiffness | 1.0 | NEEDS_JUSTIFICATION |
| Surface cell BM stiffness | 5.0 | NEEDS_JUSTIFICATION; 5× stiffer based on BM being at surface |
| Stem cell probability (factory) | 0.25 (2D), 0.25 (3D) | NEEDS_JUSTIFICATION |

---

## 7. Numerical / Infrastructure Parameters

| Parameter | Apical Constriction | 2D Threshold | 3D Threshold | Justified? |
|-----------|-------------------|-------------|-------------|------------|
| Movement threshold | — | 50.0 | 500.0 | Generous; prevents crash on large displacements |
| Birth time | $-10.0$ | Factory-set | Factory-set | Mature cells, avoids initial division |
| Cell cycle model | UniformG1 | OrganoidCellFactory | OrganoidCellFactory | — |

---

## Summary: Parameters Requiring Justification

1. **Apical constriction strength** (10.0): No direct experimental calibration;
   qualitative behaviour validated against Martin et al. (2009).
2. **BM stiffness** (2.0): Chosen for moderate confinement; literature values
   for Matrigel stiffness range 0.1–1.0 kPa (Gjorevski et al., 2016) but
   mapping to dimensionless Chaste units is non-trivial.
3. **ECM degradation rate** (0.4): Controls expansion speed; no direct
   experimental data used.
4. **Spring stiffnesses** (30.0, 35.0): Standard range for Chaste node-based
   models but not calibrated to specific tissue.
5. **Geometry** (organoid radii, cell counts): Chosen for computational
   tractability rather than quantitative biological values.
