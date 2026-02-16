# Simulation Parameters: 3D Organoid Formation

Parameters from `Test3dOrganoidFormation`, `Test3dCryptOrganoid`, and
`Test3dVertexCryptOrganoid`.

---

## 1. Geometry Parameters

### Test3dOrganoidFormation

| Parameter | Basic | Stiffness | LongTerm | Units | Justified? |
|-----------|-------|-----------|----------|-------|------------|
| `num_cells` | 30 | 30 | 40 | — | NEEDS_JUSTIFICATION; kept small for stability |
| `organoid_radius` | 20.0 | 18.0 | 18.0 | µm | NEEDS_JUSTIFICATION |
| Interaction cutoff | 5.0 | 5.0 | 5.0 | µm | NEEDS_JUSTIFICATION |
| Stem probability (factory) | 0.15 | 0.15 | 0.25 | — | NEEDS_JUSTIFICATION |

### Test3dCryptOrganoid

| Parameter | Value | Units | Justified? |
|-----------|-------|-------|------------|
| `num_initial_cells` | 100 | — | Reasonable for small organoid; NEEDS_JUSTIFICATION |
| `organoid_radius` | 25.0 | µm | Typical organoid starting radius (Sato et al., 2009) |
| `shell_thickness` | 3.0 | µm | Approximates epithelial monolayer thickness; partially justified |
| `interaction_cutoff` | 15.0 | µm | ~3 cell diameters; NEEDS_JUSTIFICATION |
| `ecm_domain_half` | 80.0 | µm | ECM domain $[-80, 80]^3$; NEEDS_JUSTIFICATION |
| `ecm_grid_spacing` | 10.0 | µm per voxel | NEEDS_JUSTIFICATION; no convergence study |

### Test3dVertexCryptOrganoid

| Parameter | Relaxation | Growth | Units | Justified? |
|-----------|-----------|--------|-------|------------|
| `num_cells` | 100 | 200 | — | Realistic for small organoid |
| `sphere_radius` (inner) | 10.0 | 10.0 | µm (dimensionless) | NEEDS_JUSTIFICATION |
| Cell height $h$ | Computed from $\gamma$ | Same | µm | Derived from surface tensions (Drozdowski & Schwarz, 2025) |
| T1 threshold $l_{T1}$ | Computed from $\gamma_L$ | Same | µm | Derived; Yes |
| `ecm_domain_half` | — | 40.0 | µm | NEEDS_JUSTIFICATION |
| `ecm_grid_spacing` | — | 5.0 | µm | Finer than node-based; NEEDS_JUSTIFICATION |

---

## 2. Time Parameters

| Parameter | OrganoidFormation (Basic) | OrganoidFormation (LongTerm) | CryptOrganoid | VertexCrypt (Relax) | VertexCrypt (Grow) | Units |
|-----------|-------------------------|----------------------------|---------------|--------------------|--------------------|-------|
| `dt` | 0.002 | 0.002 | 0.01 | 0.001 | 0.006 | hours (Chaste) |
| `end_time` | 2.0 | 100.0 | 168.0 (7 days) | 5.0 | 100.0 | hours |
| `sampling_multiple` | 50 | 50 | 20 | 20–50 | 100 | timesteps |
| Output interval | 0.1 | 0.1 | 0.2 | 0.02–0.05 | 0.6 | hours |

**Justification:** dt values empirically chosen for numerical stability.
Vertex model requires smaller dt during relaxation due to surface tension
stiffness. All NEED formal CFL analysis.

---

## 3. Cell–Cell Spring Parameters (Node-Based Only)

| Parameter | Basic | Stiffness | LongTerm | CryptOrganoid | Units | Justified? |
|-----------|-------|-----------|----------|---------------|-------|------------|
| `spring_stiffness` ($k_s$) | 30.0 | 25.0 | 35.0 | 20.0 | force/length | NEEDS_JUSTIFICATION |
| `division_resting_spring_length` | 0.5 | — | — | 0.5 | cell diameters | Standard Chaste default |
| `spring_growth_duration` | 1.0 | — | — | 1.0 | hours | Standard Chaste default |
| `cutoff_length` | default | default | default | 15.0 | µm | NEEDS_JUSTIFICATION |

---

## 4. Surface Tension Parameters (Vertex Model Only)

| Parameter | Value | Units | Source | Justified? |
|-----------|-------|-------|--------|------------|
| $\gamma_A$ (apical) | 0.85 | dimensionless | Drozdowski & Schwarz (2025) | Yes |
| $\gamma_B$ (basal) | 0.85 | dimensionless | Drozdowski & Schwarz (2025) | Yes |
| $\gamma_L$ (lateral) | 0.7 | dimensionless | Drozdowski & Schwarz (2025) | Yes |
| Protorosette formation prob. | 1.0 | — | OrganoidChaste default | Yes |
| Protorosette resolution prob. | 0.1 | per timestep | OrganoidChaste default | Yes |

### Simulated Annealing (T1 Transitions)

| Phase | Rate | Max steps | Amplitude | Active T1? |
|-------|------|-----------|-----------|-----------|
| Relaxation | 0.0 | 50–190 | 0–10 | No |
| Growth | 0.003 | 1,900,000 | 1.0 | Yes |

**Note:** Simulated annealing parameters control the stochastic T1 transition
search. Growth-phase parameters enable active neighbour rearrangements needed
for buckling.

---

## 5. Basement Membrane Parameters

| Parameter | OrganoidFormation (Basic) | LongTerm | CryptOrganoid | VertexCrypt (Relax) | VertexCrypt (Grow) | Units | Justified? |
|-----------|-------------------------|----------|---------------|--------------------|--------------------|-------|------------|
| $k_\text{BM}$ | 1.5 | 2.0 | 5.0 | 2.0 | 3.0 | force/length | NEEDS_JUSTIFICATION |
| $R_\text{BM}^0$ | 25.0 | 20.0 | 30.0 | $R + h + 1$ | $R + h + 1$ | µm | Set slightly above organoid radius |
| $\dot{R}$ (degradation rate) | 0.3 | 1.0 | 0.15 | 0.0 | 0.002 | µm/hour | NEEDS_JUSTIFICATION |
| $R_\text{max}$ | 60.0 | 70.0 | 80.0 | — | 25.0 | µm | NEEDS_JUSTIFICATION |

**Design notes:**
- Relaxation phase (vertex): no degradation ($\dot{R} = 0$) to maintain confinement
  while mesh equilibrates.
- `CryptOrganoid` uses slower degradation (0.15 µm/h) for a 7-day simulation,
  allowing gradual expansion to ~80 µm.
- `OrganoidFormation` LongTerm uses faster degradation (1.0 µm/h) for a shorter
  run.

---

## 6. 3D ECM Field Parameters

| Parameter | CryptOrganoid | VertexCrypt | Units | Justified? |
|-----------|---------------|-------------|-------|------------|
| Fiber pattern | radial | radial | — | Mimics Matrigel; partially justified (Gjorevski et al., 2016) |
| Grid spacing | 10.0 | 5.0 | µm | NEEDS_JUSTIFICATION |
| Domain extent | $[-80, 80]^3$ | $[-40, 40]^3$ | µm | Must contain organoid + expansion |
| $k_d$ (degradation rate) | 0.002 | 0.002 | per timestep | NEEDS_JUSTIFICATION |
| $k_r$ (remodeling rate) | 0.05 | 0.05 | per timestep | NEEDS_JUSTIFICATION |
| $k_p$ (deposition rate) | 0.0003 | 0.0003 | per timestep | NEEDS_JUSTIFICATION |
| Base speed ($v_0$) | 0.3 | 0.1 | µm/hour | NEEDS_JUSTIFICATION; much slower than invasion tests |
| ECM sensitivity ($s$) | 1.0 | 1.0 | dimensionless | Full sensitivity; Metzcar et al. (2025) |
| Degradation enabled | Yes | Yes | — | — |
| Remodeling enabled | Yes | Yes | — | — |
| Deposition enabled | No | No | — | Simplified; future addition |

**Design note:** Base speed is much lower in organoid formation (0.1–0.3 µm/h)
than in invasion tests (1.25 µm/min) because organoid cell migration is slow
and constrained by the epithelial monolayer, unlike freely migrating mesenchymal
cells.

---

## 7. Volume Homeostasis Parameters

### Node-based (`VolumeTrackingModifier` + `ContactInhibitionCellCycleModel`)

| Parameter | Value | Units | Justified? |
|-----------|-------|-------|------------|
| Quiescent volume fraction ($\phi_q$) | 0.7 | dimensionless | NEEDS_JUSTIFICATION |
| Equilibrium volume ($V_\text{eq}$) | 1.0 | cell volumes | Standard normalisation |

### Vertex-based (`GeometricalTargetVolumeModifier`)

| Parameter | Value | Units | Justified? |
|-----------|-------|-------|------------|
| Reference target volume | $\frac{4\pi}{3N}(R_\text{out}^3 - R_\text{in}^3)$ | µm³ | Derived from geometry; Yes |
| Growth duration | 0.0 | hours | Immediate volume target assignment |
| T1 adaptation duration | 0.1 | hours | OrganoidChaste default |

---

## 8. Cell Cycle Parameters

### CryptOrganoid (ContactInhibitionCellCycleModel)

| Parameter | Value | Units | Justified? |
|-----------|-------|-------|------------|
| Birth time desync | $-\mathcal{U}(0, 18)$ | hours | Prevents synchronsied division |

### VertexCryptOrganoid (FixedG1GenerationalCellCycleModel)

| Parameter | Stem | Transit | Differentiated | Units | Justified? |
|-----------|------|---------|----------------|-------|------------|
| G1 duration | 25.0 | 15.0 | $10^6$ | hours | NEEDS_JUSTIFICATION |
| S duration | 25.0 | 25.0 | — | hours | NEEDS_JUSTIFICATION; unusually long |
| G2 duration | 0.001 | 0.001 | — | hours | Negligible (vertex model handles shape change) |
| M duration | 0.001 | 0.001 | — | hours | Negligible |
| Max transit generations | 2 | 2 | 0 | — | NEEDS_JUSTIFICATION |
| Birth time desync | $-\mathcal{U}(0, 10)$ | — | — | hours | Prevents synchronised division |

**Note:** S-phase duration of 25 h in the vertex model is unusually long
(typical: 6–8 h). This appears to be a deliberate choice to slow overall
proliferation rate, enabling the vertex mesh to handle growth without
instability. NEEDS_JUSTIFICATION.

---

## 9. Cell Type Assignment

| Parameter | CryptOrganoid | VertexCryptOrganoid | Justified? |
|-----------|---------------|---------------------|------------|
| Crypt base threshold ($z_\text{frac}$) | $< -0.5$ | $< -0.5$ | NEEDS_JUSTIFICATION; approximates crypt anatomy |
| TA zone threshold | $[-0.5, 0.3)$ | $[-0.5, 0.3)$ | NEEDS_JUSTIFICATION |
| Differentiated threshold | $\geq 0.3$ | $\geq 0.3$ | NEEDS_JUSTIFICATION |
| Paneth cell frequency | Every 3rd stem cell | — | NEEDS_JUSTIFICATION; Paneth ~20–30% of crypt base |
| ECM degradation factor (stem) | 2.0 | — | NEEDS_JUSTIFICATION; stem cells degrade BM faster |
| ECM degradation factor (TA) | 1.0 | — | Baseline |
| ECM degradation factor (diff) | 0.5 | — | Lower; differentiated cells are less invasive |

---

## 10. Numerical / Infrastructure Parameters

| Parameter | OrganoidFormation | CryptOrganoid | VertexCryptOrganoid | Justified? |
|-----------|------------------|---------------|---------------------|------------|
| Movement threshold | 50.0 | 50.0 | `RestrictVertexMovement = false` | Generous; prevents crashes |
| Volume relaxation | — | — | `DoInitialVolumeRelaxation = true` | Yes; standard for vertex |
| Progress update interval | — | 100 timesteps | — | Monitoring |

---

## Summary: Parameters Requiring Justification

1. **Spring stiffness** (20–35): varies across tests without clear rationale;
   needs calibration or sensitivity analysis.
2. **BM stiffness** (1.5–5.0): range explored but not systematically justified.
3. **BM degradation rate** (0.002–1.0): varies by 500× across tests; the slow
   rate in CryptOrganoid vs fast rate in OrganoidFormation needs explanation.
4. **ECM dynamics rates** ($k_d$, $k_r$, $k_p$): inherited defaults from
   `DynamicECMField3d`; no experimental calibration.
5. **ECM base speed** (0.1–0.3 µm/h): much slower than invasion tests;
   justified qualitatively but not quantitatively.
6. **Cell cycle S-phase** (25 h in vertex): non-standard; appears to be a
   numerical convenience rather than biological value.
7. **Cell type thresholds** ($z_\text{frac}$): approximate anatomical zones with
   no quantitative basis.
8. **Cell counts and radii**: chosen for computational tractability; sensitivity
   analysis needed.
