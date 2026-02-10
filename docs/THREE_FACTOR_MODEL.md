# Three-Factor Mechanical Model for Organoid Morphogenesis

## Overview
Implementation of the mechanical model from Metzcar et al. 2025, which describes organoid morphogenesis through three key mechanical factors that drive tissue folding and invagination.

## Mathematical Framework

### Core Mechanical Factors

The model describes tissue deformation through three primary mechanical drivers:

1. **Apical Constriction (AC)**
   - Cells reduce their apical surface area
   - Creates local curvature by wedging cells
   - Force: $F_{AC} = k_{AC} \cdot (A_{apical} - A_{target})$

2. **Differential Adhesion (DA)**
   - Varying adhesion between cell layers
   - Drives cell sorting and layer separation
   - Force: $F_{DA} = k_{DA} \cdot \Delta E_{adhesion}$

3. **Basement Membrane (BM) / ECM Stiffness**
   - Constrains tissue expansion
   - Resists deformation
   - Force: $F_{BM} = k_{BM} \cdot (r - r_{0}) \cdot \hat{r}$

### Combined Tissue Dynamics

The total force on each cell includes:
$$\vec{F}_{total} = \vec{F}_{spring} + \vec{F}_{AC} + \vec{F}_{DA} + \vec{F}_{BM}$$

Cell movement follows overdamped dynamics:
$$\eta \frac{d\vec{r}_i}{dt} = \vec{F}_{total}(i)$$

## Implementation Plan

### Phase 1: Core Force Classes ✅
- [x] `BasementMembraneForce` - ECM constraint implemented with degradation
- [x] `ApicalConstrictionForce` - Implement apical surface area reduction
- [x] `DifferentialAdhesionForce` - Cell-type specific adhesion energies

### Phase 2: Cell Properties ✅
- [x] Cell data fields:
  - `is_apical` - Mark apical cells for constriction
  - `cell_type` - Identify apical/basal layers (0=basal, 1=apical)
  - Cell area estimation via neighbor distances
  
### Phase 3: 2D Validation Tests (In Progress)
- [x] `Test2dApicalConstriction` - Single factor: apical wedging
- [ ] `Test2dDifferentialAdhesion` - Single factor: cell sorting
- [ ] `Test2dThreeFactorModel` - All three factors combined

### Phase 4: 3D Organoid Tests (Planned)
- [ ] `Test3dThreeFactorModel` - Full 3D morphogenesis
- [ ] `TestParameterSweep` - Systematic variation of $k_{AC}$, $k_{DA}$, $k_{BM}$

## Key Parameters

| Parameter | Symbol | Typical Range | Units |
|-----------|--------|---------------|-------|
| Apical constriction strength | $k_{AC}$ | 1-10 | force/area |
| Adhesion differential | $k_{DA}$ | 0.1-5 | energy |
| BM stiffness | $k_{BM}$ | 1-20 | force/length |
| Target apical reduction | $\alpha$ | 0.3-0.7 | dimensionless |

## Minimal Working Example

```cpp
// Pseudocode for combined model
MAKE_PTR(ApicalConstrictionForce<2>, p_ac_force);
p_ac_force->SetConstrictionStrength(5.0);
p_ac_force->SetTargetReduction(0.5);  // 50% area reduction

MAKE_PTR(DifferentialAdhesionForce<2>, p_da_force);
p_da_force->SetApicalAdhesion(2.0);
p_da_force->SetBasalAdhesion(1.0);

MAKE_PTR(BasementMembraneForce<2>, p_bm_force);
p_bm_force->SetBasementMembraneParameter(3.0);
p_bm_force->SetTargetRadius(20.0);

simulator.AddForce(p_ac_force);
simulator.AddForce(p_da_force);
simulator.AddForce(p_bm_force);
```

## Expected Outcomes

### Qualitative Behaviors
1. **AC dominant**: Localized invagination pits
2. **DA dominant**: Cell layer segregation
3. **BM dominant**: Constrained spherical morphology
4. **Combined**: Realistic organoid folding patterns

### Quantitative Metrics
- Tissue curvature: $\kappa = 1/R_{curvature}$
- Invagination depth: $d_{inv}$
- Layer separation distance
- Total elastic energy

## Next Steps

1. Implement `ApicalConstrictionForce.hpp` with cell area tracking
2. Implement `DifferentialAdhesionForce.hpp` with cell-type interactions
3. Create `Test2dApicalConstriction.hpp` for validation
4. Extend to 3D with `Test3dThreeFactorModel.hpp`

## References

Metzcar et al. (2025). "Agent-based modeling framework for organoid morphogenesis"
- Supplementary Material: Three-factor mechanical model
- Focus: Minimal mechanistic model for tissue folding
