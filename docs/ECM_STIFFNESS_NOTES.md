# ECM Stiffness in Dynamic ECM Model

## Current Implementation

Currently, the `DynamicECMField` stores:
- **ECM density** (0-1): Amount of ECM present
- **ECM anisotropy** (0-1): How aligned the fibers are
- **Fiber angle** (radians): Orientation of aligned fibers

## Missing: ECM Stiffness

We are **NOT** currently storing ECM stiffness/rigidity. This could be important for:

### Biological Relevance
1. **Tissue mechanics**: Stiffness affects cell traction forces
2. **Durotaxis**: Cells migrate toward stiffer substrates
3. **Mechanosensing**: Cells respond to ECM stiffness via integrins
4. **Disease models**: Tumor stroma is often stiffer than normal tissue

### Implementation Options

If you want to add ECM stiffness:

```cpp
struct ECMGridCell
{
    double fiber_angle;    // Orientation (radians)
    double density;        // Amount of ECM (0-1)
    double anisotropy;     // Alignment (0-1)
    double stiffness;      // NEW: Rigidity (kPa or relative units)
};
```

### Stiffness Effects on Cell Migration

1. **Speed modulation**: Cells migrate faster on intermediate stiffness
   ```cpp
   double optimal_stiffness = 10.0;  // kPa
   double stiffness_factor = exp(-pow(stiffness - optimal_stiffness, 2) / (2*5*5));
   effective_speed = base_speed * stiffness_factor;
   ```

2. **Durotaxis**: Cells bias migration toward stiffer ECM
   ```cpp
   c_vector<double, 2> stiffness_gradient = ComputeStiffnessGradient(position);
   migration_direction += durotaxis_strength * stiffness_gradient;
   ```

3. **Stiffness-density coupling**: Dense ECM is typically stiffer
   ```cpp
   stiffness = base_stiffness * (1.0 + 2.0 * density);
   ```

### Typical ECM Stiffness Values (Pa)
- Soft tissue: 100-1000 Pa (0.1-1 kPa)
- Normal mammary tissue: ~200 Pa
- Breast tumor: ~4000 Pa (4 kPa)
- Collagen gel: 100-10,000 Pa
- Bone: ~10 GPa

### References
- Metzcar et al. 2025: Focus on contact guidance, not stiffness
- Painter 2009: Contact guidance model (no stiffness)
- Pathak & Kumar 2012: Durotaxis models
- Lo et al. 2000: Cell response to stiffness gradients

## Recommendation

For the **invasive front model** (Metzcar 2025, Section 3.1), ECM stiffness is **not required** because:
1. The paper focuses on contact guidance (fiber orientation)
2. All cells are the same type (no differential durotaxis)
3. Simplicity: fewer parameters to tune

Add stiffness only if you need to model:
- Durotaxis (stiffness-directed migration)
- Traction force mechanics
- Tissue stiffening in disease progression
