# Physical Units System — Parameter Configuration Guide

## Overview

The TissueMorphology simulator now supports **physical units** (kPa, Pa·s, N/m, etc.) for parameter definition, automatically converting them to simulation units. This allows you to:
- Define parameters in standard biological units
- Compare simulations directly with experimental values
- Switch between physical and non-dimensional modes with a single flag

## Quick Start

### Enable Physical Units

Edit your JSON parameter file (e.g., `params-Node2d.json`):

```json
{
  "unitConversion": {
    "usePhysicalUnits": true,
    "cellDiameterMicrometers": 10.0,
    "tissueViscosityPaS": 1.0
  },
  ...
}
```

### Define Stiffness in kPa

With `usePhysicalUnits: true`, you can now express stiffness in **kPa** instead of non-dimensional units:

```json
"forces": {
  "ECMConfinementForce": {
    "cellGhostSpringStiffness": 10.0,    # Now in kPa (was non-dimensional before)
    ...
  },
  "GhostNodeECM": {
    "ghostGhostStiffness": 5.0,        # kPa
    "ghostCellGhostStiffness": 5.0,    # kPa
    "ghostE0": 2.0,      # kPa
    "ghostE1": 1.0      # kPa
  }
}
```

## Configuration Fields

### `unitConversion` Section (JSON)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `usePhysicalUnits` | bool | `false` | Enable/disable physical unit conversion |
| `cellDiameterMicrometers` | double | `10.0` | Cell diameter (typical: 8–15 µm) |
| `tissueViscosityPaS` | double | `1.0` | Tissue viscosity (typical: 0.1–10 Pa·s) |

**Backward Compatibility:** When `usePhysicalUnits = false` (default), all parameters work as before. No code changes required to existing simulations.

---

## Unit Conversion System

### Fundamental Assumptions

The conversion system models cells as **spheres in a viscous medium** with Stokes drag:

$$F_{\text{drag}} = 6\pi \eta r v$$

where:
- η = tissue viscosity (Pa·s)
- r = cell radius (m)
- v = velocity (m/s)

### Conversion Factors (Derived Automatically)

Given:
- Cell diameter: **d_CD** (in micrometers, e.g., 10 µm)
- Cell radius: **r** = d_CD / 2 = 5 µm = 5 × 10⁻⁶ m
- Tissue viscosity: **η** = 1 Pa·s (default)

Computed quantities:

1. **Drag coefficient** (Pa·s·m):
   $$\eta_{\text{drag}} = 6\pi \eta r \approx 9.42 \times 10^{-5} \text{ Pa·s·m}$$

2. **Stiffness conversion factor** (Pa·m⁻¹):
   $$k_{\text{factor}} = \frac{\eta_{\text{drag}}}{d_{\text{CD (m)}}} = \frac{9.42 \times 10^{-5}}{10 \times 10^{-6}} \approx 9.42$$

   **Translation:** 1 kPa stiffness → 9.42 simulation units

3. **Pressure conversion factor** (Pa):
   $$p_{\text{factor}} = \frac{\eta_{\text{drag}}}{A_{\text{cell}}} \text{ where } A_{\text{cell}} = (10 \text{ µm})^2 = 10^{-10} \text{ m}^2$$

   **Translation:** 1 Pa pressure → 9.42 × 10⁻⁴ simulation units (approximately)

### Conversion Examples (10 µm cells, 1 Pa·s viscosity)

| Physical Value | Simulation Units | Description |
|----------------|------------------|-------------|
| 5 kPa | 47.1 | ECM stiffness (soft tissue) |
| 20 kPa | 188 | ECM stiffness (moderately stiff) |
| 100 kPa | 942 | ECM stiffness (stiff/fibrotic) |
| 50 N/m | 5.3 | Cell-cell adhesion (per interface) |
| 100 N/m | 10.6 | Cell-cell adhesion (stiff) |
| 10 Pa | 10.6 | Lumen pressure (small) |
| 100 Pa | 106 | Lumen pressure (moderate) |

---

## Which Parameters Support Physical Units?

### ✅ Stiffness Parameters (in kPa when `usePhysicalUnits=true`)

| Parameter | Location | Use Case |
|-----------|----------|----------|
| `cellGhostSpringStiffness` | `forces.ECMConfinementForce` | ECM elastic modulus |
| `springStiffness` | `forces.SpringForce` | Cell-cell adhesion spring |
| `bendingStiffness` | `forces.RingSmoothingForce` | Monolayer curvature bending |
| `ghostGhostStiffness` | `forces.GhostNodeECM` | ECM-ECM fiber interactions |
| `ghostCellGhostStiffness` | `forces.GhostNodeECM` | Cell-ECM interactions |
| `ghostE0` | `forces.GhostNodeECM.ViscoelasticECM` | Long-time elastic modulus (Maxwell model) |
| `ghostE1` | `forces.GhostNodeECM.ViscoelasticECM` | Transient elastic modulus (Maxwell model) |
| `nhMembraneSurface` | `forces.VertexModel` | Membrane elasticity (vertex models) |
| `gammaApical` | `forces.VertexModel` | Apical surface tension (vertex models) |
| `gammaBasal` | `forces.VertexModel` | Basal surface tension (vertex models) |
| `gammaLateral` | `forces.VertexModel` | Lateral surface tension (vertex models) |

### ✅ Pressure/Force Parameters (in Pa when `usePhysicalUnits=true`)

| Parameter | Location | Use Case |
|-----------|----------|----------|
| `lumenPressure` | `forces.LumenPressureForce` | Hydrostatic lumen pressure |
| `apicalConstrictionStrength` | `forces.ApicalConstrictionForce` | Apical surface contractility |
| `polarityBendingStrength` | `forces.PolarlityForce` | Cell polarity-driven bending |

### ⚠️ Parameters NOT Converted (Keep as Non-Dimensional)

These remain **dimensionless** even with `usePhysicalUnits=true`:

| Parameter | Reason |
|-----------|--------|
| `ghostDamping` | Damping coefficient (usually fixed at ~5) |
| All fractions (stemFraction, transitFraction, etc.) | Already dimensionless |
| Geometry (cell counts, radii in CD) | Already normalized to cell diameters |
| Rates (ecmDegradationRate, fibreRemodelingRate) | Per-unit-time rates (keep as h⁻¹) |
| Probabilities (probStemToStem, etc.) | Already normalized [0,1] |

---

## Practical Workflow

### Example: Modeling Intestinal Crypts (Healthy Tissue)

**Step 1: Set Unit Conversion Configuration**
```json
{
  "unitConversion": {
    "usePhysicalUnits": true,
    "cellDiameterMicrometers": 10.0,
    "tissueViscosityPaS": 1.0
  },
  ...
}
```

**Step 2: Consult Literature for Target Stiffness**
- Normal intestinal mucosa ECM: **5–10 kPa** (from experimental measurements)

**Step 3: Define Stiffness in Physical Units**
```json
"forces": {
  "ECMConfinementForce": {
    "cellGhostSpringStiffness": 7.5,         # 7.5 kPa — normal tissue
    ...
  }
}
```

**Step 4: Run Simulation**
```bash
./CryptBuddingApp -config params-Node2d.json
```

**Output (console):**
```
Loaded 42 parameters from: params-Node2d.json
Converting parameters from physical to simulation units...
Cell radius: 5.00e-06 m
Drag coefficient: 9.42e-05 Pa·s·m
Stiffness conversion factor: 9.42
Pressure conversion factor: 9.42e-04

Converted cellGhostSpringStiffness: 7.5 kPa → 70.7 (simulation units)
```

---

## Changing Cell Type or Tissue Properties

### Example: Simulate Larger Cells (12 µm diameter)

```json
{
  "unitConversion": {
    "usePhysicalUnits": true,
    "cellDiameterMicrometers": 12.0,    # Adjust cell diameter
    "tissueViscosityPaS": 1.0
  },
  ...
}
```

The conversion factors are **automatically recalculated**. All stiffness values remain unchanged in your config, but map to different simulation units.

### Example: More Viscous Tissue (2 Pa·s, e.g., mucus-rich environment)

```json
{
  "unitConversion": {
    "usePhysicalUnits": true,
    "cellDiameterMicrometers": 10.0,
    "tissueViscosityPaS": 2.0           # Increase viscosity
  },
  ...
}
```

Drag coefficient doubles → stiffness conversion factor doubles → all stiffness parameters become smaller in simulation units (more compliant relative to viscous damping).

---

## Experimental Data Comparison

### Intestinal ECM Stiffness (Literature)

| Tissue State | Stiffness (kPa) | Reference |
|--------------|-----------------|-----------|
| Normal mucosa | 5–10 | Stomach/colon biopsies |
| Inflamed (IBD) | 15–25 | Inflammatory response |
| Fibrotic/Scarred | 50–100+ | Post-inflammatory remodeling |
| Tumor-adjacent stroma | 30–50 | Desmoplastic reaction |

**Use these values to parameterize your simulations:**
```json
"forces": {
  "ECMConfinementForce": {
    "cellGhostSpringStiffness": 7.5       # Normal
    // "cellGhostSpringStiffness": 20.0   # Inflamed
    // "cellGhostSpringStiffness": 75.0   # Fibrotic
  }
}
```

### Cell-Cell Adhesion Stiffness (AFM Measurements)

| Junction Type | Stiffness (N/m) | Notes |
|---------------|-----------------|-------|
| Weak/immature | 10–50 | Early in development |
| Moderate | 50–200 | Typical epithelium |
| Strong/mature | 200–500 | Confluent monolayers |

Convert to simulation units: **N/m value × 0.106** (approximately, for 10 µm cells, 1 Pa·s)

---

## Troubleshooting

### Issue: Conversion factors seem wrong

**Check:** Did you set `usePhysicalUnits: true`?

The code only applies conversion if the flag is enabled. Default is `false` for backward compatibility.

### Issue: Simulation behaves differently after enabling physical units

**Expected behavior.** If your stiffness values change from (e.g.) 50 to 5.0 kPa, the conversion factor (~9.42) means:
- Old value: 50 (sim units)
- New value: 5.0 × 9.42 ≈ 47.1 (sim units)

The simulation should be similar, not identical (unless you account for all interactions).

### Issue: Parameters seem unreasonably large/small

**Check:** Are your cell diameter or viscosity assumptions realistic?

Example: If you set `cellDiameterMicrometers: 1000` (1 mm), the conversion factor becomes 9420 — all stiffness values will explode. Use typical values: 8–15 µm.

---

## Disabling Physical Units

To revert to non-dimensional mode:

```json
{
  "unitConversion": {
    "usePhysicalUnits": false,    # Disable conversion
    ...
  },
  "forces": {
    "ECMConfinementForce": {
      "cellGhostSpringStiffness": 70.7,       # Back to simulation units
      ...
    }
  }
}
```

All parameters are interpreted as **non-dimensional simulation units** when `usePhysicalUnits: false`.

---

## References

1. **Unit Conversion Guide:** See [UNIT_CONVERSION_GUIDE.md](./UNIT_CONVERSION_GUIDE.md) for detailed mathematical derivations
2. **Parameter Reference:** See [apps/PARAMS_README.md](./apps/PARAMS_README.md) for all parameters and their meanings
3. **ECM Force Models:** See [ECM_Forces_Comparison.md](./ECM_Forces_Comparison.md) for force law details

---

## FAQ

**Q: Can I use physical units for geometry (cell counts, radii)?**
A: No. Geometry is always in cell diameters (CD). Physical unit conversion applies only to stiffness/pressure parameters.

**Q: What if I don't specify `usePhysicalUnits`?**
A: Defaults to `false` (non-dimensional mode). No conversion is applied.

**Q: Can I mix physical and non-dimensional parameters in the same config?**
A: No. The flag is global. When enabled, all supported stiffness/pressure parameters are converted.

**Q: What's the best cell diameter to use?**
A: For intestinal epithelium: **10 µm** (default). For other tissues, measure or use literature values: pancreatic cells ~15 µm, hepatocytes ~20 µm, fibroblasts ~10 µm.

**Q: Why does my simulation run slower after enabling physical units?**
A: It doesn't. The conversion happens once at startup. The only difference is different stiffness values in simulation units.

---

## Advanced: Custom Conversion Factors

If you want to manually override the computed conversion factors (advanced use), you would need to modify the C++ code. The conversion logic is in `CryptBuddingParams::ComputeConversionFactors()` — see [CryptBuddingParams.hpp](./apps/src/CryptBuddingParams.hpp) for details.

For most users, the automatic Stokes-drag-based factors should suffice.

