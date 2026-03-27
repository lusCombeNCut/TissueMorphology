# Unit Conversion System Implementation Summary

## Overview

Added comprehensive unit conversion support to TissueMorphology simulation parameters, allowing researchers to define parameters in physical units (kPa, Pa·s, etc.) instead of non-dimensional simulation units. The system automatically converts these values based on cell diameter and tissue viscosity, with full backward compatibility preserved.

## Files Modified

### 1. `/home/orlando/Thesis/Chaste/projects/TissueMorphology/apps/src/CryptBuddingParams.hpp`

**Changes:**

#### Added Struct Fields (lines 39-45)
Three new fields for controlling unit conversion:
- `bool usePhysicalUnits` - Master toggle (default: false)
- `double cellDiameterMicrometers` - Cell diameter (default: 10.0 µm)
- `double tissueViscosityPaS` - Tissue viscosity (default: 1.0 Pa·s)

#### New Method: `ComputeConversionFactors()` (lines 201-270)
Comprehensive conversion method that:
- Calculates cell radius in meters
- Computes drag coefficient: `η_drag = 6π·η·r` (Pa·s·m)
- Derives stiffness conversion factor: `k_factor = η_drag / d_cell`
- Derives pressure conversion factor: `p_factor = η_drag`
- Applies conversions to all relevant parameters:
  - **Stiffness parameters** (multiplied by `k_factor`):
    - `ecmStiffness`, `springStiffness`, `bendingStiffness`
    - Ghost node ECM: `ghostGhostStiffness`, `ghostCellGhostStiffness`, `ghostRelaxedStiffness`, `ghostRelaxationModulus`
    - Nagai-Honda surface tensions: `gammaApical`, `gammaBasal`, `gammaLateral`
    - Adhesion parameters: `nhMembraneSurface`, `nhCellCellAdhesion`, `nhBoundaryAdhesion`, and all per-type adhesions (12 total)
  - **Pressure/force parameters** (multiplied by `p_factor`):
    - `lumenPressure`, `apicalConstrictionStrength`
    - `polarityBendingStrength`
- Outputs conversion factors to console for debugging

#### Updated `SetDefaults()` (lines 272-277)
Initializes new conversion fields:
```cpp
usePhysicalUnits = false;
cellDiameterMicrometers = 10.0;         // Typical mammalian cell
tissueViscosityPaS = 1.0;               // Typical tissue
```

#### Updated `ApplyConfigMap()` (lines 586-588)
- Added parameter loading for unit conversion fields (loaded before stiffness/pressure parameters)
- Added call to `ComputeConversionFactors()` at end of method to apply conversions after all parameters are loaded

#### Updated `SaveToFile()` (INI format)
- Added `[UnitConversion]` section with documentation of:
  - Unit conversion fields and typical ranges
  - Guidance on cell diameter (8-15 µm)
  - Guidance on tissue viscosity (0.1-10 Pa·s)
  - Guidance on ECM stiffness (0.1-10 kPa)
  - Guidance on lumen pressure (0.1-2.0 kPa)

#### Updated `SaveToJson()` (JSON format)
- Added `"unitConversion"` object with all three conversion parameters

### 2. `/home/orlando/Thesis/Chaste/projects/TissueMorphology/apps/config/default_params.ini`

**Changes:**

Added comprehensive `[UnitConversion]` section with:
- Documentation of the unit conversion system
- Three configurable parameters with explanations
- Typical ranges from biological literature for:
  - Cell diameter: 8-15 µm (mammalian cells, 10 µm typical)
  - Tissue viscosity: 0.1-10 Pa·s (water ≈ 0.001, cytoplasm ≈ 1-10, honey ≈ 2-10)
  - ECM stiffness: 0.1-10 kPa (biological range; collagen ≈ 1-5 kPa)
  - Lumen pressure: 0.1-2.0 kPa (osmotic/hydrostatic)

### 3. `/home/orlando/Thesis/Chaste/projects/TissueMorphology/apps/README.md`

**Changes:**

Added comprehensive "Unit Conversion System" section (before "Feature Toggles") with:
- Configuration table explaining all three parameters
- Mathematical formula for conversion factors
- List of all parameters affected by stiffness conversion
- List of all parameters affected by pressure conversion
- Three detailed examples:
  1. **Non-dimensional (backward compatible)**: Shows default behavior unchanged
  2. **Physical units with standard tissue**: Demonstrates 1 kPa ECM conversion
  3. **Softer tissue with larger cells**: Shows how cell diameter and viscosity scale effects

### 4. Created Four JSON Parameter Files

All files created in `/home/orlando/Thesis/Chaste/projects/TissueMorphology/apps/` with identical structure:

#### `/home/orlando/Thesis/Chaste/projects/TissueMorphology/apps/params-Node2d.json`
- 2D node-based model
- dt=0.005, endTime=168.0

#### `/home/orlando/Thesis/Chaste/projects/TissueMorphology/apps/params-Vertex2d.json`
- 2D vertex-based model
- dt=0.005, endTime=168.0
- `enableApicalConstriction=true`, `enableCurvatureBending=false`

#### `/home/orlando/Thesis/Chaste/projects/TissueMorphology/apps/params-Node3d.json`
- 3D node-based model
- dt=0.01, endTime=168.0

#### `/home/orlando/Thesis/Chaste/projects/TissueMorphology/apps/params-Vertex3d.json`
- 3D vertex-based model (OrganoidChaste)
- dt=0.01, endTime=168.0
- `enableApicalConstriction=true`, `enableCurvatureBending=false`

**Each file contains:**
- `unitConversion` section with default values (all disabled for backward compatibility)
- All other simulation parameters in hierarchical JSON structure
- Complete documentation structure for all force, ECM, and cell cycle parameters

## Conversion Mathematics

When `usePhysicalUnits=true`, the conversion proceeds as follows:

### Step 1: Compute drag coefficient
```
r_cell = cellDiameterMicrometers × 1e-6 / 2        (meters)
η_drag = 6π × tissueViscosityPaS × r_cell          (Pa·s·m)
```

For default parameters (10 µm, 1 Pa·s):
```
r_cell = 10×10^-6 / 2 = 5×10^-6 m
η_drag = 6π × 1.0 × 5×10^-6 ≈ 9.42×10^-5 Pa·s·m
```

### Step 2: Compute conversion factors
```
d_cell = cellDiameterMicrometers × 1e-6            (meters)
k_factor = η_drag / d_cell                          (dimensionless)
p_factor = η_drag                                   (Pa·s·m)
```

For default parameters:
```
d_cell = 10×10^-6 = 1×10^-5 m
k_factor = 9.42×10^-5 / 1×10^-5 ≈ 9.42
p_factor = 9.42×10^-5 Pa·s·m
```

### Step 3: Apply conversions
- Stiffness: `k_sim = k_physical[Pa] × k_factor`
- Pressure: `f_sim = p_physical[Pa] × p_factor`

**Example**: 1 kPa ECM stiffness becomes:
```
k_sim = 1000 Pa × 9.42 ≈ 9420 (simulation units)
```

## Backward Compatibility

The implementation maintains full backward compatibility:

1. **Default state**: `usePhysicalUnits=false` - no conversions applied
2. **Existing parameters**: All non-dimensional default values unchanged
3. **Configuration files**: New parameters are optional; existing configs work unmodified
4. **Code changes**: Non-invasive; conversion happens after all parameters are loaded
5. **No API changes**: Force classes and simulation code unchanged

## Usage

### To use physical units in a simulation:

#### In INI format:
```ini
[UnitConversion]
usePhysicalUnits = true
cellDiameterMicrometers = 10.0
tissueViscosityPaS = 1.0

[ECM]
ecmStiffness = 1000.0    # Now interpreted as 1 kPa

[Forces]
[LumenPressureForce]
lumenPressure = 500.0    # Now interpreted as 500 Pa
```

#### In JSON format:
```json
{
  "unitConversion": {
    "usePhysicalUnits": true,
    "cellDiameterMicrometers": 10.0,
    "tissueViscosityPaS": 1.0
  },
  "ECM": {
    "ecmStiffness": 1000.0
  },
  "forces": {
    "LumenPressureForce": {
      "lumenPressure": 500.0
    }
  }
}
```

### Command line usage:
```bash
./CryptBuddingApp -model node2d -config params-Node2d.json
```

## Parameters Affected by Unit Conversion

### Stiffness Parameters (×k_factor):
1. ecmStiffness
2. springStiffness
3. bendingStiffness
4. ghostGhostStiffness
5. ghostCellGhostStiffness
6. ghostRelaxedStiffness
7. ghostRelaxationModulus
8. ecmConfinementStiffness
9. gammaApical
10. gammaBasal
11. gammaLateral
12. nhMembraneSurface
13. nhCellCellAdhesion
14. nhBoundaryAdhesion
15. nhStemStemAdhesion
16. nhStemTransitAdhesion
17. nhStemDiffAdhesion
18. nhTransitTransitAdhesion
19. nhTransitDiffAdhesion
20. nhDiffDiffAdhesion
21. nhStemBoundaryAdhesion
22. nhTransitBoundaryAdhesion
23. nhDiffBoundaryAdhesion

### Pressure/Force Parameters (×p_factor):
1. lumenPressure
2. apicalConstrictionStrength
3. polarityBendingStrength

## Testing Recommendations

1. **Verify backward compatibility**: Run simulations with default INI file (usePhysicalUnits=false) and confirm identical results
2. **Test conversion factors**: Load physical unit config and verify console output matches hand-calculated values
3. **Cross-check units**: For a 1 kPa stiffness, verify simulation parameter matches k_factor × 1000
4. **Validate edge cases**: Test with different cell diameters (8, 10, 15 µm) and viscosities (0.1, 1, 10 Pa·s)

## Related Documentation

- `UNIT_CONVERSION_GUIDE.md` (if exists) - Detailed parameter-by-parameter guidance
- `apps/README.md` - Configuration reference with examples
- `apps/config/default_params.ini` - INI template with typical values

## Notes

- Unit conversion is applied **after** all parameters are loaded from config files
- Conversion factors are printed to console for transparency and debugging
- The system is fully transparent: no internal API changes needed
- All 23 stiffness parameters and 3 pressure parameters support conversion
- Typical biological ranges documented in both INI and markdown documentation
