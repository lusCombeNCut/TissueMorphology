# CryptBuddingApp Configuration Reference

This document maps configuration parameters to their corresponding C++ classes.

## Unit Conversion System

Parameters can be defined in physical units (kPa, Pa·s, etc.) and automatically converted to simulation units using the tissue viscosity and cell diameter.

### Configuration

Three parameters control unit conversion in the `unitConversion` section:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `usePhysicalUnits` | bool | false | Enable/disable physical unit conversion |
| `cellDiameterMicrometers` | double | 10.0 | Cell diameter in micrometers (typical 8-15 µm) |
| `tissueViscosityPaS` | double | 1.0 | Tissue viscosity in Pa·s (typical 0.1-10 Pa·s) |

### Conversion Formula

When `usePhysicalUnits=true`, parameters are converted from physical units to simulation units as follows:

```
cellRadiusMeters = cellDiameterMicrometers × 1e-6 / 2
dragCoeff = 6π × tissueViscosityPaS × cellRadiusMeters  (Pa·s·m)
stiffnessConversionFactor = dragCoeff / cellDiameterMeters
pressureConversionFactor = dragCoeff  (Pa·s·m)
```

All stiffness parameters (Pa) are multiplied by `stiffnessConversionFactor`:
- `ecmStiffness`, `springStiffness`, `bendingStiffness`
- `ghostGhostStiffness`, `ghostCellGhostStiffness`, `ghostRelaxedStiffness`, `ghostRelaxationModulus`
- Nagai-Honda surface tensions: `gammaApical`, `gammaBasal`, `gammaLateral`
- Adhesion parameters: `nhMembraneSurface`, `nhCellCellAdhesion`, `nhBoundaryAdhesion`, per-type adhesions

Pressure/force-like parameters (Pa) are multiplied by `pressureConversionFactor`:
- `lumenPressure`, `apicalConstrictionStrength`
- `polarityBendingStrength`

### Examples

#### Example 1: Non-dimensional (backward compatible, default)

```json
{
  "unitConversion": {
    "usePhysicalUnits": false,
    "cellDiameterMicrometers": 10.0,
    "tissueViscosityPaS": 1.0
  },
  "ECM": {
    "ecmStiffness": 5.0
  }
}
```

Simulation uses `ecmStiffness = 5.0` directly (no conversion).

#### Example 2: Physical units (10 µm cells, 1 Pa·s tissue)

```json
{
  "unitConversion": {
    "usePhysicalUnits": true,
    "cellDiameterMicrometers": 10.0,
    "tissueViscosityPaS": 1.0
  },
  "ECM": {
    "ecmStiffness": 1000.0
  }
}
```

With these parameters:
- `dragCoeff = 6π × 1.0 × (5×10^-6) ≈ 9.42×10^-5 Pa·s·m`
- `stiffnessConversionFactor = 9.42×10^-5 / 10×10^-6 ≈ 9.42`
- `ecmStiffness` in simulation units = `1000.0 × 9.42 ≈ 9420`

This represents a 1 kPa (1000 Pa) ECM in physical units.

#### Example 3: Softer tissue (softer cells, higher viscosity)

```json
{
  "unitConversion": {
    "usePhysicalUnits": true,
    "cellDiameterMicrometers": 15.0,
    "tissueViscosityPaS": 5.0
  },
  "ECM": {
    "ecmStiffness": 500.0
  }
}
```

Larger cells and higher viscosity increase the drag coefficient, which scales all parameters up.

## Feature Toggles → C++ Classes

| Parameter | C++ Class | Notes |
|-----------|-----------|-------|
| `enableLumenPressure` | `LumenPressureForce` | Outward pressure from lumen interior |
| `enableApicalConstriction` | `ApicalConstrictionForce` | Vertex models only |
| `enableEcmGuidance` | `DynamicECMContactGuidanceForce` | 3D only; enables ECM remodeling |
| `enableRelaxation` | *(no force)* | Pauses proliferation during initial equilibration |
| `enableSloughing` | `PlaneBasedCellKiller` | Via `AddBoundingBoxKillers()` |
| `enableDifferentialAdhesion` | `DifferentialAdhesionForce` | Cell-type-dependent adhesion |
| `enableCurvatureBending` | `RingSmoothingForce` (2D) / `CurvatureBendingForce` (3D) | Surface smoothing |
| `enableCellPolarity` | `CellPolarityForce` | Apicobasal polarity maintenance |
| `enableEcmConfinement` | `ECMConfinementForce` + `DynamicECMField` | Fiber-based ECM confinement |
| `enableContinuousPvd` | `ContinuousPvdModifier` | Live ParaView preview during simulation |
| `useTopologyBasedSprings` | `RingSpringForce` (true) / `GeneralisedLinearSpringForce` (false) | Node models only |

## Force Parameters → C++ Classes

| INI Section | C++ Class(es) | Description |
|-------------|---------------|-------------|
| `[SpringForce]` | `RingSpringForce`, `GeneralisedLinearSpringForce`, `DifferentialAdhesionForce` | Cell-cell spring interactions |
| `[LumenPressureForce]` | `LumenPressureForce` | Hydrostatic lumen pressure |
| `[RingSmoothingForce]` | `RingSmoothingForce` (node2d), `SurfaceSmoothingForce` (node3d) | Discrete Laplacian smoothing |
| `[ECMConfinementForce]` | `ECMConfinementForce`, `DynamicECMField` | ECM fiber field and confinement |
| `[CellPolarityForce]` | `CellPolarityForce` | Polarity-based forces |
| `[ApicalConstrictionForce]` | `ApicalConstrictionForce` | Apical surface contractility |
| `[DifferentialAdhesionForce]` | `DifferentialAdhesionForce` | Cell-type adhesion differences |
| `[ECMGuidanceForce]` | `DynamicECMContactGuidanceForce` | 3D ECM-guided migration |
| `[VertexModel]` | `FastNagaiHondaForce` | Nagai-Honda vertex model forces |

## Cell Cycle → C++ Classes

| Parameter(s) | C++ Class | Description |
|-------------|-----------|-------------|
| `enableStochasticFourType = true` | `StochasticFourTypeCellCycleModel` | Stochastic SC/TA/PC/EC transitions (Montes-Olivas et al. 2023) with contact inhibition |
| `enableGenerationalCascade = true` | `UniformContactInhibitionGenerationalCellCycleModel` | Generational cascade (Meineke et al. 2001) — legacy, disabled by default |
| *(both false)* | `UniformContactInhibitionCellCycleModel` | Simple contact inhibition only |

### Stochastic 4-Type Model Parameters

| INI Parameter | Default | Description |
|--------------|---------|-------------|
| `probStemToStem` | 0.89 | P(SC daughter = SC) |
| `probStemToPaneth` | 0.09 | P(SC daughter = PC); P(SC→TA) = 1 − p_sc − p_pc |
| `probTaToTaEarly` | 0.9 | P(TA daughter = TA) for t < `transitionTime` |
| `probTaToTaLate` | 0.7 | P(TA daughter = TA) for t ≥ `transitionTime` |
| `transitionTime` | 120.0 | Hours at which TA probabilities switch (day 5) |
| `panethFraction` | 0.09 | Initial fraction of Paneth cells |

Cell types are encoded via `AbstractCellMutationState` subclasses:
- **Stem (SC)**: `WildTypeCellMutationState` + `TransitCellProliferativeType`
- **Transit-Amplifying (TA)**: `TACellMutationState` + `TransitCellProliferativeType`
- **Paneth (PC)**: `PanethCellMutationState` + `DifferentiatedCellProliferativeType` (non-proliferative)
- **Enterocyte (EC)**: `EnterocyteCellMutationState` + `DifferentiatedCellProliferativeType` (non-proliferative)

## Other Components

| Component | C++ Class | Description |
|-----------|-----------|-------------|
| Ring topology tracking | `RingTopologyTracker` | Tracks left/right neighbors on 2D ring |
| Surface topology tracking | `SurfaceTopologyTracker` | Tracks neighbors on 3D surface |
| ECM field visualization | `ECMFieldWriter`, `ECMFieldWriter3d` | VTK output of ECM density |
| Ring outline visualization | `RingOutlineWriter` | VTP polydata for 2D ring boundary |
| Cell polarity output | `CellPolarityWriter` | Writes polarity vectors to VTK |
| Division rules (node) | `TangentialCentreBasedDivisionRule` | Tangential division on ring/surface |
| Division rules (vertex) | `LocalTangentVertexBasedDivisionRule` | Local tangent-aligned vertex division |

## Source File Locations

All custom classes are in:
- `TissueMorphology/src/` — Force classes, modifiers, writers
- `TissueMorphology/apps/src/` — App-specific utilities and runners

Chaste built-in classes are in:
- `Chaste/cell_based/src/population/forces/` — `GeneralisedLinearSpringForce`, etc.
- `Chaste/cell_based/src/population/killers/` — `PlaneBasedCellKiller`

## Parameter Autotuning (CryptBuddingParams::Finalise)

Several parameters are **auto-derived** in `Finalise()` and will silently override `SetDefaults()` values unless you explicitly set them in your `.ini` file.  Setting a parameter in the INI file prevents the override.

### dt (time step)

| Model | Default dt | Condition |
|-------|-----------|-----------|
| `node2d` | `0.005` | Always |
| `vertex2d` | `0.0002` | `ecmStiffness < 1.0` |
| `vertex2d` | `0.0005` | `ecmStiffness ≥ 1.0` |
| `node3d` | `0.01` | Always |
| `vertex3d` | `0.0001` | Always |

### endTime

Defaults to `168.0` hours for `vertex2d`, `node3d`, and `vertex3d` if not set in INI.

### t1Threshold2d (vertex T1 swap threshold)

| Condition | Default |
|-----------|---------|
| `ecmStiffness < 2.0` | `0.2` |
| `ecmStiffness ≥ 2.0` | `0.15` |

Lower values delay T1 swaps (more tolerant of short edges). Higher stiffness → lower threshold because stiffer ECM compresses edges less.

### t2Threshold2d (vertex T2 removal threshold)

Always defaults to `0.05` if not set.

### How to prevent autotuning

Set `dt`, `endTime`, `t1Threshold2d`, and `t2Threshold2d` explicitly in your `.ini` file. When present in the INI, `Finalise()` detects the `*Overridden` flags and skips the auto-derivation.
