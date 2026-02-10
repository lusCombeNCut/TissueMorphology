# Invasive Front Simulation - Three-Factor Model Application

## Overview

This implements Section 3.1 "Invasive Cellular Front Pushing into ECM" from Metzcar et al. 2025, demonstrating how the three-factor morphogenesis framework can be adapted to study tumor invasion.

## Three-Factor Model Configuration

### 1. **ECM Contact Guidance Force** (replaces BasementMembraneForce)
   - **Role**: Structured ECM provides directional guidance to migrating cells
   - **Parameters**:
     - Base speed: v₀ = 1.25 µm/min
     - ECM sensitivity: s = 1.0 (full response to ECM)
     - Anisotropy: a = 1.0 (highly aligned fibers)
   - **Biological interpretation**: ECM fibers act as "highways" for cell migration
   - **Implementation**: `ECMContactGuidanceForce.hpp`

### 2. **Differential Adhesion Force**
   - **Role**: Cell-cell interactions in the invasive front
   - **Parameters**:
     - Adhesion strength: 0.4 (all cells same type)
     - Cutoff distance: 50 µm
   - **Biological interpretation**: Cells maintain contact but can push apart
   - **Implementation**: `DifferentialAdhesionForce.hpp` (uniform adhesion)

### 3. **Apical Constriction Force** (not used)
   - **Role**: Not applicable to invasion scenario
   - **Reason**: Invasive front cells are not polarized epithelium
   - **Status**: Force present but not applied

## ECM Orientation Scenarios

Following Painter 2009, we test four ECM patterns:

### 1. Random Orientation
- **Description**: Isotropic, randomly oriented ECM fibers
- **Expected**: Diffusive invasion pattern (baseline)
- **Biological relevance**: Unstructured/degraded ECM

### 2. Parallel Orientation (θ = 0°)
- **Description**: Fibers aligned horizontally (parallel to invasion front)
- **Expected**: Slowest invasion (fibers block perpendicular movement)
- **Biological relevance**: Restrictive ECM barrier

### 3. Perpendicular Orientation (θ = 90°)
- **Description**: Fibers aligned vertically (perpendicular to front)
- **Expected**: Fastest invasion (fibers guide cells forward)
- **Biological relevance**: Pre-aligned invasion tracks

### 4. Mixed Orientation
- **Description**: Perpendicular in center, parallel on sides
- **Expected**: Heterogeneous invasion (fast center, slow sides)
- **Biological relevance**: Spatially variable ECM structure

## Key Results from Paper

From Metzcar 2025:
> "Consistent with the work by Painter, cells traveled further in matrices with perpendicular and random ECM orientations compared to the parallel setting... the fastest advances of the invasive front occurred when the ECM orientations were perpendicular to the front followed by the random orientations and finally parallel orientations."

**Invasion Speed Ranking**: Perpendicular > Random > Parallel

## Simulation Parameters

| Parameter | Value | Units | Source |
|-----------|-------|-------|--------|
| Domain size | 600 × 1000 | µm | Paper Section 3.1 |
| Initial cells | 30 | cells | Bottom edge |
| Cell influx | 30 cells / 180 min | cells/min | (Not yet implemented) |
| Simulation time | 5 days | days | (Testing with 1 day) |
| Base cell speed | 1.25 | µm/min | Paper |
| Cell-cell adhesion | 0.4 | - | Paper |
| Repulsion | 25.0 | µm | Paper |
| ECM sensitivity | 1.0 | - | Full response |
| ECM anisotropy | 1.0 | - | Highly aligned |

## Quantification Metrics

To match the paper's analysis:

1. **Invasion front position**: 95th percentile of cell y-coordinates
2. **Cell distribution**: Histogram of cell positions (40 bins of 25 µm)
3. **Invasion speed**: Distance traveled over 5 days
4. **Front shape**: Cell count contours to visualize diffusive vs superdiffusive spread

## Implementation Notes

### Completed
- ✅ ECMContactGuidanceForce class with 4 orientation patterns
- ✅ Integration with DifferentialAdhesionForce
- ✅ Test harness for all 4 scenarios
- ✅ Basic parameter setup matching paper

### To Do
- ⏳ Cell influx modifier (30 cells every 180 minutes at boundary)
- ⏳ Extended simulation to 5 days (currently 1 day for testing)
- ⏳ Quantification of invasion front position
- ⏳ Stochastic replicates for statistical analysis
- ⏳ Visualization comparing 4 scenarios side-by-side

## Differences from Three-Factor Morphogenesis

| Aspect | Organoid Model | Invasion Model |
|--------|---------------|----------------|
| **ECM role** | Radial constraint (basement membrane) | Directional guidance (fiber orientation) |
| **Cell types** | Apical vs Basal | Uniform invasive front |
| **Geometry** | Spherical/circular organoid | Linear invasion front |
| **Forces** | AC + DA + BM | ECM guidance + DA |
| **Output** | Tissue folding/buckling | Invasion depth/pattern |

## Biological Context

This simulation addresses:
- **Tumor invasion**: Cancer cells infiltrating surrounding tissue
- **ECM structure**: How matrix architecture affects invasion
- **Contact guidance**: Cells following ECM fiber alignment
- **Collective migration**: Cell groups moving together with adhesion

The framework demonstrates that the same computational infrastructure (three-factor model) can represent diverse morphogenetic phenomena by reconfiguring force parameters and spatial patterns.

## Files

- **Test**: `Test2dInvasiveFront.hpp`
- **Force**: `ECMContactGuidanceForce.hpp`
- **Output**: `testoutput/InvasiveFront2d/{random,parallel,perpendicular,mixed}/`
