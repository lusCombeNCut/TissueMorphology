# ECM Visualization Guide for Invasive Front Simulations

## Overview

This guide explains how to visualize the ECM (Extracellular Matrix) properties in the invasive front simulations, both using ParaView and custom Python scripts.

## ECM Data Stored in Simulations

The `ECMContactGuidanceForce` automatically stores ECM and cell migration information in the VTK output files:

### Cell Data Fields:
- **ecm_angle**: ECM fiber orientation angle in radians at each cell's position
- **ecm_orientation_x**: X-component of ECM fiber direction vector
- **ecm_orientation_y**: Y-component of ECM fiber direction vector  
- **ecm_orientation_z**: Z-component (0 for 2D simulations)
- **migration_direction_x**: X-component of actual cell migration direction
- **migration_direction_y**: Y-component of actual cell migration direction
- **migration_direction_z**: Z-component (0 for 2D)

## Method 1: Python Visualization (Recommended)

### Quick Start

Run the simple visualization script:

```bash
cd /home/orlando/Thesis
python3 visualize_invasive_front_simple.py
```

This generates:
- `invasion_comparison.png` - 4-panel comparison of all ECM scenarios
- `invasion_timeseries_*.png` - Time series plots for each scenario
- `invasion_*_final.png` - Final snapshots with ECM fiber field overlaid

### Features

The Python script visualizes:
1. **Background ECM fiber field** - Gray lines showing fiber orientations
2. **Cell positions** - Blue dots for individual cells
3. **Invasion depth** - Red dashed line at 95th percentile
4. **Time series** - Invasion depth and cell count over time

### ECM Patterns Visualized

- **Random**: Isotropic, randomly oriented fibers (diffusive invasion)
- **Parallel**: Horizontal fibers (restrictive, slowest invasion)
- **Perpendicular**: Vertical fibers (facilitative, fastest invasion)
- **Mixed**: Perpendicular in center, parallel on sides (heterogeneous)

## Method 2: ParaView Visualization

### Loading Data

1. Open ParaView
2. File → Open → Navigate to results directory:
   ```
   /home/orlando/Thesis/testoutput/InvasiveFront2d/perpendicular/results_from_time_0/
   ```
3. Select `results.pvd` (this loads all timepoints)
4. Click **Apply**

### Visualizing ECM Orientation Vectors

1. **Select the data** in Pipeline Browser
2. **Filters** → **Common** → **Glyph**
3. In Glyph properties:
   - **Glyph Type**: Arrow
   - **Orientation Array**: Select **ecm_orientation** (computed from x,y,z components)
   - **Scale Array**: Set to constant (e.g., 10)
   - **Glyph Mode**: All Points
4. Click **Apply**

### Visualizing Migration Directions

1. **Select the data** again
2. **Filters** → **Common** → **Glyph** (create second glyph filter)
3. In Glyph properties:
   - **Glyph Type**: Arrow
   - **Orientation Array**: **migration_direction**
   - **Scale Array**: constant (e.g., 15)
   - **Color**: Set to red (differentiate from ECM)
4. Click **Apply**

### Color-Coding by ECM Angle

1. Select the main data object (not glyphs)
2. In **Coloring** dropdown, select **ecm_angle**
3. Choose a color map (e.g., "Cool to Warm")
4. This shows spatial variation in ECM orientation

### Creating Animations

1. **View** → **Animation View**
2. Set time range (0 to 1440 minutes)
3. Click play button to animate invasion over time
4. **File** → **Save Animation** to export as video

## Method 3: Custom Analysis Scripts

### Extract ECM Data from VTU Files

To work with raw ECM data, use PyVista:

```python
import pyvista as pv
import numpy as np

# Load VTU file
mesh = pv.read('/path/to/results_1440.vtu')

# Extract data
positions = mesh.points  # Cell positions (Nx3 array)
ecm_angles = mesh.point_data['ecm_angle']  # ECM angles
ecm_x = mesh.point_data['ecm_orientation_x']
ecm_y = mesh.point_data['ecm_orientation_y']
mig_x = mesh.point_data['migration_direction_x']
mig_y = mesh.point_data['migration_direction_y']

# Calculate alignment between ECM and migration
alignment = ecm_x * mig_x + ecm_y * mig_y  # Dot product
print(f"Mean alignment: {np.mean(alignment):.3f}")  # 1.0 = perfect alignment
```

## Understanding the Visualizations

### ECM Fiber Patterns

- **Gray lines in background**: Show local ECM fiber orientation
- **Line direction**: Indicates fiber alignment
- **Pattern consistency**: 
  - Random: chaotic, no preferred direction
  - Parallel: all horizontal
  - Perpendicular: all vertical
  - Mixed: transitions between regions

### Cell Behavior

- **Blue dots**: Individual cell positions
- **Vertical spread**: Invasion depth into domain
- **Horizontal clustering**: Cell-cell adhesion effects
- **Leading edge**: Most advanced cells (top of distribution)

### Quantitative Metrics

From time series plots:

1. **95th percentile depth**: Position of leading edge
   - Perpendicular ECM: deepest invasion
   - Parallel ECM: shallowest invasion
   
2. **Mean position**: Center of mass of cell population
   - Tracks bulk invasion speed
   
3. **Cell count**: Should remain constant (30 cells, no proliferation)

## Expected Results

Based on Metzcar et al. 2025 and Painter 2009:

| ECM Type | Invasion Speed | Pattern | Biological Relevance |
|----------|---------------|---------|---------------------|
| **Perpendicular** | Fastest | Superdiffusive spread | Pre-aligned invasion tracks |
| **Random** | Medium | Diffusive (Fisher-like) | Unstructured/degraded ECM |
| **Parallel** | Slowest | Sharp boundary | ECM barrier effect |
| **Mixed** | Heterogeneous | Fast center, slow sides | Spatially variable matrix |

### Invasion Depth Ranking (at 1 day)

Expected: **Perpendicular > Random > Parallel**

Check this in the comparison plot or time series.

## Troubleshooting

### No ECM data in ParaView?

- Verify you're loading results from latest simulation run
- Check that `ECMContactGuidanceForce` was applied (not just `GeneralisedLinearSpringForce`)
- Rebuild and rerun test if force was recently modified

### Python script fails?

```bash
# Install required packages
pip install numpy matplotlib

# Check data exists
ls /home/orlando/Thesis/testoutput/InvasiveFront2d/*/results_from_time_0/results.viznodes
```

### VTU files compressed?

The VTU files use zlib compression. Use PyVista to read them:

```bash
pip install pyvista
```

Then use the advanced `visualize_invasive_front.py` script (handles compressed VTU).

## Files

- **ECMContactGuidanceForce.hpp**: Force class that computes and stores ECM data
- **visualize_invasive_front_simple.py**: Main visualization script (uses .viznodes)
- **visualize_invasive_front.py**: Advanced script (parses compressed .vtu files, requires PyVista)
- **Test2dInvasiveFront.hpp**: Simulation test that applies ECM force

## Advanced: Create Custom Metrics

Example: Calculate ECM-migration alignment for each scenario:

```python
for ecm_type in ['random', 'parallel', 'perpendicular', 'mixed']:
    mesh = pv.read(f'testoutput/InvasiveFront2d/{ecm_type}/results_from_time_0/results_1440.vtu')
    
    # Dot product of normalized vectors
    ecm_x = mesh.point_data['ecm_orientation_x']
    ecm_y = mesh.point_data['ecm_orientation_y']
    mig_x = mesh.point_data['migration_direction_x']
    mig_y = mesh.point_data['migration_direction_y']
    
    alignment = ecm_x * mig_x + ecm_y * mig_y
    
    print(f"{ecm_type:15s}: alignment = {np.mean(alignment):.3f} ± {np.std(alignment):.3f}")
```

Expected output:
- Perpendicular & Parallel: high alignment (~1.0)
- Random: moderate alignment (~0.5-0.7)
- Mixed: bimodal distribution

This quantifies how well cells follow ECM guidance!
