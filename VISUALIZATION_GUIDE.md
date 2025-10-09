# Organoid Formation Visualization Guide

This document explains the visualization capabilities added to the TissueMorphology project for analyzing organoid formation simulations.

## Overview

The organoid formation tests now generate VTK-compatible output files that can be visualized using Python scripts. The visualization system captures:

- Cell positions over time
- Cell ages (represented by colors)
- Effects of different basement membrane stiffness values
- Dynamic morphological changes during organoid formation

## Generated Files

### Test Output Data
The tests generate the following visualization data in `/home/orlando/Thesis/testoutput/OrganoidFormation/`:

- `BasicTest/`: Basic organoid formation without basement membrane forces
- `LowStiffness/`: Organoid formation with low basement membrane stiffness (0.1)
- `HighStiffness/`: Organoid formation with high basement membrane stiffness (10.0)  
- `WithBasementMembrane/`: Organoid formation with basement membrane force applied

Each test directory contains:
- `cellages.dat`: Cell position and age data over time
- `results.pvd`: ParaView data file for 3D visualization
- Additional Chaste visualization files

### Visualization Scripts

1. **`visualize_chaste_results.py`**
   - Main visualization script for Chaste output data
   - Supports both static plots and animations
   - Usage: `python3 visualize_chaste_results.py [data_file] [options]`

2. **`create_visualizations.py`**  
   - Comprehensive script creating multiple comparison figures
   - Generates side-by-side scenario comparisons
   - Creates time series evolution plots

### Generated Visualization Files

1. **`organoid_comparison.png`**
   - 2x2 grid comparing all four test scenarios
   - Shows final cell configurations at end times
   - Color-coded by cell age

2. **`organoid_time_series.png`**
   - Time evolution of basement membrane scenario
   - Shows progression from t=0.0 to t=0.2
   - Demonstrates morphological changes over time

3. **`organoid_formation.gif`**
   - Animated visualization of organoid formation
   - Shows dynamic cell movement and aging
   - 5-frame animation covering simulation time course

4. **Individual scenario plots**
   - `low_stiffness.png`: Low basement membrane stiffness results
   - `high_stiffness.png`: High basement membrane stiffness results  
   - `organoid_t0.2.png`: Basement membrane scenario at t=0.2

## Visualization Features

### Color Coding
- **Cell Age**: Represented by color gradient (viridis colormap)
  - Darker colors: Younger cells
  - Brighter colors: Older cells
  - Colorbar shows age scale in simulation time units

### Plot Elements
- **Scatter points**: Individual cells positioned by (x,y) coordinates
- **Grid**: Background grid for spatial reference
- **Annotations**: Cell count and time point information
- **Equal aspect ratio**: Maintains proper spatial proportions

## Scientific Insights

### Basement Membrane Effects
The visualizations reveal several key morphological features:

1. **Spatial Organization**: Basement membrane forces promote more compact, organized cell arrangements
2. **Age Distribution**: Older cells tend to remain in central positions while newer cells are positioned peripherally
3. **Stiffness Effects**: Higher basement membrane stiffness leads to tighter cell clustering
4. **Dynamic Evolution**: Cell positions change significantly over time, showing active morphological remodeling

### Comparison Analysis
- **Basic vs. Membrane Tests**: Clear differences in final organization patterns
- **Stiffness Variations**: Dose-dependent effects of basement membrane stiffness on morphology
- **Time Evolution**: Progressive changes in cell arrangement and tissue architecture

## Usage Instructions

### Running Individual Visualizations
```bash
# Static plot of specific time point
python3 visualize_chaste_results.py [data_file] --time 0.2 --output plot.png

# Create animation  
python3 visualize_chaste_results.py [data_file] --animate --output animation.gif

# Show available time points
python3 visualize_chaste_results.py [data_file]
```

### Creating Comparison Figures
```bash
# Generate all comparison visualizations
python3 create_visualizations.py
```

### Viewing Results
The visualization files can be viewed with any image viewer or incorporated into presentations and publications. The PNG files are high-resolution (300 DPI) suitable for scientific publication.

## Technical Notes

### Data Format
- Chaste cellages.dat format: `time cell_id x_position y_position age`
- Multiple time points in single file separated by time stamps
- Parser handles mixed line formats automatically

### Dependencies
- Python 3 with matplotlib, numpy
- PIL (for GIF animation export)
- Standard Python libraries (os, sys, argparse)

### Customization
The visualization scripts can be modified to:
- Change color schemes (modify colormap parameter)
- Adjust plot sizes and layouts
- Add additional cell properties for visualization
- Export different file formats