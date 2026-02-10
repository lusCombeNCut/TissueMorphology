# ParaView ECM Grid Visualization Guide

## Quick Start (TL;DR)

**The ECM data is in a separate PVD file called `ecm_results.pvd`!**

1. **Load cells**: File → Open → `results.pvd` → Apply
2. **Load ECM**: File → Open → `ecm_results.pvd` → Apply
3. **Add arrows**: Select ecm_results.pvd → Filters → Glyph
   - Orientation Array: **ecm_orientation** (it's already a vector!)
   - Scale Factor: 30-50
   - Apply
4. Done! Both time series should advance together.

---

## Invasion Analysis Plots

To analyze and plot cell invasion over time for all ECM configurations, run:

```bash
cd Chaste/projects/TissueMorphology
python3 plot_invasion_analysis.py
```

This generates:
- **invasion_analysis.png**: Comprehensive plots showing:
  - Invasion front position over time
  - 90th percentile cell depth
  - Invasion speed over time
  - Cell population growth
  - Final invasion comparison bar chart
- **invasion_summary.txt**: Quantitative metrics for each ECM type

The script compares perpendicular, parallel, random, and mixed ECM fiber orientations.

---

## ECM Simulation Snapshots

To export first and last timestep comparison images:

```bash
cd Chaste/projects/TissueMorphology
python3 export_ecm_snapshots.py
```

This generates:
- **Individual comparison images** (e.g., `perpendicular_comparison.png`): Side-by-side initial vs final state for each ECM type
- **all_ecm_comparison.png**: 4×2 grid showing all ECM configurations together
- Images include both cell positions (colored dots) and ECM fiber orientation (gray arrows)
- Output saved to `ecm_comparison_output/`

---

## Problem: Why Can't I See ECM Fields?

The ECM fiber orientation data is written to **separate VTK files** grouped in `ecm_results.pvd` by the `ECMFieldWriter` class, not to the main `results.pvd` file that contains cell data. 

**Common mistake**: Looking for ECM fields in `results.pvd` - they're not there! The ECM orientation is stored as a **VECTORS** field called `ecm_orientation` in the separate `ecm_results.pvd` file.

## Solution: Load ECM Grid Files Directly

The ECM data is written to separate `.vtk` files (not in `results.pvd`) by the `ECMFieldWriter`. You need to load these files separately.

### Method 1: Load ECM Grid Files (Easiest)

The simulation writes ECM grid files to the output directory as `ecm_grid_*.vtk` and groups them in `ecm_results.pvd`.

#### Step-by-Step:

1. **Load cell data**
   - File → Open → `results.pvd`
   - Apply
   - This shows cells only

2. **Load ECM grid data**
   - File → Open → Navigate to same output directory
   - Select `ecm_results.pvd`
   - Click OK
   - Apply
   - This shows ECM grid data as a time series

3. **Verify ECM data fields**
   - Select the ecm_results.pvd in Pipeline Browser
   - Look at Information tab - you should see:
     - **ecm_density** (scalar)
     - **ecm_anisotropy** (scalar)  
     - **ecm_orientation** (vector) ← This is what we need for arrows!

4. **Add arrow glyphs**
   - Select the ecm_results.pvd
   - Filters → Alphabetical → **Glyph**
   - In Properties:
     - **Glyph Type**: Arrow (or 2D Glyph → Arrow)
     - **Orientation Array**: **ecm_orientation** (NOT ecm_orientation_x!)
     - **Scale Array**: No scale array (constant size)
       - Or use **ecm_density** or **ecm_anisotropy** to scale by ECM strength
     - **Scale Factor**: 30-50 (adjust to see arrows clearly)
     - **Glyph Mode**: All Points
     - **Coloring**: 
       - Solid Color (gray) for simple view
       - Or **ecm_density** to show ECM degradation
       - Or **ecm_anisotropy** to show fiber alignment
   - Click **Apply**

5. **Synchronize time with cells**
   - Both `results.pvd` and `ecm_results.pvd` should advance together
   - Use the time controls at top to animate
   - NOTE: ECM grid files are written every 120 timesteps (60 minutes)
     - This matches the VTU file output frequency

---

### Alternative Method: If ECM Fields Don't Appear

If you don't see the ECM grid files or fields, the ECM data might not have been written. Check:

1. **Verify ECM files exist**
   ```bash
   ls testoutput/DynamicECMInvasion2d/perpendicular/ecm_results.pvd
   ls testoutput/DynamicECMInvasion2d/perpendicular/ecm_grid_*.vtk
   ```
   You should see `ecm_results.pvd` and files like: `ecm_grid_0.vtk`, `ecm_grid_120.vtk`, etc.

2. **If files don't exist**, make sure your simulation includes:
   ```cpp
   // In your test file
   MAKE_PTR_ARGS(ECMFieldWriter<2>, p_ecm_writer, (p_ecm_field, 120));
   simulator.AddSimulationModifier(p_ecm_writer);
   ```

3. **Check one of the VTK files** to verify content:
   ```bash
   head -50 testoutput/DynamicECMInvasion2d/perpendicular/ecm_grid_0.vtk
   ```
   You should see data sections for:
   - `SCALARS ecm_density`
   - `SCALARS ecm_anisotropy`
   - `VECTORS ecm_orientation`

---

### Method 2: Using Resample Filter (If you want ECM on cell grid)

**Note:** This method only works if ECM fields are stored WITH the cell data in `results.pvd`. 
Since your ECM is in separate files, use Method 1 instead.

<details>
<summary>Click to expand legacy instructions</summary>

This creates a regular 3D grid and samples ECM values onto it.

#### Step-by-Step:

1. **Load your data**
   - File → Open → `results.pvd`
   - Apply

2. **Create a uniform grid**
   - Select your data in Pipeline Browser
   - Filters → Alphabetical → **Resample To Image**
   - In Properties panel:
     - **Sampling Bounds**: Use Input Bounds (should be automatic)
     - **Sampling Dimensions**: Set to `[13, 21, 1]` for 2D
       - This creates grid with ~50 µm spacing (600/13 ≈ 46 µm, 1000/21 ≈ 48 µm)
       - For finer grid: `[21, 41, 1]` (~30 µm spacing)
       - For coarser: `[7, 13, 1]` (~85 µm spacing)
   - Click **Apply**

3. **Create vectors from components (ONLY if needed)**
   - **Check first**: Does your data have `ecm_orientation` as a vector field?
   - If YES: Skip this step - go directly to step 4
   - If NO (you have separate x/y/z components):
     - Select the ResampleToImage output
     - Filters → Alphabetical → **Calculator**
     - **Result Array Name**: `ECM_Vector`
     - **Expression**: `ecm_orientation_x*iHat + ecm_orientation_y*jHat + ecm_orientation_z*kHat`
     - Click **Apply**

4. **Add arrow glyphs**
   - Select the Calculator output
   - Filters → Alphabetical → **Glyph**
   - In Properties:
     - **Glyph Type**: Arrow (or 2D Glyph → Arrow)
     - **Orientation Array**: ECM_Vector
     - **Scale Array**: No scale array (constant size)
     - **Scale Factor**: 30-50 (adjust to see arrows clearly)
     - **Glyph Mode**: All Points
     - **Coloring**: Solid Color (gray) or by ecm_angle
   - Click **Apply**

5. **Overlay cells**
   - Make original data visible again (eye icon in Pipeline Browser)
   - Color cells differently (e.g., blue)
   - Adjust cell size/representation

</details>

---

### Method 3: Using Programmable Source (Advanced)

For more control over grid placement:

1. **Create custom grid**
   - Sources → **Programmable Source**
   - In **Script** field, paste:

```python
import numpy as np

# Grid parameters
x_min, x_max = 0, 600
y_min, y_max = 0, 1000
spacing = 50  # Grid spacing in µm

# Create grid
x = np.arange(x_min, x_max + spacing, spacing)
y = np.arange(y_min, y_max + spacing, spacing)
X, Y = np.meshgrid(x, y)

# Flatten for VTK
points = np.column_stack([X.ravel(), Y.ravel(), np.zeros(X.size)])
n_points = len(points)

# Create output
output.SetPoints(vtk.vtkPoints())
for p in points:
    output.GetPoints().InsertNextPoint(p)

# Add ECM orientation vectors
# (You'd need to compute these based on ECM type)
ecm_type = 'perpendicular'  # Change as needed

if ecm_type == 'perpendicular':
    vectors = np.column_stack([np.zeros(n_points), 
                               np.ones(n_points), 
                               np.zeros(n_points)])
elif ecm_type == 'parallel':
    vectors = np.column_stack([np.ones(n_points), 
                               np.zeros(n_points), 
                               np.zeros(n_points)])
# etc.

arr = vtk.vtkFloatArray()
arr.SetName("ECM_Direction")
arr.SetNumberOfComponents(3)
for v in vectors:
    arr.InsertNextTuple(v)
output.GetPointData().AddArray(arr)
```

2. Then apply **Glyph** filter to this source

---

### Method 4: Export Grid from Python and Load in ParaView

Create a separate VTK file with just the ECM grid:

```python
import pyvista as pv
import numpy as np

# Create uniform grid
x = np.arange(0, 601, 50)
y = np.arange(0, 1001, 50)
z = [0]

grid = pv.RectilinearGrid(x, y, z)

# Compute ECM orientation at each point
points = grid.points
ecm_vectors = np.zeros((len(points), 3))

ecm_type = 'perpendicular'  # Change as needed

for i, p in enumerate(points):
    if ecm_type == 'perpendicular':
        ecm_vectors[i] = [0, 1, 0]
    elif ecm_type == 'parallel':
        ecm_vectors[i] = [1, 0, 0]
    elif ecm_type == 'random':
        angle = np.random.uniform(0, 2*np.pi)
        ecm_vectors[i] = [np.cos(angle), np.sin(angle), 0]
    elif ecm_type == 'mixed':
        x_pos = p[0]
        if 200 < x_pos < 400:  # Center region
            ecm_vectors[i] = [0, 1, 0]  # Perpendicular
        else:
            ecm_vectors[i] = [1, 0, 0]  # Parallel

grid['ECM_Direction'] = ecm_vectors

# Save
grid.save('ecm_grid_perpendicular.vtu')
```

Then in ParaView:
- Open `ecm_grid_perpendicular.vtu`
- Apply Glyph filter with ECM_Direction as orientation

## Recommended Settings

### For Best Visualization:

- **Grid spacing**: 40-60 µm (balances detail vs clutter)
- **Arrow scale**: 30-50 µm
- **Arrow color**: Light gray (so cells stand out)
- **Cell representation**: Points or Spheres, blue color
- **Background**: White or light gray

### Multiple Representations:

Show both:
1. **ECM grid** (gray arrows, regular grid)
2. **Cells** (blue spheres at actual positions)
3. **Cells with migration vectors** (red arrows from cell positions)

To do this, create separate Glyph filters for ECM grid and cell migration.

## Quick Comparison Setup

To compare all 4 ECM types side-by-side:

1. Load all 4 scenarios in separate ParaView instances
2. Use **Link Views** (Tools → Link Camera)
3. Arrange windows in 2×2 grid
4. Synchronize time slider

Or use Python trace to automate this!

## Troubleshooting

**Arrows too small/large?**
- Adjust Scale Factor in Glyph properties
- Try values 20-100 depending on domain size

**Grid too dense?**
- Reduce Sampling Dimensions in ResampleToImage
- Or increase spacing in programmable source

**Can't see ECM vectors?**
- Check that ecm_orientation_x/y/z fields exist in data
- Verify Calculator expression is correct
- Make sure Glyph orientation is set to ECM_Vector

**Want colored arrows by angle?**
- In Glyph properties, set Coloring to ecm_angle
- Apply color map (e.g., HSV for cyclic data)
