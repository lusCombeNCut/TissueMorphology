#!/usr/bin/env python3
"""
Export first and last timestep images from Dynamic ECM Invasion simulations.

This script creates a 4x2 grid showing:
- 4 rows (one for each ECM type)
- 2 columns (left = initial state, right = final state)
- Each row shows the same ECM configuration at different times
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import xml.etree.ElementTree as ET
from typing import Dict, Tuple

# Font size constants for easy customization
FONT_SIZE_TITLE = 20
FONT_SIZE_LABEL = 18
FONT_SIZE_TICK = 18
FONT_SIZE_COLORBAR = 20
FONT_SIZE_LEGEND = 18
FONT_SIZE_SUPTITLE = 20

Y_LIM = 400
X_LIM = 600

def read_vtu_file(filepath: Path) -> Tuple[np.ndarray, Dict]:
    """Read cell positions from VTU file using VTK library."""
    try:
        import vtk
        from vtk.util.numpy_support import vtk_to_numpy
        
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(str(filepath))
        reader.Update()
        
        output = reader.GetOutput()
        points = output.GetPoints()
        
        if points:
            positions = vtk_to_numpy(points.GetData())[:, :2]
            return positions, {}
        else:
            return np.array([]), {}
            
    except ImportError:
        print("Warning: python3-vtk not available, cannot read VTU files")
        return np.array([]), {}
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return np.array([]), {}

def read_ecm_file(filepath: Path):
    """Read ECM grid data from VTK file (handles both legacy .vtk and .vti formats)."""
    try:
        import vtk
        from vtk.util.numpy_support import vtk_to_numpy
        
        # Determine reader based on file extension
        if filepath.suffix == '.vti':
            reader = vtk.vtkXMLImageDataReader()
        else:  # .vtk legacy format
            reader = vtk.vtkStructuredPointsReader()
        
        reader.SetFileName(str(filepath))
        reader.Update()
        
        output = reader.GetOutput()
        dims = output.GetDimensions()
        
        # Get ECM data arrays
        orientation = output.GetPointData().GetArray('ecm_orientation')
        density = output.GetPointData().GetArray('ecm_density')
        
        if orientation and density:
            # Reshape data to 2D grid
            orient_data = vtk_to_numpy(orientation).reshape(dims[1], dims[0], 3)
            density_data = vtk_to_numpy(density).reshape(dims[1], dims[0])
            
            # Get grid spacing and origin
            spacing = output.GetSpacing()
            origin = output.GetOrigin()
            
            # Create coordinate arrays
            x = np.linspace(origin[0], origin[0] + (dims[0]-1)*spacing[0], dims[0])
            y = np.linspace(origin[1], origin[1] + (dims[1]-1)*spacing[1], dims[1])
            X, Y = np.meshgrid(x, y)
            
            return X, Y, orient_data, density_data
        else:
            return None, None, None, None
            
    except ImportError:
        print("Warning: python3-vtk not available, cannot read VTK files")
        return None, None, None, None
    except Exception as e:
        print(f"Error reading VTK file {filepath}: {e}")
        return None, None, None, None

def plot_ecm_from_file(ax, ecm_file):
    """Plot ECM grid from VTK file with density heatmap and fiber orientations."""
    X, Y, orient, density = read_ecm_file(ecm_file)
    
    if X is None:
        return None
    
    # Plot density heatmap (background) - yellow to red colormap
    im = ax.pcolormesh(X, Y, density, cmap='viridis', alpha=0.2, 
                      vmin=0.0, vmax=1.0, shading='auto', zorder=0)
    
    # Subsample for arrow visualization
    skip = 2
    U = orient[::skip, ::skip, 0]
    V = orient[::skip, ::skip, 1]
    
    # Plot fiber orientations as arrows (shorter with scale=35)
    ax.quiver(X[::skip, ::skip], Y[::skip, ::skip], U, V,
             color='black', alpha=0.8, 
             scale=35, width=0.003,
             headwidth=4, headlength=5,
             zorder=1)
    
    return im

def read_pvd_file(pvd_path: Path) -> list:
    """Read PVD file and return list of (time, filepath) tuples."""
    tree = ET.parse(pvd_path)
    root = tree.getroot()
    
    timesteps = []
    for dataset in root.findall('.//DataSet'):
        time = float(dataset.get('timestep'))
        filename = dataset.get('file')
        file_path = pvd_path.parent / filename
        timesteps.append((time, file_path))
    
    return sorted(timesteps, key=lambda x: x[0])

def create_combined_comparison(base_dir: Path, output_dir: Path):
    """Create a 4x2 grid showing 4 ECM types (each row is a configuration, column 0 is t=0, column 1 is t=96h)."""
    
    ecm_types = ['perpendicular', 'parallel', 'random', 'mixed']
    # Use consistent blue color for all cells
    cell_color = 'blue'
    titles = {
        'perpendicular': 'Perpendicular',
        'parallel': 'Parallel',
        'random': 'Random',
        'mixed': 'Mixed'
    }
    
    # Create figure - 4 rows x 2 columns (each row is one configuration)
    fig = plt.figure(figsize=(12, 20))
    gs = fig.add_gridspec(4, 2, hspace=0.25, wspace=0.25)
    
    ecm_im = None
    for idx, ecm_type in enumerate(ecm_types):
        results_dir = base_dir / ecm_type
        
        # Check if we have multiple runs (run_12345, run_13345, etc.)
        run_dirs = sorted([d for d in results_dir.iterdir() if d.is_dir() and d.name.startswith('run_')])
        if run_dirs:
            # Use the first run directory
            results_dir = run_dirs[0]
        
        if (results_dir / "results_from_time_0").exists():
            results_dir = results_dir / "results_from_time_0"
        
        # Read data
        cell_pvd = results_dir / "results.pvd"
        if not cell_pvd.exists():
            continue
        
        cell_timesteps = read_pvd_file(cell_pvd)
        if len(cell_timesteps) < 2:
            continue
        
        ecm_pvd = results_dir / "ecm_results.pvd"
        ecm_timesteps = []
        if ecm_pvd.exists():
            ecm_timesteps = read_pvd_file(ecm_pvd)
        
        first_cells, _ = read_vtu_file(cell_timesteps[0][1])
        last_cells, _ = read_vtu_file(cell_timesteps[-1][1])
        
        # Each configuration gets a row (idx = 0, 1, 2, 3)
        row = idx
        
        # Initial state (column 0)
        ax = fig.add_subplot(gs[row, 0])
        ax.set_facecolor('#f8f9fa')
        
        # Plot ECM grid from VTI file
        if ecm_timesteps and len(ecm_timesteps) > 0:
            plot_ecm_from_file(ax, ecm_timesteps[0][1])
        
        # Plot cells
        if len(first_cells) > 0:
            ax.scatter(first_cells[:, 0], first_cells[:, 1], 
                      c=cell_color, s=30, alpha=0.7, edgecolors='darkblue', linewidth=1, zorder=3)
            
            # Add invasion front line (maximum)
            invasion_depth = np.max(first_cells[:, 1])
            ax.axhline(invasion_depth, color='red', linestyle='--', 
                      linewidth=1.5, alpha=0.8, zorder=5)
        
        ax.set_xlim(0, X_LIM)
        ax.set_ylim(0, Y_LIM)
        if idx == 3:
            ax.set_xlabel('X (μm)', fontsize=FONT_SIZE_LABEL, fontweight='bold')
        if idx == 0:
            ax.set_title(f't=0', fontsize=FONT_SIZE_TITLE, fontweight='bold')
        ax.set_ylabel(f'{titles[ecm_type]}\nY (μm)', fontsize=FONT_SIZE_LABEL, fontweight='bold')
        ax.tick_params(labelsize=FONT_SIZE_TICK)
        ax.grid(True, alpha=0.2)
        
        # Final state (column 1)
        ax = fig.add_subplot(gs[row, 1])
        ax.set_facecolor('#f8f9fa')
        
        # Plot ECM grid from VTI file
        if ecm_timesteps and len(ecm_timesteps) > 0:
            ecm_im = plot_ecm_from_file(ax, ecm_timesteps[-1][1])
        
        # Plot cells
        if len(last_cells) > 0:
            ax.scatter(last_cells[:, 0], last_cells[:, 1], 
                      c=cell_color, s=30, alpha=0.7, edgecolors='darkblue', linewidth=1, zorder=3)
            
            # Add invasion front line (maximum)
            invasion_depth = np.max(last_cells[:, 1])
            ax.axhline(invasion_depth, color='red', linestyle='--', 
                      linewidth=1.5, alpha=0.8, zorder=5)
        
        ax.set_xlim(0, X_LIM)
        ax.set_ylim(0, Y_LIM)
        if idx == 3:
            ax.set_xlabel('X (μm)', fontsize=FONT_SIZE_LABEL, fontweight='bold')
        if idx == 0:
            # final_hours = cell_timesteps[-1][0] / 60
            final_hours = 5*24  # 5 days
            ax.set_title(f't={final_hours:.0f}h', 
                        fontsize=FONT_SIZE_TITLE, fontweight='bold')
        ax.tick_params(labelsize=FONT_SIZE_TICK)
        ax.grid(True, alpha=0.2)
    
    # Add colorbar for ECM density
    if ecm_im is not None:
        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(ecm_im, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=FONT_SIZE_TICK)
        cbar.set_label('ECM Density (volume fraction)', rotation=270, labelpad=20, fontsize=FONT_SIZE_COLORBAR, fontweight='bold')
    
    output_file = output_dir / 'all_ecm_comparison_vertical.pdf'
    plt.savefig(output_file, format='pdf', bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def main():
    """Main function."""
    import argparse
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Export ECM invasion comparison images (4 rows x 2 columns grid)')
    parser.add_argument('--data-dir', type=str, default=None,
                       help='Base data directory (default: ~/Thesis/testoutput/DynamicECMInvasion2d)')
    parser.add_argument('--output-dir', type=str, default=None,
                       help='Output directory for images (default: script_dir/ecm_comparison_output)')
    args = parser.parse_args()
    
    # Set up directories
    script_dir = Path(__file__).parent.resolve()
    thesis_dir = script_dir.parent.parent.parent
    
    if args.data_dir:
        base_dir = Path(args.data_dir).expanduser()
    else:
        base_dir = thesis_dir / "testoutput" / "DynamicECMInvasion2d"
    
    if args.output_dir:
        output_dir = Path(args.output_dir).expanduser()
    else:
        output_dir = script_dir / "ecm_comparison_output"
    
    output_dir.mkdir(exist_ok=True)
    
    print("Creating ECM invasion comparison images (4 rows x 2 columns layout)...")
    print(f"Base directory: {base_dir}")
    print(f"Output directory: {output_dir}")
    
    # Combined comparison
    print("\nCreating combined comparison...")
    create_combined_comparison(base_dir, output_dir)
    
    print(f"\n✓ Complete! Images saved to: {output_dir}")

if __name__ == '__main__':
    main()
