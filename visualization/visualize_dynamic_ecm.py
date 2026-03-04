#!/usr/bin/env python3
"""
Visualize Dynamic ECM Invasion Simulations

Shows:
1. Cell positions and migration
2. ECM fiber orientation (from cell data)
3. ECM density heatmap
4. ECM remodeling over time
5. Invasion depth progression

Compares static vs dynamic ECM
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import os
import glob
import vtk
from scipy.interpolate import griddata

def parse_viznodes_file(filepath):
    """Parse Chaste .viznodes file format"""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    timepoints = {}
    for line in lines:
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) < 3:
            continue
        
        time = float(parts[0])
        coords = [float(x) for x in parts[1:]]
        
        # Coords are: x1 y1 x2 y2 x3 y3 ...
        num_cells = len(coords) // 2
        cells = np.array(coords).reshape(num_cells, 2)
        
        timepoints[time] = cells
    
    return timepoints

def read_vtu_file(filepath):
    """Read VTU file and extract cell positions and ECM data"""
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(filepath)
    reader.Update()
    output = reader.GetOutput()
    
    # Get cell positions
    points = output.GetPoints()
    num_points = points.GetNumberOfPoints()
    
    positions = np.zeros((num_points, 2))
    for i in range(num_points):
        pos = points.GetPoint(i)
        positions[i] = [pos[0], pos[1]]
    
    # Get cell data
    point_data = output.GetPointData()
    
    cell_data = {}
    
    # Get is_apical if available
    if point_data.HasArray('is_apical'):
        cell_data['is_apical'] = np.array([point_data.GetArray('is_apical').GetValue(i) 
                                            for i in range(num_points)])
    
    # Get ECM data
    ecm_data = {}
    if point_data.HasArray('ecm_angle'):
        ecm_data['angle'] = np.array([point_data.GetArray('ecm_angle').GetValue(i) 
                                       for i in range(num_points)])
    if point_data.HasArray('ecm_density'):
        ecm_data['density'] = np.array([point_data.GetArray('ecm_density').GetValue(i) 
                                         for i in range(num_points)])
    if point_data.HasArray('ecm_orientation_x'):
        ecm_data['orientation_x'] = np.array([point_data.GetArray('ecm_orientation_x').GetValue(i) 
                                               for i in range(num_points)])
    if point_data.HasArray('ecm_orientation_y'):
        ecm_data['orientation_y'] = np.array([point_data.GetArray('ecm_orientation_y').GetValue(i) 
                                               for i in range(num_points)])
    
    return positions, ecm_data, cell_data

def read_ecm_grid_vtk(filepath):
    """Read ECM grid VTK file (structured points)"""
    reader = vtk.vtkStructuredPointsReader()
    reader.SetFileName(filepath)
    reader.Update()
    output = reader.GetOutput()
    
    # Get dimensions and spacing
    dims = output.GetDimensions()
    origin = output.GetOrigin()
    spacing = output.GetSpacing()
    
    # Create coordinate arrays
    nx, ny = dims[0], dims[1]
    x = origin[0] + np.arange(nx) * spacing[0]
    y = origin[1] + np.arange(ny) * spacing[1]
    X, Y = np.meshgrid(x, y, indexing='ij')
    
    # Get ECM data
    point_data = output.GetPointData()
    
    ecm_grid = {}
    
    if point_data.HasArray('ecm_density'):
        density_array = point_data.GetArray('ecm_density')
        density = np.array([density_array.GetValue(i) for i in range(density_array.GetNumberOfTuples())])
        ecm_grid['density'] = density.reshape((nx, ny), order='F')  # Fortran order for VTK
    
    if point_data.HasArray('ecm_anisotropy'):
        anisotropy_array = point_data.GetArray('ecm_anisotropy')
        anisotropy = np.array([anisotropy_array.GetValue(i) for i in range(anisotropy_array.GetNumberOfTuples())])
        ecm_grid['anisotropy'] = anisotropy.reshape((nx, ny), order='F')
    
    if point_data.HasArray('ecm_orientation'):
        orientation_array = point_data.GetArray('ecm_orientation')
        orientation_x = np.array([orientation_array.GetTuple3(i)[0] for i in range(orientation_array.GetNumberOfTuples())])
        orientation_y = np.array([orientation_array.GetTuple3(i)[1] for i in range(orientation_array.GetNumberOfTuples())])
        ecm_grid['orientation_x'] = orientation_x.reshape((nx, ny), order='F')
        ecm_grid['orientation_y'] = orientation_y.reshape((nx, ny), order='F')
    
    ecm_grid['X'] = X
    ecm_grid['Y'] = Y
    
    return ecm_grid

def visualize_comparison(static_dir, dynamic_dir, scenario_name, output_dir):
    """Compare static vs dynamic ECM for a scenario"""
    
    # Parse viznodes for static
    static_viznodes = os.path.join(static_dir, f"{scenario_name}/results_from_time_0/results.viznodes")
    # Use VTU for dynamic cells
    dynamic_vtu = os.path.join(dynamic_dir, f"{scenario_name}/results_from_time_0/results_1440.vtu")
    # Use ECM grid VTK for ECM field
    ecm_grid_vtk = os.path.join(dynamic_dir, f"{scenario_name}/results_from_time_0/ecm_grid_1440.vtk")
    
    if not os.path.exists(static_viznodes):
        print(f"Static file not found: {static_viznodes}")
        return
    if not os.path.exists(dynamic_vtu):
        print(f"Dynamic file not found: {dynamic_vtu}")
        return
    if not os.path.exists(ecm_grid_vtk):
        print(f"ECM grid file not found: {ecm_grid_vtk}")
        return
    
    # Parse static data
    static_data = parse_viznodes_file(static_viznodes)
    static_times = sorted(static_data.keys())
    final_time_static = static_times[-1]
    cells_static = static_data[final_time_static]
    
    # Parse dynamic data
    cells_dynamic, _, cell_data = read_vtu_file(dynamic_vtu)
    
    # Read ECM grid
    ecm_grid = read_ecm_grid_vtk(ecm_grid_vtk)
    
    # Create comparison plot
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    
    # Static ECM
    ax = axes[0]
    ax.scatter(cells_static[:, 0], cells_static[:, 1], c='blue', s=50, alpha=0.6, label='Cells')
    ax.set_xlim(0, 600)
    ax.set_ylim(0, 1000)
    ax.set_xlabel('X (µm)')
    ax.set_ylabel('Y (µm)')
    ax.set_title(f'Static ECM - {scenario_name}\n{len(cells_static)} cells at t={final_time_static:.0f} min')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    
    static_depth = np.max(cells_static[:, 1])
    ax.axhline(static_depth, color='red', linestyle='--', label=f'Invasion depth: {static_depth:.1f} µm')
    ax.legend()
    
    # Dynamic ECM with density heatmap and fiber arrows FROM GRID
    ax = axes[1]
    
    # Plot ECM density heatmap
    if 'density' in ecm_grid:
        im = ax.imshow(ecm_grid['density'].T, extent=[0, 600, 0, 1000], origin='lower', 
                       cmap='YlOrRd', alpha=0.6, vmin=0, vmax=1.5)
        plt.colorbar(im, ax=ax, label='ECM Density')
    
    # Plot ECM fiber arrows from grid
    if 'orientation_x' in ecm_grid and 'orientation_y' in ecm_grid:
        # Subsample for clarity
        skip = 3
        ax.quiver(ecm_grid['X'][::skip, ::skip], ecm_grid['Y'][::skip, ::skip],
                 ecm_grid['orientation_x'][::skip, ::skip], ecm_grid['orientation_y'][::skip, ::skip],
                 color='black', alpha=0.5, scale=25, width=0.003, headwidth=3, headlength=4,
                 label='ECM fibers')
    
    # Plot cells on top - color by cell type
    if 'is_apical' in cell_data and len(cell_data['is_apical']) > 0:
        # Color by apical (red) vs basal (blue)
        colors = ['red' if is_apical > 0.5 else 'blue' for is_apical in cell_data['is_apical']]
        ax.scatter(cells_dynamic[:, 0], cells_dynamic[:, 1], c=colors, s=50, alpha=0.8, 
                   edgecolors='black', linewidth=0.5, label='Cells (red=apical, blue=basal)', zorder=10)
    else:
        # Default green if no cell type data
        ax.scatter(cells_dynamic[:, 0], cells_dynamic[:, 1], c='green', s=50, alpha=0.8, 
                   edgecolors='darkgreen', linewidth=0.5, label='Cells', zorder=10)
    
    ax.set_xlim(0, 600)
    ax.set_ylim(0, 1000)
    ax.set_xlabel('X (µm)')
    ax.set_ylabel('Y (µm)')
    ax.set_title(f'Dynamic ECM - {scenario_name}\n{len(cells_dynamic)} cells at t=1440 min')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    
    dynamic_depth = np.max(cells_dynamic[:, 1])
    ax.axhline(dynamic_depth, color='red', linestyle='--', label=f'Invasion depth: {dynamic_depth:.1f} µm')
    ax.legend(loc='upper right')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, f'dynamic_vs_static_{scenario_name}.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()
    
    # Print statistics
    print(f"\n{scenario_name.upper()} COMPARISON:")
    print(f"  Static ECM:  {len(cells_static)} cells, depth {static_depth:.1f} µm")
    print(f"  Dynamic ECM: {len(cells_dynamic)} cells, depth {dynamic_depth:.1f} µm")
    if 'density' in ecm_grid:
        mean_density = np.mean(ecm_grid['density'])
        print(f"  Mean ECM density: {mean_density:.3f}")
    print(f"  Difference:  {dynamic_depth - static_depth:+.1f} µm ({100*(dynamic_depth-static_depth)/static_depth:+.1f}%)")

def create_dynamic_animation(dynamic_dir, scenario_name, output_dir):
    """Create animation showing dynamic ECM evolution"""
    
    results_dir = os.path.join(dynamic_dir, f"{scenario_name}/results_from_time_0")
    
    if not os.path.exists(results_dir):
        print(f"Directory not found: {results_dir}")
        return
    
    # Get all VTU files (for cells) and ECM grid files
    vtu_files_unsorted = glob.glob(os.path.join(results_dir, "results_*.vtu"))
    ecm_grid_files_unsorted = glob.glob(os.path.join(results_dir, "ecm_grid_*.vtk"))
    
    if len(vtu_files_unsorted) == 0:
        print(f"No VTU files found in {results_dir}")
        return
    if len(ecm_grid_files_unsorted) == 0:
        print(f"No ECM grid files found in {results_dir}")
        return
    
    # Extract times from VTU filenames and sort by time
    vtu_time_pairs = []
    for vtu_file in vtu_files_unsorted:
        basename = os.path.basename(vtu_file)
        time_str = basename.replace('results_', '').replace('.vtu', '')
        try:
            time = float(time_str)
            vtu_time_pairs.append((time, vtu_file))
        except:
            pass
    
    # Sort by time
    vtu_time_pairs.sort(key=lambda x: x[0])
    times = [t for t, _ in vtu_time_pairs]
    vtu_files = [f for _, f in vtu_time_pairs]
    
    # Sort ECM grid files by time as well
    ecm_time_pairs = []
    for ecm_file in ecm_grid_files_unsorted:
        basename = os.path.basename(ecm_file)
        time_str = basename.replace('ecm_grid_', '').replace('.vtk', '')
        try:
            time = float(time_str)
            ecm_time_pairs.append((time, ecm_file))
        except:
            pass
    ecm_time_pairs.sort(key=lambda x: x[0])
    ecm_grid_files = [f for _, f in ecm_time_pairs]
    
    # Create figure with space for colorbar
    fig, ax = plt.subplots(figsize=(10, 11))
    
    # Reserve space for colorbar
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    
    def update(frame):
        ax.clear()
        cax.clear()
        time = times[frame]
        vtu_file = vtu_files[frame]
        
        # Find corresponding ECM grid file - try exact match first, then closest
        ecm_grid_file = None
        
        # Try exact match
        target_time = int(time)
        for ecm_file in ecm_grid_files:
            if f"ecm_grid_{target_time}.vtk" in ecm_file:
                ecm_grid_file = ecm_file
                break
        
        # If no exact match, find closest time
        if ecm_grid_file is None and len(ecm_grid_files) > 0:
            ecm_times = []
            for ecm_file in ecm_grid_files:
                basename = os.path.basename(ecm_file)
                ecm_time_str = basename.replace('ecm_grid_', '').replace('.vtk', '')
                try:
                    ecm_times.append((float(ecm_time_str), ecm_file))
                except:
                    pass
            
            if ecm_times:
                # Find closest time
                closest = min(ecm_times, key=lambda x: abs(x[0] - time))
                if abs(closest[0] - time) < 61:  # Within 61 minutes (allow one timestep mismatch)
                    ecm_grid_file = closest[1]
        
        # Read cell data
        cells, _, cell_data = read_vtu_file(vtu_file)
        
        # Read ECM grid
        mean_density = 0
        if ecm_grid_file and os.path.exists(ecm_grid_file):
            try:
                ecm_grid = read_ecm_grid_vtk(ecm_grid_file)
                
                # Plot ECM density heatmap
                if 'density' in ecm_grid:
                    im = ax.imshow(ecm_grid['density'].T, extent=[0, 600, 0, 1000], origin='lower', 
                                   cmap='YlOrRd', alpha=0.5, vmin=0, vmax=1.5)
                    # Add colorbar
                    plt.colorbar(im, cax=cax, label='ECM Density')
                
                # Plot ECM fiber arrows from grid
                if 'orientation_x' in ecm_grid and 'orientation_y' in ecm_grid:
                    skip = 3
                    ax.quiver(ecm_grid['X'][::skip, ::skip], ecm_grid['Y'][::skip, ::skip],
                             ecm_grid['orientation_x'][::skip, ::skip], ecm_grid['orientation_y'][::skip, ::skip],
                             color='black', alpha=0.6, scale=25, width=0.004, headwidth=3, headlength=4)
                
                mean_density = np.mean(ecm_grid['density']) if 'density' in ecm_grid else 0
            except Exception as e:
                print(f"Warning: Could not read ECM grid file {ecm_grid_file}: {e}")
                mean_density = 0
        else:
            mean_density = 0
        
        # Plot cells on top - color by cell type
        if 'is_apical' in cell_data and len(cell_data['is_apical']) > 0:
            # Color by apical (red) vs basal (blue)
            colors = ['red' if is_apical > 0.5 else 'blue' for is_apical in cell_data['is_apical']]
            ax.scatter(cells[:, 0], cells[:, 1], c=colors, s=60, alpha=0.8, 
                       edgecolors='black', linewidth=0.8, zorder=10)
        else:
            # Default green if no cell type data
            ax.scatter(cells[:, 0], cells[:, 1], c='green', s=60, alpha=0.8, 
                       edgecolors='darkgreen', linewidth=0.8, zorder=10)
        
        # Plot settings
        ax.set_xlim(0, 600)
        ax.set_ylim(0, 1000)
        ax.set_xlabel('X (µm)', fontsize=12)
        ax.set_ylabel('Y (µm)', fontsize=12)
        
        title = f'Dynamic ECM - {scenario_name}\nTime: {time:.0f} min, Cells: {len(cells)}'
        title += f', Mean Density: {mean_density:.2f}'
        ax.set_title(title, fontsize=13)
        
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        
        # Show invasion depth
        depth = np.max(cells[:, 1])
        ax.axhline(depth, color='red', linestyle='--', alpha=0.6, linewidth=2,
                   label=f'Depth: {depth:.1f} µm')
        ax.legend(loc='upper right', fontsize=10)
    
    anim = FuncAnimation(fig, update, frames=len(vtu_files), interval=200, repeat=True)
    
    output_file = os.path.join(output_dir, f'dynamic_ecm_{scenario_name}_animation.gif')
    writer = PillowWriter(fps=5)
    anim.save(output_file, writer=writer)
    print(f"✓ Saved: {output_file}")
    plt.close()

def main():
    """Main visualization routine"""
    
    testoutput = "/home/orlando/Thesis/testoutput"
    static_dir = os.path.join(testoutput, "InvasiveFront2d")
    dynamic_dir = os.path.join(testoutput, "DynamicECMInvasion2d")
    
    # Create output directory for visualizations
    output_dir = "/home/orlando/Thesis/Chaste/projects/TissueMorphology/visualization_output"
    os.makedirs(output_dir, exist_ok=True)
    
    scenarios = ["random", "parallel", "perpendicular", "mixed"]
    
    print("=" * 60)
    print("DYNAMIC ECM VISUALIZATION")
    print("=" * 60)
    print(f"Output directory: {output_dir}")
    
    # Create comparisons
    print("\n1. Creating static vs dynamic comparisons...")
    for scenario in scenarios:
        visualize_comparison(static_dir, dynamic_dir, scenario, output_dir)
    
    # Create animations
    print("\n2. Creating dynamic ECM animations...")
    for scenario in scenarios:
        create_dynamic_animation(dynamic_dir, scenario, output_dir)
    
    print("\n" + "=" * 60)
    print("✓ Visualization complete!")
    print("=" * 60)
    print(f"\nGenerated files in {output_dir}:")
    print("  - dynamic_vs_static_*.png (4 files)")
    print("  - dynamic_ecm_*_animation.gif (4 files)")

if __name__ == "__main__":
    main()
