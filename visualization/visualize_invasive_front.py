#!/usr/bin/env python3
"""
visualize_invasive_front.py

Visualize the invasive front simulations with ECM orientation patterns.
Shows cells, ECM fiber orientations, and migration directions.

Usage:
    python visualize_invasive_front.py <scenario>
    
    where <scenario> is one of: random, parallel, perpendicular, mixed
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import LineCollection
import xml.etree.ElementTree as ET
import glob
import os
import sys

def parse_vtu_file(vtu_file):
    """
    Parse VTU file and extract cell positions and data.
    
    Returns:
        positions: Nx2 array of cell positions
        data: dictionary of cell data arrays
    """
    tree = ET.parse(vtu_file)
    root = tree.getroot()
    
    # Find the Points data
    points_data = None
    for piece in root.findall('.//Piece'):
        for points in piece.findall('.//Points'):
            for data_array in points.findall('.//DataArray'):
                if data_array.get('Name') == 'Points':
                    points_data = data_array.text
                    break
    
    if points_data is None:
        return None, None
    
    # Parse points
    points = np.fromstring(points_data, sep=' ').reshape(-1, 3)
    positions = points[:, :2]  # Just x, y for 2D
    
    # Extract cell data
    data = {}
    for piece in root.findall('.//Piece'):
        for point_data in piece.findall('.//PointData'):
            for data_array in point_data.findall('.//DataArray'):
                name = data_array.get('Name')
                if name and data_array.text:
                    values = np.fromstring(data_array.text, sep=' ')
                    data[name] = values
    
    return positions, data


def plot_ecm_field(ax, ecm_type, domain_width=600.0, domain_height=1000.0, grid_spacing=50.0):
    """
    Plot background ECM fiber orientation field.
    
    Args:
        ax: matplotlib axis
        ecm_type: 'random', 'parallel', 'perpendicular', or 'mixed'
        domain_width: width in µm
        domain_height: height in µm
        grid_spacing: spacing between ECM fiber indicators
    """
    # Create grid of points for ECM visualization
    x = np.arange(0, domain_width, grid_spacing)
    y = np.arange(0, domain_height, grid_spacing)
    X, Y = np.meshgrid(x, y)
    
    # Calculate ECM orientation at each point
    if ecm_type == 'random':
        # Random orientations - show as light background pattern
        np.random.seed(42)  # For reproducibility
        angles = np.random.uniform(0, 2*np.pi, X.shape)
        
    elif ecm_type == 'parallel':
        # Horizontal fibers (θ = 0)
        angles = np.zeros(X.shape)
        
    elif ecm_type == 'perpendicular':
        # Vertical fibers (θ = π/2)
        angles = np.full(X.shape, np.pi/2)
        
    elif ecm_type == 'mixed':
        # Perpendicular in center, parallel on sides
        center_start = domain_width / 3.0
        center_end = 2.0 * domain_width / 3.0
        angles = np.where((X > center_start) & (X < center_end), 
                         np.pi/2, 0.0)
    else:
        angles = np.zeros(X.shape)
    
    # Calculate ECM fiber direction vectors
    U = np.cos(angles)
    V = np.sin(angles)
    
    # Plot ECM fibers as light gray lines
    fiber_length = grid_spacing * 0.4
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            x_start = X[i, j] - fiber_length * U[i, j] / 2
            x_end = X[i, j] + fiber_length * U[i, j] / 2
            y_start = Y[i, j] - fiber_length * V[i, j] / 2
            y_end = Y[i, j] + fiber_length * V[i, j] / 2
            ax.plot([x_start, x_end], [y_start, y_end], 
                   color='lightgray', linewidth=1.5, alpha=0.5, zorder=1)


def plot_invasion_snapshot(vtu_file, ecm_type, output_file=None, show_ecm=True, 
                          show_migration=True):
    """
    Create a visualization of the invasive front at a single timepoint.
    
    Args:
        vtu_file: path to VTU file
        ecm_type: ECM orientation type
        output_file: path to save figure (optional)
        show_ecm: whether to show ECM fiber field
        show_migration: whether to show cell migration directions
    """
    positions, data = parse_vtu_file(vtu_file)
    
    if positions is None:
        print(f"Could not parse {vtu_file}")
        return
    
    # Extract timepoint from filename
    timepoint = os.path.basename(vtu_file).replace('results_', '').replace('.vtu', '')
    
    fig, ax = plt.subplots(figsize=(8, 13))
    
    # Plot ECM field in background
    if show_ecm:
        plot_ecm_field(ax, ecm_type)
    
    # Plot cells
    ax.scatter(positions[:, 0], positions[:, 1], c='blue', s=50, 
              alpha=0.6, edgecolors='darkblue', linewidth=1, zorder=3,
              label='Cells')
    
    # Plot ECM orientation vectors at cell positions (if available)
    if 'ecm_orientation_x' in data and 'ecm_orientation_y' in data:
        ecm_x = data['ecm_orientation_x']
        ecm_y = data['ecm_orientation_y']
        
        # Plot as arrows
        scale = 15  # Arrow length scale
        ax.quiver(positions[:, 0], positions[:, 1], 
                 ecm_x, ecm_y,
                 color='green', alpha=0.4, scale=20, width=0.003,
                 label='ECM orientation', zorder=2)
    
    # Plot migration directions (if available and requested)
    if show_migration and 'migration_direction_x' in data and 'migration_direction_y' in data:
        mig_x = data['migration_direction_x']
        mig_y = data['migration_direction_y']
        
        ax.quiver(positions[:, 0], positions[:, 1],
                 mig_x, mig_y,
                 color='red', alpha=0.7, scale=15, width=0.004,
                 label='Migration direction', zorder=4)
    
    # Formatting
    ax.set_xlim(-50, 650)
    ax.set_ylim(-50, 1050)
    ax.set_xlabel('X (µm)', fontsize=12)
    ax.set_ylabel('Y (µm)', fontsize=12)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    
    title = f'Invasive Front - {ecm_type.capitalize()} ECM\nTime: {timepoint} min'
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    ax.legend(loc='upper right', fontsize=10)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_file}")
    else:
        plt.show()
    
    plt.close()


def plot_invasion_comparison(results_dirs, timepoint='1440', output_file=None):
    """
    Create a 2x2 comparison of all four ECM scenarios at the same timepoint.
    
    Args:
        results_dirs: dictionary mapping ecm_type to results directory
        timepoint: timepoint to visualize (default: final)
        output_file: path to save figure
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 20))
    axes = axes.flatten()
    
    ecm_types = ['random', 'parallel', 'perpendicular', 'mixed']
    
    for idx, ecm_type in enumerate(ecm_types):
        ax = axes[idx]
        
        if ecm_type not in results_dirs:
            continue
        
        # Find VTU file for this timepoint
        vtu_pattern = os.path.join(results_dirs[ecm_type], 
                                   f'results_{timepoint}.vtu')
        vtu_files = glob.glob(vtu_pattern)
        
        if not vtu_files:
            print(f"No file found for {ecm_type} at timepoint {timepoint}")
            continue
        
        vtu_file = vtu_files[0]
        positions, data = parse_vtu_file(vtu_file)
        
        if positions is None:
            continue
        
        # Plot ECM field
        plot_ecm_field(ax, ecm_type)
        
        # Plot cells
        ax.scatter(positions[:, 0], positions[:, 1], c='blue', s=30,
                  alpha=0.6, edgecolors='darkblue', linewidth=0.5, zorder=3)
        
        # Calculate invasion depth (95th percentile of y-coordinates)
        invasion_depth = np.percentile(positions[:, 1], 95)
        ax.axhline(invasion_depth, color='red', linestyle='--', linewidth=2,
                  label=f'95th percentile: {invasion_depth:.1f} µm', zorder=5)
        
        # Formatting
        ax.set_xlim(-50, 650)
        ax.set_ylim(-50, 1050)
        ax.set_xlabel('X (µm)', fontsize=11)
        ax.set_ylabel('Y (µm)', fontsize=11)
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        ax.set_title(f'{ecm_type.capitalize()} ECM\n{len(positions)} cells, depth: {invasion_depth:.1f} µm',
                    fontsize=12, fontweight='bold')
        ax.legend(loc='upper right', fontsize=9)
    
    plt.suptitle(f'Invasive Front Comparison - Time: {timepoint} min',
                fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_file}")
    else:
        plt.show()
    
    plt.close()


def plot_invasion_timeseries(results_dir, ecm_type, output_dir=None):
    """
    Create time series analysis of invasion depth.
    
    Args:
        results_dir: path to results directory
        ecm_type: ECM type name
        output_dir: directory to save plots
    """
    # Get all VTU files
    vtu_files = sorted(glob.glob(os.path.join(results_dir, 'results_*.vtu')))
    
    if not vtu_files:
        print(f"No VTU files found in {results_dir}")
        return
    
    timepoints = []
    invasion_depths = []
    cell_counts = []
    mean_y = []
    
    for vtu_file in vtu_files:
        positions, data = parse_vtu_file(vtu_file)
        
        if positions is None:
            continue
        
        # Extract timepoint
        basename = os.path.basename(vtu_file)
        try:
            t = int(basename.replace('results_', '').replace('.vtu', ''))
        except:
            t = 0
        
        timepoints.append(t)
        invasion_depths.append(np.percentile(positions[:, 1], 95))
        cell_counts.append(len(positions))
        mean_y.append(np.mean(positions[:, 1]))
    
    # Convert to arrays and sort by time
    timepoints = np.array(timepoints)
    invasion_depths = np.array(invasion_depths)
    cell_counts = np.array(cell_counts)
    mean_y = np.array(mean_y)
    
    sort_idx = np.argsort(timepoints)
    timepoints = timepoints[sort_idx]
    invasion_depths = invasion_depths[sort_idx]
    cell_counts = cell_counts[sort_idx]
    mean_y = mean_y[sort_idx]
    
    # Create figure
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    
    # Plot invasion depth over time
    ax = axes[0]
    ax.plot(timepoints, invasion_depths, 'b-o', linewidth=2, markersize=4,
           label='95th percentile')
    ax.plot(timepoints, mean_y, 'g--o', linewidth=2, markersize=4,
           label='Mean position')
    ax.set_xlabel('Time (min)', fontsize=12)
    ax.set_ylabel('Y Position (µm)', fontsize=12)
    ax.set_title(f'Invasion Depth - {ecm_type.capitalize()} ECM', 
                fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Plot cell count over time
    ax = axes[1]
    ax.plot(timepoints, cell_counts, 'r-o', linewidth=2, markersize=4)
    ax.set_xlabel('Time (min)', fontsize=12)
    ax.set_ylabel('Cell Count', fontsize=12)
    ax.set_title('Cell Population', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_dir:
        output_file = os.path.join(output_dir, f'invasion_timeseries_{ecm_type}.png')
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_file}")
    else:
        plt.show()
    
    plt.close()


if __name__ == '__main__':
    # Base output directory
    base_dir = '/home/orlando/Thesis/testoutput/InvasiveFront2d'
    
    if len(sys.argv) > 1:
        # Single scenario specified
        ecm_type = sys.argv[1]
        results_dir = os.path.join(base_dir, ecm_type, 'results_from_time_0')
        
        if not os.path.exists(results_dir):
            print(f"Results directory not found: {results_dir}")
            sys.exit(1)
        
        # Plot final snapshot
        vtu_file = os.path.join(results_dir, 'results_1440.vtu')
        if os.path.exists(vtu_file):
            plot_invasion_snapshot(vtu_file, ecm_type, 
                                 output_file=f'invasion_{ecm_type}_final.png')
        
        # Plot time series
        plot_invasion_timeseries(results_dir, ecm_type, output_dir='.')
        
    else:
        # Compare all scenarios
        results_dirs = {}
        for ecm_type in ['random', 'parallel', 'perpendicular', 'mixed']:
            results_dir = os.path.join(base_dir, ecm_type, 'results_from_time_0')
            if os.path.exists(results_dir):
                results_dirs[ecm_type] = results_dir
        
        if results_dirs:
            # Create comparison plot
            plot_invasion_comparison(results_dirs, timepoint='1440',
                                   output_file='invasion_comparison_all.png')
            
            # Create time series for each scenario
            for ecm_type, results_dir in results_dirs.items():
                plot_invasion_timeseries(results_dir, ecm_type, output_dir='.')
        else:
            print("No results found. Please run the simulations first.")
