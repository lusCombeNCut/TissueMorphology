#!/usr/bin/env python3
"""
visualize_invasive_front_simple.py

Simple visualization of invasive front using .viznodes files (uncompressed text format).
This avoids dealing with compressed VTU files.

Usage:
    python visualize_invasive_front_simple.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import glob

def read_viznodes_file(viznodes_file):
    """
    Read .viznodes file (Chaste format: timepoint followed by x y pairs on same line).
    
    Returns:
        timepoints: array of timepoint values
        positions_dict: dictionary mapping timepoint to Nx3 array of positions
    """
    with open(viznodes_file, 'r') as f:
        lines = f.readlines()
    
    positions_dict = {}
    
    for line in lines:
        parts = line.strip().split()
        if not parts or len(parts) < 3:
            continue
        
        # First value is timepoint, rest are x y pairs
        timepoint = float(parts[0])
        coords = list(map(float, parts[1:]))
        
        # Group into x, y pairs (z=0 for 2D)
        positions = []
        for i in range(0, len(coords), 2):
            if i+1 < len(coords):
                x = coords[i]
                y = coords[i+1]
                positions.append([x, y, 0.0])
        
        if positions:
            positions_dict[timepoint] = np.array(positions)
    
    timepoints = sorted(positions_dict.keys())
    return timepoints, positions_dict


def plot_ecm_field(ax, ecm_type, domain_width=600.0, domain_height=1000.0, grid_spacing=50.0):
    """Plot background ECM fiber orientation field."""
    x = np.arange(0, domain_width, grid_spacing)
    y = np.arange(0, domain_height, grid_spacing)
    X, Y = np.meshgrid(x, y)
    
    if ecm_type == 'random':
        np.random.seed(42)
        angles = np.random.uniform(0, 2*np.pi, X.shape)
    elif ecm_type == 'parallel':
        angles = np.zeros(X.shape)
    elif ecm_type == 'perpendicular':
        angles = np.full(X.shape, np.pi/2)
    elif ecm_type == 'mixed':
        center_start = domain_width / 3.0
        center_end = 2.0 * domain_width / 3.0
        angles = np.where((X > center_start) & (X < center_end), np.pi/2, 0.0)
    else:
        angles = np.zeros(X.shape)
    
    U = np.cos(angles)
    V = np.sin(angles)
    
    fiber_length = grid_spacing * 0.4
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            x_start = X[i, j] - fiber_length * U[i, j] / 2
            x_end = X[i, j] + fiber_length * U[i, j] / 2
            y_start = Y[i, j] - fiber_length * V[i, j] / 2
            y_end = Y[i, j] + fiber_length * V[i, j] / 2
            ax.plot([x_start, x_end], [y_start, y_end], 
                   color='lightgray', linewidth=1.5, alpha=0.5, zorder=1)


def plot_invasion_snapshot(positions, ecm_type, timepoint, output_file=None):
    """Create visualization at a single timepoint."""
    fig, ax = plt.subplots(figsize=(8, 13))
    
    # Plot ECM field
    plot_ecm_field(ax, ecm_type)
    
    # Plot cells
    ax.scatter(positions[:, 0], positions[:, 1], c='blue', s=50, 
              alpha=0.7, edgecolors='darkblue', linewidth=1, zorder=3,
              label=f'Cells (n={len(positions)})')
    
    # Calculate invasion depth
    invasion_depth = np.percentile(positions[:, 1], 95)
    ax.axhline(invasion_depth, color='red', linestyle='--', linewidth=2,
              label=f'95th percentile: {invasion_depth:.1f} µm', zorder=5)
    
    # Formatting
    ax.set_xlim(-50, 650)
    ax.set_ylim(-50, 1050)
    ax.set_xlabel('X (µm)', fontsize=12)
    ax.set_ylabel('Y (µm)', fontsize=12)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    
    title = f'Invasive Front - {ecm_type.capitalize()} ECM\nTime: {timepoint:.0f} min'
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_file}")
    else:
        plt.show()
    
    plt.close()


def plot_invasion_comparison(base_dir, output_file=None):
    """Create 2x2 comparison of all four ECM scenarios."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 20))
    axes = axes.flatten()
    
    ecm_types = ['random', 'parallel', 'perpendicular', 'mixed']
    
    for idx, ecm_type in enumerate(ecm_types):
        ax = axes[idx]
        
        viznodes_file = os.path.join(base_dir, ecm_type, 
                                    'results_from_time_0', 'results.viznodes')
        
        if not os.path.exists(viznodes_file):
            print(f"File not found: {viznodes_file}")
            continue
        
        timepoints, positions_dict = read_viznodes_file(viznodes_file)
        
        if not timepoints:
            continue
        
        # Use final timepoint
        final_time = timepoints[-1]
        positions = positions_dict[final_time]
        
        # Plot ECM field
        plot_ecm_field(ax, ecm_type)
        
        # Plot cells
        ax.scatter(positions[:, 0], positions[:, 1], c='blue', s=30,
                  alpha=0.6, edgecolors='darkblue', linewidth=0.5, zorder=3)
        
        # Calculate invasion depth
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
    
    plt.suptitle(f'Invasive Front Comparison - Final Timepoint',
                fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_file}")
    else:
        plt.show()
    
    plt.close()


def plot_invasion_timeseries(viznodes_file, ecm_type, output_file=None):
    """Create time series analysis."""
    timepoints, positions_dict = read_viznodes_file(viznodes_file)
    
    if not timepoints:
        print(f"No data found in {viznodes_file}")
        return
    
    timepoints = np.array(timepoints)
    invasion_depths = []
    cell_counts = []
    mean_y = []
    
    for t in timepoints:
        positions = positions_dict[t]
        invasion_depths.append(np.percentile(positions[:, 1], 95))
        cell_counts.append(len(positions))
        mean_y.append(np.mean(positions[:, 1]))
    
    invasion_depths = np.array(invasion_depths)
    cell_counts = np.array(cell_counts)
    mean_y = np.array(mean_y)
    
    # Create figure
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    
    # Plot invasion depth
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
    
    # Plot cell count
    ax = axes[1]
    ax.plot(timepoints, cell_counts, 'r-o', linewidth=2, markersize=4)
    ax.set_xlabel('Time (min)', fontsize=12)
    ax.set_ylabel('Cell Count', fontsize=12)
    ax.set_title('Cell Population', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_file}")
    else:
        plt.show()
    
    plt.close()


def create_invasion_animation_frames(viznodes_file, ecm_type, output_dir):
    """Create animation frames for each timepoint."""
    timepoints, positions_dict = read_viznodes_file(viznodes_file)
    
    if not timepoints:
        return
    
    os.makedirs(output_dir, exist_ok=True)
    
    for i, t in enumerate(timepoints):
        positions = positions_dict[t]
        output_file = os.path.join(output_dir, f'frame_{i:04d}.png')
        plot_invasion_snapshot(positions, ecm_type, t, output_file)


if __name__ == '__main__':
    base_dir = '/home/orlando/Thesis/testoutput/InvasiveFront2d'
    
    # Create comparison plot
    plot_invasion_comparison(base_dir, output_file='invasion_comparison.png')
    
    # Create time series for each scenario
    for ecm_type in ['random', 'parallel', 'perpendicular', 'mixed']:
        viznodes_file = os.path.join(base_dir, ecm_type, 
                                     'results_from_time_0', 'results.viznodes')
        if os.path.exists(viznodes_file):
            # Time series plot
            plot_invasion_timeseries(viznodes_file, ecm_type,
                                    output_file=f'invasion_timeseries_{ecm_type}.png')
            
            # Final snapshot with ECM field
            timepoints, positions_dict = read_viznodes_file(viznodes_file)
            if timepoints:
                final_time = timepoints[-1]
                plot_invasion_snapshot(positions_dict[final_time], ecm_type, 
                                      final_time, 
                                      output_file=f'invasion_{ecm_type}_final.png')
    
    print("\nVisualization complete!")
    print("Generated files:")
    print("  - invasion_comparison.png (4-panel comparison)")
    print("  - invasion_timeseries_*.png (time series for each scenario)")
    print("  - invasion_*_final.png (final snapshot with ECM field)")
