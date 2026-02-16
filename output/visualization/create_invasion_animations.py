#!/usr/bin/env python3
"""
create_invasion_animations.py

Create GIF animations showing invasion front dynamics over time.
Also creates ECM grid visualization for better understanding of the fiber field.

Usage:
    python create_invasion_animations.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation, PillowWriter
import os
import sys

def read_viznodes_file(viznodes_file):
    """Read Chaste .viznodes file format."""
    with open(viznodes_file, 'r') as f:
        lines = f.readlines()
    
    positions_dict = {}
    
    for line in lines:
        parts = line.strip().split()
        if not parts or len(parts) < 3:
            continue
        
        timepoint = float(parts[0])
        coords = list(map(float, parts[1:]))
        
        positions = []
        for i in range(0, len(coords), 2):
            if i+1 < len(coords):
                positions.append([coords[i], coords[i+1], 0.0])
        
        if positions:
            positions_dict[timepoint] = np.array(positions)
    
    return sorted(positions_dict.keys()), positions_dict


def plot_ecm_grid(ax, ecm_type, domain_width=600.0, domain_height=1000.0, grid_spacing=50.0):
    """Plot ECM fiber field on a regular grid."""
    x = np.arange(0, domain_width + grid_spacing, grid_spacing)
    y = np.arange(0, domain_height + grid_spacing, grid_spacing)
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
    
    # Plot as quiver (arrow glyphs)
    ax.quiver(X, Y, U, V, 
             color='lightgray', alpha=0.6, 
             scale=25, width=0.003,
             headwidth=4, headlength=5,
             zorder=1)
    
    return X, Y, U, V


def create_frame(ax, positions, ecm_type, timepoint, show_ecm_grid=True):
    """Create a single frame for animation."""
    ax.clear()
    
    # Plot ECM grid
    if show_ecm_grid:
        plot_ecm_grid(ax, ecm_type, grid_spacing=50.0)
    
    # Plot cells
    if len(positions) > 0:
        ax.scatter(positions[:, 0], positions[:, 1], 
                  c='blue', s=50, alpha=0.7, 
                  edgecolors='darkblue', linewidth=1, 
                  zorder=3, label=f'Cells (n={len(positions)})')
        
        # Calculate and show invasion depth
        invasion_depth = np.percentile(positions[:, 1], 95)
        ax.axhline(invasion_depth, color='red', linestyle='--', 
                  linewidth=2, alpha=0.8,
                  label=f'Front: {invasion_depth:.0f} µm', zorder=5)
    
    # Formatting
    ax.set_xlim(-50, 650)
    ax.set_ylim(-50, 1050)
    ax.set_xlabel('X (µm)', fontsize=11)
    ax.set_ylabel('Y (µm)', fontsize=11)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_title(f'{ecm_type.capitalize()} ECM - Time: {timepoint:.0f} min',
                fontsize=12, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)


def create_animation_gif(viznodes_file, ecm_type, output_file, 
                         frame_skip=1, show_ecm_grid=True):
    """
    Create animated GIF of invasion over time.
    
    Args:
        viznodes_file: path to .viznodes file
        ecm_type: ECM scenario name
        output_file: output GIF filename
        frame_skip: use every Nth frame (1=all frames, 2=every other, etc.)
        show_ecm_grid: whether to show ECM grid
    """
    print(f"Creating animation for {ecm_type}...")
    
    timepoints, positions_dict = read_viznodes_file(viznodes_file)
    
    if not timepoints:
        print(f"  No data found in {viznodes_file}")
        return
    
    # Subsample frames if requested
    timepoints = timepoints[::frame_skip]
    
    print(f"  Using {len(timepoints)} frames")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 13))
    
    def update(frame_idx):
        t = timepoints[frame_idx]
        positions = positions_dict[t]
        create_frame(ax, positions, ecm_type, t, show_ecm_grid)
        return ax,
    
    # Create animation
    anim = FuncAnimation(fig, update, frames=len(timepoints),
                        interval=200, blit=False, repeat=True)
    
    # Save as GIF
    writer = PillowWriter(fps=5)
    anim.save(output_file, writer=writer, dpi=100)
    print(f"  Saved: {output_file}")
    
    plt.close()


def create_comparison_gif(base_dir, output_file, frame_skip=1):
    """Create 2x2 comparison GIF of all scenarios."""
    print("Creating comparison animation...")
    
    ecm_types = ['random', 'parallel', 'perpendicular', 'mixed']
    
    # Load all data
    data = {}
    all_timepoints = None
    
    for ecm_type in ecm_types:
        viznodes_file = os.path.join(base_dir, ecm_type, 
                                     'results_from_time_0', 'results.viznodes')
        if os.path.exists(viznodes_file):
            timepoints, positions_dict = read_viznodes_file(viznodes_file)
            data[ecm_type] = (timepoints, positions_dict)
            if all_timepoints is None:
                all_timepoints = timepoints
    
    if not data:
        print("  No data found")
        return
    
    all_timepoints = all_timepoints[::frame_skip]
    print(f"  Using {len(all_timepoints)} frames")
    
    # Create 2x2 figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 18))
    axes = axes.flatten()
    
    def update(frame_idx):
        t = all_timepoints[frame_idx]
        
        for idx, ecm_type in enumerate(ecm_types):
            ax = axes[idx]
            ax.clear()
            
            if ecm_type not in data:
                continue
            
            timepoints, positions_dict = data[ecm_type]
            
            # Find closest timepoint
            closest_t = min(timepoints, key=lambda x: abs(x - t))
            positions = positions_dict[closest_t]
            
            # Plot ECM grid
            plot_ecm_grid(ax, ecm_type, grid_spacing=60.0)
            
            # Plot cells
            ax.scatter(positions[:, 0], positions[:, 1], 
                      c='blue', s=25, alpha=0.6,
                      edgecolors='darkblue', linewidth=0.5, zorder=3)
            
            # Invasion depth
            invasion_depth = np.percentile(positions[:, 1], 95)
            ax.axhline(invasion_depth, color='red', linestyle='--',
                      linewidth=2, alpha=0.8, zorder=5)
            
            # Formatting
            ax.set_xlim(-50, 650)
            ax.set_ylim(-50, 1050)
            ax.set_xlabel('X (µm)', fontsize=10)
            ax.set_ylabel('Y (µm)', fontsize=10)
            ax.set_aspect('equal')
            ax.grid(True, alpha=0.3)
            ax.set_title(f'{ecm_type.capitalize()}\nDepth: {invasion_depth:.0f} µm',
                        fontsize=11, fontweight='bold')
        
        fig.suptitle(f'Invasive Front Comparison - Time: {t:.0f} min',
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        return axes
    
    # Create animation
    anim = FuncAnimation(fig, update, frames=len(all_timepoints),
                        interval=200, blit=False, repeat=True)
    
    # Save as GIF
    writer = PillowWriter(fps=5)
    anim.save(output_file, writer=writer, dpi=80)
    print(f"  Saved: {output_file}")
    
    plt.close()


def create_ecm_grid_only(ecm_type, output_file):
    """Create static image showing just the ECM grid pattern."""
    fig, ax = plt.subplots(figsize=(8, 13))
    
    # Plot ECM with finer grid for clarity
    X, Y, U, V = plot_ecm_grid(ax, ecm_type, grid_spacing=40.0)
    
    # Add colored background regions for mixed
    if ecm_type == 'mixed':
        center_start = 600.0 / 3.0
        center_end = 2.0 * 600.0 / 3.0
        ax.axvspan(center_start, center_end, alpha=0.1, color='green',
                  label='Perpendicular region')
        ax.axvspan(-50, center_start, alpha=0.1, color='red', label='Parallel region')
        ax.axvspan(center_end, 650, alpha=0.1, color='red')
    
    ax.set_xlim(-50, 650)
    ax.set_ylim(-50, 1050)
    ax.set_xlabel('X (µm)', fontsize=12)
    ax.set_ylabel('Y (µm)', fontsize=12)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_title(f'ECM Fiber Field - {ecm_type.capitalize()}',
                fontsize=14, fontweight='bold')
    if ecm_type == 'mixed':
        ax.legend(loc='upper right')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved ECM grid: {output_file}")
    plt.close()


if __name__ == '__main__':
    base_dir = '/home/orlando/Thesis/testoutput/InvasiveFront2d'
    
    # Check if required packages are available
    try:
        from matplotlib.animation import PillowWriter
    except ImportError:
        print("Error: Pillow is required for GIF creation")
        print("Install with: pip install Pillow")
        sys.exit(1)
    
    print("Creating invasion front animations...")
    print("This may take a few minutes...\n")
    
    # Create individual scenario animations
    for ecm_type in ['random', 'parallel', 'perpendicular', 'mixed']:
        viznodes_file = os.path.join(base_dir, ecm_type,
                                     'results_from_time_0', 'results.viznodes')
        
        if os.path.exists(viznodes_file):
            # Animation with ECM grid
            output_file = f'invasion_{ecm_type}_animation.gif'
            create_animation_gif(viznodes_file, ecm_type, output_file,
                               frame_skip=1, show_ecm_grid=True)
            
            # Static ECM grid reference
            ecm_grid_file = f'ecm_grid_{ecm_type}.png'
            create_ecm_grid_only(ecm_type, ecm_grid_file)
    
    # Create comparison animation (all 4 scenarios)
    print("\nCreating comparison animation...")
    create_comparison_gif(base_dir, 'invasion_comparison_animation.gif',
                         frame_skip=1)
    
    print("\n✓ Animation creation complete!")
    print("\nGenerated files:")
    print("  - invasion_*_animation.gif (individual scenarios)")
    print("  - invasion_comparison_animation.gif (4-panel comparison)")
    print("  - ecm_grid_*.png (static ECM field patterns)")
