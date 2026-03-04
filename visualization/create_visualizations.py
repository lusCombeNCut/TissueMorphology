#!/usr/bin/env python3
"""
Create comparison visualizations for different organoid scenarios.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from visualize_chaste_results import parse_chaste_data

def create_comparison_figure():
    """
    Create a figure comparing different organoid formation scenarios.
    """
    
    # Define the scenarios and their data files
    scenarios = {
        'Basic Test': '/home/orlando/Thesis/testoutput/OrganoidFormation/BasicTest/results_from_time_0/cellages.dat',
        'Low Stiffness': '/home/orlando/Thesis/testoutput/OrganoidFormation/LowStiffness/results_from_time_0/cellages.dat',
        'High Stiffness': '/home/orlando/Thesis/testoutput/OrganoidFormation/HighStiffness/results_from_time_0.2/cellages.dat',
        'With Basement Membrane': '/home/orlando/Thesis/testoutput/OrganoidFormation/WithBasementMembrane/results_from_time_0/cellages.dat'
    }
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    axes = axes.flatten()
    
    for i, (scenario_name, data_file) in enumerate(scenarios.items()):
        if not os.path.exists(data_file):
            print(f"Warning: Data file not found for {scenario_name}: {data_file}")
            continue
            
        data = parse_chaste_data(data_file)
        if not data:
            print(f"Warning: No data parsed for {scenario_name}")
            continue
        
        # Get the last time point for comparison
        time_points = sorted(data.keys())
        last_time = time_points[-1]
        cells = data[last_time]
        
        # Extract coordinates and ages
        x_coords = [cell['x'] for cell in cells]
        y_coords = [cell['y'] for cell in cells]
        ages = [cell['age'] for cell in cells]
        
        ax = axes[i]
        
        # Create scatter plot colored by cell age
        scatter = ax.scatter(x_coords, y_coords, c=ages, s=60, alpha=0.7, 
                           cmap='viridis', edgecolors='black', linewidth=0.3)
        
        # Set labels and title
        ax.set_xlabel('X Position', fontsize=10)
        ax.set_ylabel('Y Position', fontsize=10)
        ax.set_title(f'{scenario_name}\\n(t={last_time:.2f}, {len(cells)} cells)', fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal', adjustable='box')
        
        # Add colorbar for the last subplot
        if i == len(scenarios) - 1:
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('Cell Age', fontsize=10)
    
    plt.suptitle('Organoid Formation: Comparison of Different Scenarios', fontsize=16)
    plt.tight_layout()
    plt.savefig('organoid_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("Comparison figure saved as 'organoid_comparison.png'")

def create_time_series_comparison():
    """
    Create a time series comparison for basement membrane scenario.
    """
    data_file = '/home/orlando/Thesis/testoutput/OrganoidFormation/WithBasementMembrane/results_from_time_0/cellages.dat'
    
    if not os.path.exists(data_file):
        print(f"Data file not found: {data_file}")
        return
    
    data = parse_chaste_data(data_file)
    if not data:
        print("No data found")
        return
    
    time_points = sorted(data.keys())
    n_times = len(time_points)
    
    fig, axes = plt.subplots(1, n_times, figsize=(4 * n_times, 4))
    if n_times == 1:
        axes = [axes]
    
    # Get global bounds for consistent scaling
    all_x = []
    all_y = []
    all_ages = []
    
    for time_point in time_points:
        cells = data[time_point]
        all_x.extend([cell['x'] for cell in cells])
        all_y.extend([cell['y'] for cell in cells])
        all_ages.extend([cell['age'] for cell in cells])
    
    x_min, x_max = min(all_x), max(all_x)
    y_min, y_max = min(all_y), max(all_y)
    age_min, age_max = min(all_ages), max(all_ages)
    
    # Add margins
    x_margin = (x_max - x_min) * 0.1 if x_max != x_min else 0.5
    y_margin = (y_max - y_min) * 0.1 if y_max != y_min else 0.5
    
    for i, time_point in enumerate(time_points):
        cells = data[time_point]
        
        x_coords = [cell['x'] for cell in cells]
        y_coords = [cell['y'] for cell in cells]
        ages = [cell['age'] for cell in cells]
        
        ax = axes[i]
        
        scatter = ax.scatter(x_coords, y_coords, c=ages, s=80, alpha=0.7, 
                           cmap='viridis', vmin=age_min, vmax=age_max,
                           edgecolors='black', linewidth=0.5)
        
        ax.set_xlim(x_min - x_margin, x_max + x_margin)
        ax.set_ylim(y_min - y_margin, y_max + y_margin)
        ax.set_xlabel('X Position', fontsize=10)
        ax.set_ylabel('Y Position', fontsize=10)
        ax.set_title(f't = {time_point:.2f}', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal', adjustable='box')
        
        # Add cell count annotation
        ax.text(0.02, 0.98, f'{len(cells)} cells', 
                transform=ax.transAxes, fontsize=9, 
                verticalalignment='top', 
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=axes[-1])
    cbar.set_label('Cell Age', fontsize=10)
    
    plt.suptitle('Organoid Formation Over Time (With Basement Membrane)', fontsize=14)
    plt.tight_layout()
    plt.savefig('organoid_time_series.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("Time series figure saved as 'organoid_time_series.png'")

if __name__ == "__main__":
    print("Creating organoid visualization comparisons...")
    create_comparison_figure()
    create_time_series_comparison()
    print("\\nVisualization files created:")
    print("- organoid_comparison.png: Comparison of different scenarios")
    print("- organoid_time_series.png: Time evolution with basement membrane")
    print("- organoid_formation.gif: Animation of formation process")
    print("- Individual scenario plots: low_stiffness.png, high_stiffness.png, organoid_t0.2.png")