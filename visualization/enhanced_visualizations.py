#!/usr/bin/env python3
"""
Enhanced visualization for large organoid datasets with 2D and 3D support
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def parse_chaste_data_memory_efficient(data_file, max_timepoints=None, skip_timepoints=None):
    """
    Memory-efficient parser that can skip timepoints or limit total timepoints
    """
    data = {}
    
    with open(data_file, 'r') as f:
        content = f.read().strip()
    
    tokens = content.replace('\n', ' ').replace('\t', ' ').split()
    
    i = 0
    timepoint_count = 0
    
    while i < len(tokens):
        try:
            time_val = float(tokens[i])
            
            # Check if this looks like a timepoint header
            if i + 1 < len(tokens) and tokens[i + 1] == '0':
                # Skip if we've reached max timepoints
                if max_timepoints and timepoint_count >= max_timepoints:
                    break
                
                # Skip if we're skipping timepoints
                if skip_timepoints and (timepoint_count % (skip_timepoints + 1)) != 0:
                    timepoint_count += 1
                    # Skip ahead to find next timepoint
                    i += 1
                    continue
                
                # Parse this timepoint
                timepoint_count += 1
                current_time = time_val
                data[current_time] = []
                i += 1  # Move past time
                
                # Parse cells for this timepoint
                while i < len(tokens):
                    try:
                        if i + 4 >= len(tokens):
                            break
                            
                        # Try to parse next timepoint
                        next_time = float(tokens[i + 4])
                        if i + 5 < len(tokens) and tokens[i + 5] == '0':
                            # This is the start of next timepoint
                            break
                        
                        # Parse cell data
                        cell_id = int(tokens[i])
                        x = float(tokens[i + 1])
                        y = float(tokens[i + 2]) 
                        age = float(tokens[i + 3])
                        
                        data[current_time].append({
                            'id': cell_id, 'x': x, 'y': y, 'age': age
                        })
                        i += 4
                        
                    except (ValueError, IndexError):
                        i += 1
                        break
            else:
                i += 1
                
        except (ValueError, IndexError):
            i += 1
    
    return data

def parse_3d_chaste_data(data_file, max_timepoints=4):
    """
    Parse 3D Chaste data (expects cell_id x y z age format)
    """
    data = {}
    
    with open(data_file, 'r') as f:
        content = f.read().strip()
    
    # Split into tokens
    tokens = content.replace('\n', ' ').replace('\t', ' ').split()
    
    i = 0
    timepoint_count = 0
    
    while i < len(tokens):
        try:
            # Try to parse as time value
            time_val = float(tokens[i])
            
            # Check if this looks like a timepoint header (followed by cell_id 0)
            if i + 1 < len(tokens) and tokens[i + 1] == '0':
                if max_timepoints and timepoint_count >= max_timepoints:
                    break
                    
                current_time = time_val
                data[current_time] = []
                timepoint_count += 1
                i += 1  # Move past time
                
                # Parse cells for this timepoint
                while i < len(tokens):
                    try:
                        if i + 4 >= len(tokens):
                            break
                            
                        # Check if we're at the next timepoint
                        if i + 5 < len(tokens):
                            try:
                                next_time = float(tokens[i + 5])
                                if tokens[i + 6] == '0':  # Next timepoint starts
                                    break
                            except (ValueError, IndexError):
                                pass
                        
                        # Parse 3D cell data: cell_id x y z age
                        cell_id = int(tokens[i])
                        x = float(tokens[i + 1])
                        y = float(tokens[i + 2])
                        z = float(tokens[i + 3])
                        age = float(tokens[i + 4])
                        
                        data[current_time].append({
                            'id': cell_id, 'x': x, 'y': y, 'z': z, 'age': age
                        })
                        i += 5
                        
                    except (ValueError, IndexError):
                        i += 1
                        break
            else:
                i += 1
                
        except (ValueError, IndexError):
            i += 1
    
    return data

def create_2d_organoid_visualization():
    """Create visualization for 2D organoid tests"""
    
    scenarios = {
        'Basic Test': '/home/orlando/Thesis/testoutput/OrganoidFormation/BasicTest/results_from_time_0/cellages.dat',
        'Low Stiffness': '/home/orlando/Thesis/testoutput/OrganoidFormation/LowStiffness/results_from_time_0/cellages.dat', 
        'High Stiffness': '/home/orlando/Thesis/testoutput/OrganoidFormation/HighStiffness/results_from_time_0/cellages.dat',
        'With Basement Membrane': '/home/orlando/Thesis/testoutput/OrganoidFormation/WithBasementMembrane/results_from_time_0/cellages.dat'
    }
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    axes = axes.flatten()
    
    print("Creating 2D organoid comparison...")
    
    for i, (scenario_name, data_file) in enumerate(scenarios.items()):
        if not os.path.exists(data_file):
            print(f"Warning: Data file not found for {scenario_name}: {data_file}")
            continue
        
        # Parse only final timepoint to save memory
        data = parse_chaste_data_memory_efficient(data_file, max_timepoints=1)
        if not data:
            print(f"Warning: No data parsed for {scenario_name}")
            continue
        
        # Get final time point
        final_time = max(data.keys())
        cells = data[final_time]
        
        # Extract coordinates and ages
        x_coords = [cell['x'] for cell in cells]
        y_coords = [cell['y'] for cell in cells]
        ages = [cell['age'] for cell in cells]
        
        ax = axes[i]
        
        # Create scatter plot
        scatter = ax.scatter(x_coords, y_coords, c=ages, s=40, alpha=0.8,
                           cmap='viridis', edgecolors='black', linewidth=0.2)
        
        ax.set_xlabel('X Position', fontsize=10)
        ax.set_ylabel('Y Position', fontsize=10) 
        ax.set_title(f'{scenario_name}\\n(t={final_time:.2f}, {len(cells)} cells)', fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal', adjustable='box')
        
        # Add colorbar for last subplot
        if i == len(scenarios) - 1:
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('Cell Age', fontsize=10)
    
    plt.suptitle('2D Organoid Formation: Comparison of Different Scenarios', fontsize=16)
    plt.tight_layout()
    plt.savefig('organoid_comparison_2d.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("2D comparison saved as 'organoid_comparison_2d.png'")

def create_3d_organoid_visualization():
    """Create visualization for 3D organoid tests"""
    
    scenarios_3d = {
        'Spherical Formation': '/home/orlando/Thesis/testoutput/Organoid3d/SphericalFormation/results_from_time_0/cellages.dat',
        'Stiffness Test': '/home/orlando/Thesis/testoutput/Organoid3d/StiffnessTest/results_from_time_0/cellages.dat',
        'Threshold ECM Force': '/home/orlando/Thesis/testoutput/Organoid3d/ThresholdECMForce/results_from_time_0/cellages.dat',
        '2D Threshold ECM Force': '/home/orlando/Thesis/testoutput/Organoid2d/ThresholdECMForce/results_from_time_0/cellages.dat'
    }
    
    print("Creating 3D organoid visualizations...")
    
    for scenario_name, data_file in scenarios_3d.items():
        if not os.path.exists(data_file):
            print(f"Warning: 3D data file not found: {data_file}")
            continue
        
        # Parse 3D data
        data = parse_3d_chaste_data(data_file, max_timepoints=2)
        if not data:
            print(f"Warning: No 3D data parsed for {scenario_name}")
            continue
        
        # Create 3D plot
        fig = plt.figure(figsize=(12, 5))
        
        time_points = sorted(data.keys())
        
        for i, time_point in enumerate(time_points[:2]):  # Show first and last
            cells = data[time_point]
            
            if len(cells) == 0:
                continue
                
            # Check if data has z coordinate (3D)
            if 'z' in cells[0]:
                # 3D plot
                ax = fig.add_subplot(1, 2, i+1, projection='3d')
                
                x_coords = [cell['x'] for cell in cells]
                y_coords = [cell['y'] for cell in cells] 
                z_coords = [cell['z'] for cell in cells]
                ages = [cell['age'] for cell in cells]
                
                scatter = ax.scatter(x_coords, y_coords, z_coords, c=ages, s=30,
                                   cmap='plasma', alpha=0.7)
                
                ax.set_xlabel('X Position')
                ax.set_ylabel('Y Position')
                ax.set_zlabel('Z Position')
                ax.set_title(f'{scenario_name}\\nt={time_point:.1f} ({len(cells)} cells)')
                
                if i == len(time_points[:2]) - 1:
                    plt.colorbar(scatter, ax=ax, label='Cell Age', shrink=0.6)
            else:
                # 2D plot fallback
                ax = fig.add_subplot(1, 2, i+1)
                
                x_coords = [cell['x'] for cell in cells]
                y_coords = [cell['y'] for cell in cells]
                ages = [cell['age'] for cell in cells]
                
                scatter = ax.scatter(x_coords, y_coords, c=ages, s=40, cmap='plasma')
                ax.set_xlabel('X Position')
                ax.set_ylabel('Y Position') 
                ax.set_title(f't={time_point:.1f} ({len(cells)} cells)')
                ax.set_aspect('equal')
        
        plt.tight_layout()
        filename = f'organoid_3d_{scenario_name.lower().replace(" ", "_")}.png'
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"3D visualization saved as '{filename}'")

def main():
    """Main function to create all visualizations"""
    print("Creating enhanced organoid visualizations...")
    
    # Create 2D visualizations
    create_2d_organoid_visualization()
    
    # Create 3D visualizations  
    create_3d_organoid_visualization()
    
    print("\\nVisualization files created:")
    print("- organoid_comparison_2d.png: 2D organoid comparison")
    print("- organoid_3d_*.png: Individual 3D organoid plots")
    print("\\nFor interactive 3D analysis, use ParaView with the .pvd files:")
    print("paraview /path/to/results.pvd")

if __name__ == "__main__":
    main()