#!/usr/bin/env python3
"""
Simple visualization for large organoid datasets
"""

import numpy as np
import matplotlib.pyplot as plt
import os

def parse_cellages_simple(filename):
    """Parse cellages.dat and return final timepoint only"""
    data = {}
    current_time = None
    
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 0:
                continue
                
            # Check if this is a time header
            try:
                time_val = float(parts[0])
                if len(parts) == 1:  # Time header line
                    current_time = time_val
                    data[current_time] = []
                else:  # Data line with time
                    if len(parts) >= 5:
                        cell_id = int(parts[1])
                        x = float(parts[2])
                        y = float(parts[3])
                        age = float(parts[4])
                        
                        if current_time not in data:
                            data[current_time] = []
                        data[current_time].append({'id': cell_id, 'x': x, 'y': y, 'age': age})
            except:
                continue
    
    return data

def plot_final_state():
    """Plot just the final state of the organoid"""
    
    data_file = '/home/orlando/Thesis/testoutput/OrganoidFormation/WithBasementMembrane/results_from_time_0/cellages.dat'
    
    if not os.path.exists(data_file):
        print(f"Data file not found: {data_file}")
        return
    
    print("Parsing organoid data...")
    data = parse_cellages_simple(data_file)
    
    if not data:
        print("No data found!")
        return
    
    # Get final timepoint
    final_time = max(data.keys())
    cells = data[final_time]
    
    print(f"Found {len(cells)} cells at final time {final_time}")
    
    # Extract data
    x_coords = [cell['x'] for cell in cells]
    y_coords = [cell['y'] for cell in cells]
    ages = [cell['age'] for cell in cells]
    
    # Create plot
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(x_coords, y_coords, c=ages, cmap='viridis', s=50, alpha=0.7)
    plt.colorbar(scatter, label='Cell Age')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    plt.title(f'Organoid Formation Final State (t={final_time})\n{len(cells)} cells, 100 timesteps')
    plt.grid(True, alpha=0.3)
    plt.axis('equal')
    
    # Save plot
    plt.savefig('organoid_large_final.png', dpi=150, bbox_inches='tight')
    print("Saved: organoid_large_final.png")
    plt.close()

if __name__ == "__main__":
    plot_final_state()