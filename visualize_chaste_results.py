#!/usr/bin/env python3
"""
Visualization script for Chaste organoid simulation results.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import argparse

def parse_chaste_data(data_file):
    """
    Parse Chaste cell data file (like cellages.dat).
    Format: time_1 [cell_id x y age] [cell_id x y age] ... time_2 [cell_id x y age] ...
    """
    data = {}
    
    with open(data_file, 'r') as f:
        content = f.read().strip()
    
    # Replace newlines and tabs with spaces and split by whitespace
    tokens = content.replace('\n', ' ').replace('\t', ' ').split()
    
    i = 0
    while i < len(tokens):
        try:
            # Try to parse as time value
            time_val = float(tokens[i])
            
            # Look ahead to see if next token is '0' (first cell_id)
            if i + 1 < len(tokens) and tokens[i + 1] == '0':
                # This is a time stamp
                current_time = time_val
                cells = []
                i += 1  # Move to first cell_id
                
                # Parse cells until we hit next time stamp
                while i < len(tokens):
                    try:
                        # Try to parse cell data: cell_id x y age
                        if i + 3 < len(tokens):
                            cell_id = int(tokens[i])
                            x = float(tokens[i + 1])
                            y = float(tokens[i + 2])
                            age = float(tokens[i + 3])
                            
                            cells.append({
                                'id': cell_id,
                                'x': x,
                                'y': y,
                                'age': age
                            })
                            
                            i += 4
                            
                            # Check if next token could be a new time stamp
                            if i < len(tokens):
                                try:
                                    next_val = float(tokens[i])
                                    # If next token is a small decimal, likely a time stamp
                                    if next_val < 1.0 and i + 1 < len(tokens) and tokens[i + 1] == '0':
                                        # Next time stamp found, finish current time
                                        break
                                except ValueError:
                                    pass
                        else:
                            break
                            
                    except (ValueError, IndexError):
                        break
                
                if cells:
                    data[current_time] = cells
                    
            else:
                i += 1
                
        except (ValueError, IndexError):
            i += 1
    
    return data

def plot_organoid_snapshot(data, time_point, title="Organoid Formation"):
    """
    Plot a single time point of the organoid simulation.
    """
    if time_point not in data:
        print(f"Time point {time_point} not found in data")
        return
    
    cells = data[time_point]
    
    # Extract coordinates and ages
    x_coords = [cell['x'] for cell in cells]
    y_coords = [cell['y'] for cell in cells]
    ages = [cell['age'] for cell in cells]
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create scatter plot colored by cell age
    scatter = ax.scatter(x_coords, y_coords, c=ages, s=100, alpha=0.7, 
                        cmap='viridis', edgecolors='black', linewidth=0.5)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Cell Age', fontsize=12)
    
    # Set labels and title
    ax.set_xlabel('X Position', fontsize=12)
    ax.set_ylabel('Y Position', fontsize=12)
    ax.set_title(f'{title} at t={time_point}', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal', adjustable='box')
    
    # Add cell count annotation
    ax.text(0.02, 0.98, f'Cells: {len(cells)}', 
            transform=ax.transAxes, fontsize=10, 
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    return fig, ax

def create_organoid_animation(data, output_file="organoid_animation.gif"):
    """
    Create an animation showing organoid formation over time.
    """
    if not data:
        print("No data to animate")
        return
    
    time_points = sorted(data.keys())
    
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
    x_margin = (x_max - x_min) * 0.1
    y_margin = (y_max - y_min) * 0.1
    
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('X Position', fontsize=12)
    ax.set_ylabel('Y Position', fontsize=12)
    
    # Initialize scatter plot
    scat = ax.scatter([], [], c=[], s=100, alpha=0.7, 
                     cmap='viridis', vmin=age_min, vmax=age_max,
                     edgecolors='black', linewidth=0.5)
    
    # Add colorbar
    cbar = plt.colorbar(scat, ax=ax)
    cbar.set_label('Cell Age', fontsize=12)
    
    # Title and cell count text
    title = ax.set_title('', fontsize=14)
    cell_count_text = ax.text(0.02, 0.98, '', transform=ax.transAxes, fontsize=10,
                             verticalalignment='top', 
                             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    def animate(frame):
        time_point = time_points[frame]
        cells = data[time_point]
        
        if cells:
            x_coords = [cell['x'] for cell in cells]
            y_coords = [cell['y'] for cell in cells]
            ages = [cell['age'] for cell in cells]
            
            # Update scatter plot
            scat.set_offsets(np.column_stack([x_coords, y_coords]))
            scat.set_array(np.array(ages))
            
            # Update title and cell count
            title.set_text(f'Organoid Formation at t={time_point:.2f}')
            cell_count_text.set_text(f'Cells: {len(cells)}')
        
        return scat, title, cell_count_text
    
    # Create animation
    anim = FuncAnimation(fig, animate, frames=len(time_points), 
                        interval=500, blit=False, repeat=True)
    
    # Save animation
    print(f"Saving animation to {output_file}...")
    anim.save(output_file, writer='pillow', fps=2)
    
    plt.tight_layout()
    return fig, anim

def main():
    parser = argparse.ArgumentParser(description='Visualize Chaste organoid simulation results')
    parser.add_argument('data_file', help='Path to cellages.dat or similar Chaste output file')
    parser.add_argument('--animate', action='store_true', help='Create animation instead of static plot')
    parser.add_argument('--output', '-o', default=None, help='Output file name')
    parser.add_argument('--time', '-t', type=float, default=None, help='Time point to plot (for static plots)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.data_file):
        print(f"Error: File {args.data_file} not found")
        return 1
    
    print(f"Loading data from {args.data_file}...")
    data = parse_chaste_data(args.data_file)
    
    if not data:
        print("No data found in file")
        return 1
    
    time_points = sorted(data.keys())
    print(f"Found {len(time_points)} time points: {time_points}")
    
    if args.animate:
        output_file = args.output or "organoid_animation.gif"
        create_organoid_animation(data, output_file)
        print(f"Animation saved to {output_file}")
    else:
        # Static plot
        if args.time is not None:
            time_point = args.time
        else:
            time_point = time_points[-1]  # Use last time point by default
        
        fig, ax = plot_organoid_snapshot(data, time_point)
        
        if args.output:
            plt.savefig(args.output, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {args.output}")
        else:
            plt.show()
    
    return 0

if __name__ == "__main__":
    sys.exit(main())