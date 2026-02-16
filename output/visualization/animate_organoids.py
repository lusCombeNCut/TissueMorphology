#!/usr/bin/env python3
"""
Real-time visualization of basement membrane effects in organoid formation.
This script creates an animated view of how cells move under basement membrane constraints.
"""

import vtk
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import glob
from vtk.util.numpy_support import vtk_to_numpy

def create_basement_membrane_animation(sim_dir="testoutput/OrganoidFormation/WithBasementMembrane"):
    """Create an animation showing basement membrane effects."""
    
    # Find VTK files
    pattern = os.path.join(sim_dir, "results_from_time_*", "results.vtu")
    vtk_files = sorted(glob.glob(pattern))
    
    if not vtk_files:
        print(f"No VTK files found in {sim_dir}")
        return
    
    print(f"Found {len(vtk_files)} time points for animation")
    
    # Read all data first
    all_points = []
    all_cell_data = []
    times = []
    
    for vtk_file in vtk_files:
        try:
            # Read VTK file
            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(vtk_file)
            reader.Update()
            
            output = reader.GetOutput()
            points = vtk_to_numpy(output.GetPoints().GetData())
            
            # Get cell data
            cell_data = {}
            point_data = output.GetPointData()
            
            for i in range(point_data.GetNumberOfArrays()):
                array_name = point_data.GetArrayName(i)
                array_data = vtk_to_numpy(point_data.GetArray(i))
                cell_data[array_name] = array_data
            
            # Extract time
            time_str = vtk_file.split('results_from_time_')[1].split('/')[0]
            time = float(time_str)
            
            all_points.append(points)
            all_cell_data.append(cell_data)
            times.append(time)
            
        except Exception as e:
            print(f"Error reading {vtk_file}: {e}")
            continue
    
    if not all_points:
        print("No valid data found")
        return
    
    # Create figure and animation
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Calculate consistent axis limits
    all_x = np.concatenate([pts[:, 0] for pts in all_points])
    all_y = np.concatenate([pts[:, 1] for pts in all_points])
    x_margin = (all_x.max() - all_x.min()) * 0.1
    y_margin = (all_y.max() - all_y.min()) * 0.1
    
    ax1.set_xlim(all_x.min() - x_margin, all_x.max() + x_margin)
    ax1.set_ylim(all_y.min() - y_margin, all_y.max() + y_margin)
    ax2.set_xlim(all_x.min() - x_margin, all_x.max() + x_margin)
    ax2.set_ylim(all_y.min() - y_margin, all_y.max() + y_margin)
    
    def animate(frame):
        ax1.clear()
        ax2.clear()
        
        if frame >= len(all_points):
            return
        
        points = all_points[frame]
        cell_data = all_cell_data[frame]
        time = times[frame]
        
        # Plot 1: Cell positions colored by basement membrane stiffness
        if 'basement_membrane_stiffness' in cell_data:
            stiffness = cell_data['basement_membrane_stiffness']
            scatter1 = ax1.scatter(points[:, 0], points[:, 1], c=stiffness, 
                                 s=60, alpha=0.8, cmap='viridis',
                                 edgecolors='black', linewidth=0.5)
            ax1.set_title(f'Basement Membrane Stiffness\nt = {time:.3f}')
        else:
            ax1.scatter(points[:, 0], points[:, 1], c='lightblue', 
                       s=60, alpha=0.7, edgecolors='black', linewidth=0.5)
            ax1.set_title(f'Cell Positions\nt = {time:.3f}')
        
        # Plot 2: Cell types
        if 'cell_type' in cell_data:
            cell_type = cell_data['cell_type']
            colors = ['red' if ct == 0.0 else 'blue' for ct in cell_type]
            ax2.scatter(points[:, 0], points[:, 1], c=colors, 
                       s=60, alpha=0.8, edgecolors='black', linewidth=0.5)
            ax2.set_title(f'Cell Types (Red: Stem, Blue: Diff)\nt = {time:.3f}')
        else:
            ax2.scatter(points[:, 0], points[:, 1], c='lightgreen', 
                       s=60, alpha=0.7, edgecolors='black', linewidth=0.5)
            ax2.set_title(f'Cell Positions\nt = {time:.3f}')
        
        # Draw basement membrane boundary (approximate)
        if frame == 0:  # Only show on first frame for reference
            center_x, center_y = np.mean(points[:, 0]), np.mean(points[:, 1])
            radius = 3.0  # BasementMembraneRadius from test
            circle = plt.Circle((center_x, center_y), radius, fill=False, 
                              color='red', linestyle='--', linewidth=2, alpha=0.7)
            ax1.add_patch(circle)
            circle2 = plt.Circle((center_x, center_y), radius, fill=False, 
                               color='red', linestyle='--', linewidth=2, alpha=0.7)
            ax2.add_patch(circle2)
        
        # Formatting
        for ax in [ax1, ax2]:
            ax.set_xlabel('X position')
            ax.set_ylabel('Y position')
            ax.set_aspect('equal')
            ax.grid(True, alpha=0.3)
            ax.set_xlim(all_x.min() - x_margin, all_x.max() + x_margin)
            ax.set_ylim(all_y.min() - y_margin, all_y.max() + y_margin)
    
    # Create animation
    anim = animation.FuncAnimation(fig, animate, frames=len(all_points), 
                                 interval=500, repeat=True, blit=False)
    
    plt.suptitle('Organoid Formation with Basement Membrane Forces', fontsize=14)
    plt.tight_layout()
    
    # Save animation as gif
    try:
        anim.save('organoid_basement_membrane_animation.gif', writer='pillow', fps=2)
        print("Animation saved as: organoid_basement_membrane_animation.gif")
    except:
        print("Could not save animation (pillow not available)")
    
    plt.show()
    return anim

def plot_cell_trajectory_analysis(sim_dir="testoutput/OrganoidFormation/WithBasementMembrane"):
    """Analyze and plot cell movement trajectories."""
    
    pattern = os.path.join(sim_dir, "results_from_time_*", "results.vtu")
    vtk_files = sorted(glob.glob(pattern))
    
    if len(vtk_files) < 2:
        print("Need at least 2 time points for trajectory analysis")
        return
    
    # Track cell movements
    all_positions = []
    times = []
    
    for vtk_file in vtk_files:
        try:
            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(vtk_file)
            reader.Update()
            
            output = reader.GetOutput()
            points = vtk_to_numpy(output.GetPoints().GetData())
            
            time_str = vtk_file.split('results_from_time_')[1].split('/')[0]
            time = float(time_str)
            
            all_positions.append(points)
            times.append(time)
            
        except Exception as e:
            continue
    
    if len(all_positions) < 2:
        return
    
    # Calculate displacement vectors
    initial_pos = all_positions[0]
    final_pos = all_positions[-1]
    
    displacement = final_pos - initial_pos
    displacement_mag = np.linalg.norm(displacement, axis=1)
    
    # Calculate center of mass movement
    center_initial = np.mean(initial_pos, axis=0)
    center_final = np.mean(final_pos, axis=0)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot 1: Cell trajectories
    for i in range(0, len(initial_pos), max(1, len(initial_pos)//20)):  # Sample trajectories
        traj_x = [pos[i, 0] for pos in all_positions]
        traj_y = [pos[i, 1] for pos in all_positions]
        ax1.plot(traj_x, traj_y, alpha=0.6, linewidth=1)
    
    ax1.scatter(initial_pos[:, 0], initial_pos[:, 1], c='green', s=30, 
               alpha=0.7, label='Initial', edgecolors='black', linewidth=0.5)
    ax1.scatter(final_pos[:, 0], final_pos[:, 1], c='red', s=30, 
               alpha=0.7, label='Final', edgecolors='black', linewidth=0.5)
    
    # Draw basement membrane boundary
    radius = 3.0
    circle = plt.Circle(center_initial, radius, fill=False, 
                       color='blue', linestyle='--', linewidth=2, alpha=0.7)
    ax1.add_patch(circle)
    
    ax1.set_title('Cell Trajectories')
    ax1.set_xlabel('X position')
    ax1.set_ylabel('Y position')
    ax1.legend()
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Displacement magnitude
    scatter = ax2.scatter(initial_pos[:, 0], initial_pos[:, 1], 
                         c=displacement_mag, s=50, alpha=0.8, 
                         cmap='plasma', edgecolors='black', linewidth=0.5)
    plt.colorbar(scatter, ax=ax2, label='Displacement Magnitude')
    
    ax2.set_title('Cell Displacement Magnitude')
    ax2.set_xlabel('X position')
    ax2.set_ylabel('Y position')
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    
    plt.suptitle(f'Organoid Cell Movement Analysis\nTime: {times[0]:.3f} to {times[-1]:.3f}', 
                fontsize=14)
    plt.tight_layout()
    
    plt.savefig("organoid_trajectory_analysis.png", dpi=150, bbox_inches='tight')
    print("Trajectory analysis saved as: organoid_trajectory_analysis.png")
    plt.show()

if __name__ == "__main__":
    print("Creating basement membrane animation...")
    create_basement_membrane_animation()
    
    print("\nAnalyzing cell trajectories...")
    plot_cell_trajectory_analysis()