#!/usr/bin/env python3
"""
Simple visualization script for TissueMorphology test results.
This script reads VTK files from Chaste simulations and creates basic plots.
"""

import vtk
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from vtk.util.numpy_support import vtk_to_numpy

class OrganoidVisualizer:
    def __init__(self, results_dir="testoutput"):
        self.results_dir = results_dir
        
    def find_simulation_directories(self):
        """Find all OrganoidFormation simulation directories."""
        pattern = os.path.join(self.results_dir, "OrganoidFormation", "*")
        dirs = glob.glob(pattern)
        return [d for d in dirs if os.path.isdir(d)]
    
    def get_vtk_files(self, sim_dir):
        """Get all VTK files from a simulation directory."""
        pattern = os.path.join(sim_dir, "results_from_time_*", "results.vtu")
        vtk_files = glob.glob(pattern)
        return sorted(vtk_files)
    
    def read_vtk_file(self, filename):
        """Read a VTK file and extract cell data."""
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(filename)
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
            
        return points, cell_data
    
    def plot_simulation_overview(self, sim_dir):
        """Create an overview plot of a simulation."""
        vtk_files = self.get_vtk_files(sim_dir)
        
        if not vtk_files:
            print(f"No VTK files found in {sim_dir}")
            return
            
        sim_name = os.path.basename(sim_dir)
        print(f"Visualizing simulation: {sim_name}")
        
        # Create subplots for different time points
        n_files = min(len(vtk_files), 4)  # Show max 4 time points
        fig, axes = plt.subplots(1, n_files, figsize=(4*n_files, 4))
        if n_files == 1:
            axes = [axes]
            
        for i, vtk_file in enumerate(vtk_files[:n_files]):
            try:
                points, cell_data = self.read_vtk_file(vtk_file)
                
                # Extract time from filename
                time_str = vtk_file.split('results_from_time_')[1].split('/')[0]
                time = float(time_str)
                
                # Plot cell positions
                ax = axes[i]
                scatter = ax.scatter(points[:, 0], points[:, 1], c='lightblue', 
                                   s=50, alpha=0.7, edgecolors='black', linewidth=0.5)
                
                # Color by basement membrane stiffness if available
                if 'basement_membrane_stiffness' in cell_data:
                    stiffness = cell_data['basement_membrane_stiffness']
                    scatter = ax.scatter(points[:, 0], points[:, 1], c=stiffness, 
                                       s=50, alpha=0.8, cmap='viridis', 
                                       edgecolors='black', linewidth=0.5)
                    plt.colorbar(scatter, ax=ax, label='Basement Membrane Stiffness')
                
                # Color by cell type if available
                elif 'cell_type' in cell_data:
                    cell_type = cell_data['cell_type']
                    colors = ['red' if ct == 0.0 else 'blue' for ct in cell_type]
                    ax.scatter(points[:, 0], points[:, 1], c=colors, 
                             s=50, alpha=0.8, edgecolors='black', linewidth=0.5)
                
                ax.set_title(f't = {time:.3f}')
                ax.set_xlabel('X position')
                ax.set_ylabel('Y position')
                ax.set_aspect('equal')
                ax.grid(True, alpha=0.3)
                
            except Exception as e:
                print(f"Error processing {vtk_file}: {e}")
                continue
        
        plt.suptitle(f'Organoid Formation - {sim_name}', fontsize=14)
        plt.tight_layout()
        
        # Save plot
        output_file = f"organoid_visualization_{sim_name.lower()}.png"
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved visualization: {output_file}")
        plt.show()
    
    def create_summary_plot(self):
        """Create a summary plot comparing all simulations."""
        sim_dirs = self.find_simulation_directories()
        
        if not sim_dirs:
            print("No simulation directories found")
            return
            
        print(f"Found {len(sim_dirs)} simulation directories")
        
        fig, axes = plt.subplots(1, len(sim_dirs), figsize=(4*len(sim_dirs), 4))
        if len(sim_dirs) == 1:
            axes = [axes]
            
        for i, sim_dir in enumerate(sim_dirs):
            vtk_files = self.get_vtk_files(sim_dir)
            sim_name = os.path.basename(sim_dir)
            
            if not vtk_files:
                continue
                
            # Use the last time point
            try:
                points, cell_data = self.read_vtk_file(vtk_files[-1])
                
                ax = axes[i]
                
                # Plot based on available data
                if 'basement_membrane_stiffness' in cell_data:
                    stiffness = cell_data['basement_membrane_stiffness']
                    scatter = ax.scatter(points[:, 0], points[:, 1], c=stiffness, 
                                       s=50, alpha=0.8, cmap='viridis',
                                       edgecolors='black', linewidth=0.5)
                    plt.colorbar(scatter, ax=ax, label='Stiffness')
                elif 'cell_type' in cell_data:
                    cell_type = cell_data['cell_type']
                    colors = ['red' if ct == 0.0 else 'blue' for ct in cell_type]
                    ax.scatter(points[:, 0], points[:, 1], c=colors, 
                             s=50, alpha=0.8, edgecolors='black', linewidth=0.5)
                else:
                    ax.scatter(points[:, 0], points[:, 1], c='lightblue', 
                             s=50, alpha=0.7, edgecolors='black', linewidth=0.5)
                
                ax.set_title(sim_name)
                ax.set_xlabel('X position')
                ax.set_ylabel('Y position')
                ax.set_aspect('equal')
                ax.grid(True, alpha=0.3)
                
            except Exception as e:
                print(f"Error processing {sim_dir}: {e}")
                continue
        
        plt.suptitle('Organoid Formation - All Simulations (Final States)', fontsize=14)
        plt.tight_layout()
        
        # Save summary plot
        plt.savefig("organoid_summary.png", dpi=150, bbox_inches='tight')
        print("Saved summary visualization: organoid_summary.png")
        plt.show()


def main():
    """Main function to visualize organoid formation results."""
    visualizer = OrganoidVisualizer()
    
    # Create summary plot
    visualizer.create_summary_plot()
    
    # Create individual simulation plots
    sim_dirs = visualizer.find_simulation_directories()
    for sim_dir in sim_dirs:
        visualizer.plot_simulation_overview(sim_dir)


if __name__ == "__main__":
    main()