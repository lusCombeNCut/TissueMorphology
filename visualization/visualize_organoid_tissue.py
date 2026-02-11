#!/usr/bin/env python3
"""
Enhanced 3D Organoid Visualization
Creates tissue-like surfaces from discrete cell populations
"""

import vtk
import numpy as np
import os
import glob
import argparse

class EnhancedOrganoid3DVisualizer:
    def __init__(self, data_dir):
        self.data_dir = data_dir
        self.renderer = vtk.vtkRenderer()
        self.render_window = vtk.vtkRenderWindow()
        self.render_window_interactor = vtk.vtkRenderWindowInteractor()
        
        # Set up the visualization pipeline
        self.setup_renderer()
        
    def setup_renderer(self):
        """Initialize the VTK renderer and window"""
        self.render_window.AddRenderer(self.renderer)
        self.render_window_interactor.SetRenderWindow(self.render_window)
        
        # Set background color
        self.renderer.SetBackground(0.1, 0.1, 0.2)  # Dark blue background
        
        # Set window size
        self.render_window.SetSize(1200, 800)
        self.render_window.SetWindowName("3D Organoid Tissue Visualization")
        
    def load_vtk_file(self, filename):
        """Load a VTU file and return the dataset"""
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(filename)
        reader.Update()
        return reader.GetOutput()
        
    def create_tissue_surface(self, dataset):
        """Create a continuous tissue surface from discrete cells"""
        points = dataset.GetPoints()
        if not points or points.GetNumberOfPoints() == 0:
            return None
            
        # Convert to numpy array
        coords = []
        for i in range(points.GetNumberOfPoints()):
            coords.append(points.GetPoint(i))
        coords = np.array(coords)
        
        if len(coords) < 4:
            return None
            
        # Create VTK points
        vtk_points = vtk.vtkPoints()
        for coord in coords:
            vtk_points.InsertNextPoint(coord)
            
        # Create polydata
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(vtk_points)
        
        # Method 1: Convex hull (simple but not always realistic)
        hull = vtk.vtkDelaunay3D()
        hull.SetInputData(polydata)
        hull.Update()
        
        # Extract surface
        surface_filter = vtk.vtkDataSetSurfaceFilter()
        surface_filter.SetInputData(hull.GetOutput())
        surface_filter.Update()
        
        # Smooth the surface
        smoother = vtk.vtkSmoothPolyDataFilter()
        smoother.SetInputData(surface_filter.GetOutput())
        smoother.SetNumberOfIterations(50)
        smoother.SetRelaxationFactor(0.1)
        smoother.FeatureEdgeSmoothingOff()
        smoother.BoundarySmoothingOn()
        smoother.Update()
        
        return smoother.GetOutput()
        
    def create_metaball_surface(self, dataset, radius=8.0):
        """Create smooth organoid surface using metaballs/implicit surface"""
        points = dataset.GetPoints()
        if not points or points.GetNumberOfPoints() == 0:
            return None
            
        # Create implicit function (sum of spheres)
        implicit_sum = vtk.vtkImplicitSum()
        
        for i in range(points.GetNumberOfPoints()):
            point = points.GetPoint(i)
            
            # Create sphere for each cell
            sphere = vtk.vtkSphere()
            sphere.SetCenter(point)
            sphere.SetRadius(radius)
            
            implicit_sum.AddFunction(sphere)
            
        # Create isosurface
        sample = vtk.vtkSampleFunction()
        sample.SetImplicitFunction(implicit_sum)
        
        # Calculate bounds
        bounds = [0, 0, 0, 0, 0, 0]
        points.GetBounds(bounds)
        
        # Expand bounds
        margin = 20.0
        sample.SetModelBounds(
            bounds[0] - margin, bounds[1] + margin,
            bounds[2] - margin, bounds[3] + margin, 
            bounds[4] - margin, bounds[5] + margin
        )
        sample.SetSampleDimensions(50, 50, 50)
        sample.Update()
        
        # Extract surface
        contour = vtk.vtkContourFilter()
        contour.SetInputData(sample.GetOutput())
        contour.SetValue(0, points.GetNumberOfPoints() * 0.3)  # Adjust threshold
        contour.Update()
        
        # Smooth result
        smoother = vtk.vtkSmoothPolyDataFilter()
        smoother.SetInputData(contour.GetOutput())
        smoother.SetNumberOfIterations(30)
        smoother.SetRelaxationFactor(0.1)
        smoother.Update()
        
        return smoother.GetOutput()
        
    def create_organoid_actor(self, surface_data, style="tissue"):
        """Create actor for organoid surface"""
        if not surface_data:
            return None
            
        # Calculate normals for better lighting
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputData(surface_data)
        normals.ComputePointNormalsOn()
        normals.ComputeCellNormalsOn()
        normals.Update()
        
        # Create mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(normals.GetOutput())
        
        # Create actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        
        if style == "tissue":
            # Tissue-like appearance
            actor.GetProperty().SetColor(0.9, 0.6, 0.6)  # Pink tissue color
            actor.GetProperty().SetSpecular(0.4)
            actor.GetProperty().SetSpecularPower(20)
            actor.GetProperty().SetOpacity(0.8)
        elif style == "membrane":
            # Semi-transparent membrane
            actor.GetProperty().SetColor(0.8, 0.9, 0.6)  # Light green
            actor.GetProperty().SetSpecular(0.6)
            actor.GetProperty().SetSpecularPower(40)
            actor.GetProperty().SetOpacity(0.6)
        else:
            # Default
            actor.GetProperty().SetColor(0.7, 0.7, 0.9)
            
        return actor
        
    def add_cell_centers(self, dataset, show_centers=True):
        """Optionally show individual cell centers as small spheres"""
        if not show_centers:
            return []
            
        points = dataset.GetPoints()
        actors = []
        
        for i in range(min(points.GetNumberOfPoints(), 50)):  # Limit to 50 for clarity
            point = points.GetPoint(i)
            
            sphere = vtk.vtkSphereSource()
            sphere.SetCenter(point)
            sphere.SetRadius(1.0)  # Small spheres
            sphere.SetPhiResolution(8)
            sphere.SetThetaResolution(8)
            
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(sphere.GetOutputPort())
            
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor(1.0, 1.0, 0.0)  # Yellow centers
            actor.GetProperty().SetOpacity(0.7)
            
            actors.append(actor)
            
        return actors
        
    def add_coordinate_axes(self):
        """Add coordinate axes"""
        axes = vtk.vtkAxesActor()
        axes.SetTotalLength(30, 30, 30)
        axes.SetShaftType(0)
        axes.SetAxisLabels(1)
        
        self.renderer.AddActor(axes)
        
    def add_text_overlay(self, metrics, time_point=0):
        """Add text overlay with organoid metrics"""
        if not metrics:
            return
            
        text = f"3D Organoid Tissue View\n"
        text += f"Time: {time_point}\n"
        text += f"Cells: {metrics['num_cells']}\n"
        text += f"Radius: {metrics['mean_radius']:.1f}μm\n"
        text += f"Sphericity: {metrics['sphericity']:.3f}"
        
        text_actor = vtk.vtkTextActor()
        text_actor.SetInput(text)
        text_actor.GetTextProperty().SetFontSize(16)
        text_actor.GetTextProperty().SetColor(1.0, 1.0, 1.0)
        text_actor.SetPosition(10, 10)
        
        self.renderer.AddActor2D(text_actor)
        
    def calculate_organoid_metrics(self, dataset):
        """Calculate morphology metrics"""
        points = dataset.GetPoints()
        if not points or points.GetNumberOfPoints() == 0:
            return None
            
        coords = []
        for i in range(points.GetNumberOfPoints()):
            coords.append(points.GetPoint(i))
        coords = np.array(coords)
        
        centroid = np.mean(coords, axis=0)
        distances = np.linalg.norm(coords - centroid, axis=1)
        
        return {
            'num_cells': len(coords),
            'centroid': centroid,
            'mean_radius': float(np.mean(distances)),
            'sphericity': max(0.0, 1.0 - (np.std(distances) / np.mean(distances)))
        }
        
    def visualize_organoid(self, pattern="results_*.vtu", surface_method="metaball", show_cells=False):
        """Visualize organoid with tissue-like appearance"""
        # Find VTU files
        search_path = os.path.join(self.data_dir, "**", pattern)
        vtu_files = glob.glob(search_path, recursive=True)
        
        if not vtu_files:
            print(f"No VTU files found in {self.data_dir}")
            return
            
        # Sort and get final time step
        def extract_time_step(filename):
            try:
                base = os.path.basename(filename)
                return int(base.split('_')[1].split('.')[0])
            except:
                return 0
                
        vtu_files.sort(key=extract_time_step)
        final_file = vtu_files[-1]
        
        print(f"Visualizing: {os.path.basename(final_file)}")
        
        # Load data
        dataset = self.load_vtk_file(final_file)
        metrics = self.calculate_organoid_metrics(dataset)
        
        if not metrics:
            print("No data to visualize")
            return
            
        # Create surface
        if surface_method == "metaball":
            surface = self.create_metaball_surface(dataset, radius=6.0)
        else:
            surface = self.create_tissue_surface(dataset)
            
        if surface:
            # Main organoid surface
            organoid_actor = self.create_organoid_actor(surface, "tissue")
            if organoid_actor:
                self.renderer.AddActor(organoid_actor)
                
            # Optional: show cell centers
            if show_cells:
                cell_actors = self.add_cell_centers(dataset, True)
                for actor in cell_actors:
                    self.renderer.AddActor(actor)
                    
            # Add coordinate axes
            self.add_coordinate_axes()
            
            # Add text overlay
            self.add_text_overlay(metrics, extract_time_step(final_file))
            
            # Set up camera
            camera = self.renderer.GetActiveCamera()
            camera.SetPosition(
                metrics['centroid'][0] + metrics['mean_radius'] * 2,
                metrics['centroid'][1] + metrics['mean_radius'] * 2,
                metrics['centroid'][2] + metrics['mean_radius'] * 1.5
            )
            camera.SetFocalPoint(metrics['centroid'])
            camera.SetViewUp(0, 0, 1)
            
            # Add lighting
            light = vtk.vtkLight()
            light.SetPosition(100, 100, 100)
            light.SetFocalPoint(0, 0, 0)
            self.renderer.AddLight(light)
            
            print(f"\nOrganoid Tissue Visualization:")
            print(f"  Cells: {metrics['num_cells']}")
            print(f"  Radius: {metrics['mean_radius']:.1f}μm")
            print(f"  Sphericity: {metrics['sphericity']:.3f}")
            
            # Render
            self.render_window.Render()
            
            print("\nControls:")
            print("  Rotate: Left mouse drag")
            print("  Zoom: Right mouse drag")
            print("  Close: Press 'q'")
            
            self.render_window_interactor.Start()
        else:
            print("Failed to create organoid surface")


def main():
    parser = argparse.ArgumentParser(description='Enhanced 3D Organoid Tissue Visualization')
    parser.add_argument('--dir', '-d', default='testoutput/Organoid3d', 
                       help='Base directory containing simulation results')
    parser.add_argument('--simulation', '-s', default='ThresholdECMForce',
                       help='Simulation to visualize')
    parser.add_argument('--method', '-m', choices=['metaball', 'hull'], default='metaball',
                       help='Surface creation method')
    parser.add_argument('--show-cells', action='store_true',
                       help='Show individual cell centers')
    
    args = parser.parse_args()
    
    sim_dir = os.path.join(args.dir, args.simulation)
    
    if not os.path.exists(sim_dir):
        print(f"Simulation directory not found: {sim_dir}")
        return 1
        
    visualizer = EnhancedOrganoid3DVisualizer(sim_dir)
    visualizer.visualize_organoid(
        surface_method=args.method,
        show_cells=args.show_cells
    )
    
    return 0


if __name__ == "__main__":
    exit(main())