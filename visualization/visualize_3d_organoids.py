#!/usr/bin/env python3
"""
3D Organoid Visualization for Folding Analysis
Inspired by cardiac VTK visualization approaches
"""

import vtk
import numpy as np
import os
import sys
import argparse

def create_3d_organoid_visualization(results_directory, time_step=None):
    """
    Create 3D visualization of organoid formation with folding analysis
    
    Args:
        results_directory: Path to simulation output directory
        time_step: Specific time step to visualize (None for latest)
    """
    
    # Find VTK files in results directory
    vtu_files = []
    for root, dirs, files in os.walk(results_directory):
        for file in files:
            if file.endswith('.vtu'):
                vtu_files.append(os.path.join(root, file))
    
    if not vtu_files:
        print(f"No VTU files found in {results_directory}")
        return
    
    vtu_files.sort()
    print(f"Found {len(vtu_files)} VTU files")
    
    # Select file to visualize
    if time_step is not None and time_step < len(vtu_files):
        selected_file = vtu_files[time_step]
    else:
        selected_file = vtu_files[-1]  # Latest time step
    
    print(f"Visualizing: {selected_file}")
    
    # Read VTU file
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(selected_file)
    reader.Update()
    
    data = reader.GetOutput()
    print(f"Loaded {data.GetNumberOfPoints()} points, {data.GetNumberOfCells()} cells")
    
    # Analyze available data arrays
    point_data = data.GetPointData()
    print(f"Available point data arrays:")
    for i in range(point_data.GetNumberOfArrays()):
        array = point_data.GetArray(i)
        array_name = array.GetName()
        array_range = array.GetRange()
        print(f"  {array_name}: {array_range[0]:.3f} to {array_range[1]:.3f}")
    
    # Use cell age for coloring (or cell ID as fallback)
    color_array = "Cell ages"
    if not point_data.GetArray(color_array):
        color_array = "Cell IDs"
    if not point_data.GetArray(color_array):
        color_array = point_data.GetArrayName(0)  # Use first available
    
    data.GetPointData().SetActiveScalars(color_array)
    print(f"Using '{color_array}' for coloring")
    
    # Create organoid-specific lookup table
    color_range = data.GetPointData().GetArray(color_array).GetRange()
    
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(256)
    lut.SetRange(color_range)
    
    # Organoid development color scheme
    for i in range(256):
        val = i / 255.0
        
        if "age" in color_array.lower():
            # Age-based coloring: young (blue) to old (red)
            if val < 0.3:  # Young cells - blue to cyan
                r, g, b = 0.0, val/0.3, 1.0
            elif val < 0.7:  # Mature cells - cyan to yellow
                t = (val - 0.3) / 0.4
                r, g, b = t, 1.0, 1.0 - t
            else:  # Old cells - yellow to red
                t = (val - 0.7) / 0.3
                r, g, b = 1.0, 1.0 - t, 0.0
        else:
            # Generic coloring - rainbow
            import colorsys
            r, g, b = colorsys.hsv_to_rgb(val * 0.8, 0.8, 1.0)
        
        lut.SetTableValue(i, r, g, b, 1.0)
    
    lut.Build()
    
    # Create geometry filter for better rendering
    geometry_filter = vtk.vtkGeometryFilter()
    geometry_filter.SetInputData(data)
    geometry_filter.Update()
    
    # Create mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(geometry_filter.GetOutputPort())
    mapper.SetScalarRange(color_range)
    mapper.SetLookupTable(lut)
    mapper.SetScalarModeToUsePointData()
    
    # Create actor with organoid-appropriate appearance
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetSpecular(0.3)
    actor.GetProperty().SetSpecularPower(30)
    
    # Create renderer
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(0.1, 0.1, 0.2)  # Dark blue background
    
    # Add lighting for 3D depth perception
    light1 = vtk.vtkLight()
    light1.SetPosition(1, 1, 1)
    light1.SetColor(1, 1, 1)
    renderer.AddLight(light1)
    
    light2 = vtk.vtkLight()
    light2.SetPosition(-1, -1, 0.5)
    light2.SetColor(0.5, 0.5, 0.8)
    renderer.AddLight(light2)
    
    # Calculate organoid center and size for camera positioning
    bounds = data.GetBounds()
    center_x = (bounds[0] + bounds[1]) / 2
    center_y = (bounds[2] + bounds[3]) / 2  
    center_z = (bounds[4] + bounds[5]) / 2
    
    size_x = bounds[1] - bounds[0]
    size_y = bounds[3] - bounds[2]
    size_z = bounds[5] - bounds[4]
    max_size = max(size_x, size_y, size_z)
    
    print(f"Organoid center: ({center_x:.1f}, {center_y:.1f}, {center_z:.1f})")
    print(f"Organoid size: {size_x:.1f} x {size_y:.1f} x {size_z:.1f}")
    
    # Position camera for optimal 3D organoid view
    camera = renderer.GetActiveCamera()
    camera_distance = max_size * 3
    camera.SetPosition(center_x + camera_distance*0.6, 
                      center_y - camera_distance*0.8, 
                      center_z + camera_distance*0.4)
    camera.SetFocalPoint(center_x, center_y, center_z)
    camera.SetViewUp(0, 0, 1)
    renderer.ResetCamera()
    
    # Create render window
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(1200, 900)
    render_window.SetWindowName("3D Organoid Formation - Folding Analysis")
    
    # Add color bar
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(lut)
    scalar_bar.SetTitle(color_array)
    scalar_bar.SetNumberOfLabels(6)
    scalar_bar.SetPosition(0.85, 0.15)
    scalar_bar.SetWidth(0.1)
    scalar_bar.SetHeight(0.7)
    scalar_bar.GetLabelTextProperty().SetColor(1, 1, 1)
    scalar_bar.GetTitleTextProperty().SetColor(1, 1, 1)
    renderer.AddActor2D(scalar_bar)
    
    # Add informative text overlay
    info_text = vtk.vtkTextActor()
    file_name = os.path.basename(selected_file)
    info_text.SetInput(f"3D Organoid Formation\\n"
                      f"File: {file_name}\\n"
                      f"Points: {data.GetNumberOfPoints():,}\\n"
                      f"Coloring: {color_array}\\n"
                      f"Size: {max_size:.1f} units")
    info_text.SetPosition(10, 10)
    info_text.GetTextProperty().SetColor(1, 1, 1)
    info_text.GetTextProperty().SetFontSize(12)
    info_text.GetTextProperty().SetBackgroundColor(0, 0, 0)
    info_text.GetTextProperty().SetBackgroundOpacity(0.7)
    renderer.AddActor2D(info_text)
    
    # Add analysis text for folding studies
    analysis_text = vtk.vtkTextActor()
    analysis_text.SetInput("Folding Analysis:\\n"
                          "- Observe surface irregularities\\n"
                          "- Look for invaginations\\n"
                          "- Check cell distribution patterns\\n"
                          "- Compare with experimental data")
    analysis_text.SetPosition(10, 120)
    analysis_text.GetTextProperty().SetColor(1, 1, 0.5)
    analysis_text.GetTextProperty().SetFontSize(10)
    analysis_text.GetTextProperty().SetBackgroundColor(0, 0, 0)
    analysis_text.GetTextProperty().SetBackgroundOpacity(0.5)
    renderer.AddActor2D(analysis_text)
    
    # Create interactor with enhanced controls
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)
    
    # Add mouse interaction style for 3D navigation
    style = vtk.vtkInteractorStyleTrackballCamera()
    interactor.SetInteractorStyle(style)
    
    print(f"\\n=== 3D Organoid Visualization Ready ===")
    print(f"Controls:")
    print(f"  Left mouse + drag: Rotate view")
    print(f"  Right mouse + drag: Zoom in/out") 
    print(f"  Middle mouse + drag: Pan view")
    print(f"  'r': Reset camera view")
    print(f"  'w': Wireframe mode")
    print(f"  's': Surface mode")
    print(f"  'q': Quit")
    print(f"\\nLook for folding patterns, surface topology, and cell organization!")
    
    # Start visualization
    render_window.Render()
    interactor.Start()

def create_organoid_time_series(results_directory, output_images=False):
    """
    Create time series visualization of organoid development
    
    Args:
        results_directory: Path to simulation output directory
        output_images: Whether to save images for each time step
    """
    
    # Find all VTU files
    vtu_files = []
    for root, dirs, files in os.walk(results_directory):
        for file in files:
            if file.endswith('.vtu'):
                vtu_files.append(os.path.join(root, file))
    
    vtu_files.sort()
    
    if len(vtu_files) < 2:
        print("Need at least 2 time steps for time series analysis")
        return
    
    print(f"Creating time series from {len(vtu_files)} time steps")
    
    # Process each time step
    for i, vtu_file in enumerate(vtu_files):
        print(f"Processing time step {i+1}/{len(vtu_files)}: {os.path.basename(vtu_file)}")
        
        # Read data
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(vtu_file)
        reader.Update()
        data = reader.GetOutput()
        
        # Calculate organoid metrics
        bounds = data.GetBounds()
        size_x = bounds[1] - bounds[0]
        size_y = bounds[3] - bounds[2]
        size_z = bounds[5] - bounds[4]
        
        volume_estimate = size_x * size_y * size_z
        surface_area_estimate = 2 * (size_x*size_y + size_y*size_z + size_x*size_z)
        sphericity_estimate = (volume_estimate**(2/3)) / surface_area_estimate if surface_area_estimate > 0 else 0
        
        print(f"  Size: {size_x:.1f} x {size_y:.1f} x {size_z:.1f}")
        print(f"  Volume estimate: {volume_estimate:.1f}")
        print(f"  Sphericity estimate: {sphericity_estimate:.3f}")
        
        # Optionally save images
        if output_images:
            # Create a quick visualization and save
            pass  # Implementation for batch image saving

def main():
    parser = argparse.ArgumentParser(description='3D Organoid Visualization for Folding Analysis')
    parser.add_argument('results_dir', help='Path to simulation results directory')
    parser.add_argument('--time-step', type=int, help='Specific time step to visualize')
    parser.add_argument('--time-series', action='store_true', help='Analyze time series')
    parser.add_argument('--save-images', action='store_true', help='Save images for time series')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.results_dir):
        print(f"Error: Directory {args.results_dir} not found")
        return 1
    
    if args.time_series:
        create_organoid_time_series(args.results_dir, args.save_images)
    else:
        create_3d_organoid_visualization(args.results_dir, args.time_step)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())