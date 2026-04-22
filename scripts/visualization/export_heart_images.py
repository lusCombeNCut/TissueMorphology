#!/usr/bin/env python3
"""
Export images from the cardiac simulation for visualization
"""

import vtk
import os

def export_heart_images(vtu_file, output_dir="heart_images"):
    """Export images showing the electrical propagation"""
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Read the VTU file
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_file)
    reader.Update()
    
    data = reader.GetOutput()
    point_data = data.GetPointData()
    
    # Find all voltage arrays
    voltage_arrays = []
    for i in range(point_data.GetNumberOfArrays()):
        array_name = point_data.GetArrayName(i)
        if array_name.startswith("V_"):
            voltage_arrays.append(array_name)
    
    voltage_arrays.sort()
    print(f"Found {len(voltage_arrays)} time steps")
    
    # Set up visualization pipeline
    # Extract surface for better visualization
    surface = vtk.vtkDataSetSurfaceFilter()
    surface.SetInputData(data)
    surface.Update()
    
    # Create mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(surface.GetOutputPort())
    mapper.SetScalarModeToUsePointFieldData()
    
    # Create actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    
    # Create renderer
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(0.0, 0.0, 0.0)  # Black background
    
    # Create render window
    render_window = vtk.vtkRenderWindow()
    render_window.SetOffScreenRendering(1)  # Don't show window
    render_window.AddRenderer(renderer)
    render_window.SetSize(800, 600)
    
    # Set up camera
    renderer.ResetCamera()
    camera = renderer.GetActiveCamera()
    camera.Azimuth(45)
    camera.Elevation(30)
    camera.Zoom(1.2)
    
    # Export images for each time step
    for i, array_name in enumerate(voltage_arrays[::2]):  # Every other timestep to save space
        # Set the voltage array to visualize
        mapper.SelectColorArray(array_name)
        
        # Get voltage range for this timestep
        array = point_data.GetArray(array_name)
        voltage_range = array.GetRange()
        mapper.SetScalarRange(voltage_range)
        
        # Update color map - use a cardiac-appropriate colormap
        lut = vtk.vtkLookupTable()
        lut.SetRange(voltage_range)
        lut.SetHueRange(0.67, 0.0)  # Blue to red
        lut.Build()
        mapper.SetLookupTable(lut)
        
        # Render
        render_window.Render()
        
        # Save image
        writer = vtk.vtkPNGWriter()
        window_to_image = vtk.vtkWindowToImageFilter()
        window_to_image.SetInput(render_window)
        window_to_image.Update()
        
        output_file = os.path.join(output_dir, f"heart_voltage_{i:03d}.png")
        writer.SetFileName(output_file)
        writer.SetInputConnection(window_to_image.GetOutputPort())
        writer.Write()
        
        print(f"Exported: {output_file} (Voltage range: {voltage_range[0]:.1f} to {voltage_range[1]:.1f} mV)")
    
    print(f"\nImages exported to: {os.path.abspath(output_dir)}")
    print("You can view these to see the electrical wave propagation!")

if __name__ == "__main__":
    vtu_file = "/home/orlando/Thesis/testoutput/Monodomain3dRabbitHeart/vtk_output/results.vtu"
    export_heart_images(vtu_file)