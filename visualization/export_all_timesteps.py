#!/usr/bin/env python3
"""
Export images of all time steps for offline viewing
"""

import vtk
import os

def export_all_timesteps():
    """Export images of all time steps"""
    
    vtu_file = "/home/orlando/Thesis/testoutput/Monodomain3dRabbitHeart/vtk_output/results.vtu"
    output_dir = "/home/orlando/Thesis/heart_timesteps"
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # All time steps
    time_steps = [f"V_{i:06d}" for i in range(11)]
    time_values = [i * 0.2 for i in range(11)]
    
    print(f"Exporting {len(time_steps)} time steps to {output_dir}")
    
    # Create shared lookup table
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(256)
    min_range = -100
    max_range = 110
    lut.SetRange(min_range, max_range)
    
    # Cardiac color scheme
    for j in range(256):
        val = j / 255.0
        voltage = min_range + val * (max_range - min_range)
        
        if voltage < -85:
            r, g, b = 0.0, 0.0, 0.4
        elif voltage < -70:
            t = (voltage + 85) / 15
            r, g, b = 0.0, 0.2*t, 0.4 + 0.4*t
        elif voltage < -50:
            t = (voltage + 70) / 20
            r, g, b = 0.0, 0.2 + 0.6*t, 0.8 + 0.2*t
        elif voltage < -20:
            t = (voltage + 50) / 30
            r, g, b = 0.0, 0.8 + 0.2*t, 1.0
        elif voltage < 0:
            t = (voltage + 20) / 20
            r, g, b = t*0.6, 1.0, 1.0 - t*0.4
        elif voltage < 30:
            t = voltage / 30
            r, g, b = 0.6 + 0.4*t, 1.0 - 0.3*t, 0.6*(1-t)
        elif voltage < 60:
            t = (voltage - 30) / 30
            r, g, b = 1.0, 0.7 - 0.4*t, 0.6*t
        else:
            t = min(1.0, (voltage - 60) / 50)
            r, g, b = 1.0, 0.3*(1-t), 0.6 + 0.4*t
        
        lut.SetTableValue(j, r, g, b, 1.0)
    
    lut.Build()
    
    for i, (time_step, time_ms) in enumerate(zip(time_steps, time_values)):
        print(f"Processing {time_step} (t={time_ms:.1f} ms)...")
        
        # Read data
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(vtu_file)
        reader.Update()
        
        data = reader.GetOutput()
        
        if not data.GetPointData().GetArray(time_step):
            print(f"  Warning: {time_step} not found, skipping...")
            continue
        
        data.GetPointData().SetActiveScalars(time_step)
        voltage_array = data.GetPointData().GetArray(time_step)
        v_range = voltage_array.GetRange()
        
        # Create surface filter for better visualization
        surface = vtk.vtkDataSetSurfaceFilter()
        surface.SetInputData(data)
        surface.Update()
        
        # Create mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(surface.GetOutputPort())
        mapper.SetScalarRange(min_range, max_range)
        mapper.SetLookupTable(lut)
        
        # Create actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        
        # Create renderer
        renderer = vtk.vtkRenderer()
        renderer.AddActor(actor)
        renderer.SetBackground(0.0, 0.0, 0.0)
        
        # Add color bar
        scalar_bar = vtk.vtkScalarBarActor()
        scalar_bar.SetLookupTable(lut)
        scalar_bar.SetTitle("Voltage (mV)")
        scalar_bar.SetNumberOfLabels(6)
        scalar_bar.SetPosition(0.85, 0.15)
        scalar_bar.SetWidth(0.12)
        scalar_bar.SetHeight(0.7)
        scalar_bar.GetLabelTextProperty().SetColor(1, 1, 1)
        scalar_bar.GetTitleTextProperty().SetColor(1, 1, 1)
        renderer.AddActor2D(scalar_bar)
        
        # Add title
        title_text = vtk.vtkTextActor()
        title_text.SetInput(f"Rabbit Heart Electrophysiology\\n" +
                           f"Time: {time_ms:.1f} ms ({time_step})\\n" +
                           f"Voltage Range: {v_range[0]:.1f} to {v_range[1]:.1f} mV")
        title_text.SetPosition(10, 10)
        title_text.GetTextProperty().SetColor(1, 1, 1)
        title_text.GetTextProperty().SetFontSize(14)
        title_text.GetTextProperty().SetBackgroundColor(0, 0, 0)
        title_text.GetTextProperty().SetBackgroundOpacity(0.8)
        renderer.AddActor2D(title_text)
        
        # Create render window (off-screen)
        render_window = vtk.vtkRenderWindow()
        render_window.SetOffScreenRendering(1)
        render_window.AddRenderer(renderer)
        render_window.SetSize(1200, 900)
        
        # Position camera
        camera = renderer.GetActiveCamera()
        bounds = data.GetBounds()
        center_x = (bounds[0] + bounds[1]) / 2
        center_y = (bounds[2] + bounds[3]) / 2
        center_z = (bounds[4] + bounds[5]) / 2
        
        distance = max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4]) * 1.5
        camera.SetPosition(center_x + distance*0.4, center_y - distance*0.7, center_z + distance*0.3)
        camera.SetFocalPoint(center_x, center_y, center_z)
        camera.SetViewUp(0, 0, 1)
        renderer.ResetCamera()
        
        # Render and save
        render_window.Render()
        
        # Save image
        writer = vtk.vtkPNGWriter()
        window_to_image = vtk.vtkWindowToImageFilter()
        window_to_image.SetInput(render_window)
        window_to_image.Update()
        
        output_file = os.path.join(output_dir, f"heart_timestep_{i:02d}_{time_ms:03.1f}ms.png")
        writer.SetFileName(output_file)
        writer.SetInputConnection(window_to_image.GetOutputPort())
        writer.Write()
        
        print(f"  Exported: {output_file}")
        print(f"    Voltage range: {v_range[0]:.1f} to {v_range[1]:.1f} mV")
    
    print(f"\\nAll images exported to: {output_dir}")
    print("You can now view these images to see the electrical wave propagation!")
    
    # Create a summary file
    summary_file = os.path.join(output_dir, "README.txt")
    with open(summary_file, 'w') as f:
        f.write("Rabbit Heart Electrical Wave Propagation\\n")
        f.write("=======================================\\n\\n")
        f.write("These images show the electrical activity in a 3D rabbit heart\\n")
        f.write("over time from 0.0 to 2.0 milliseconds.\\n\\n")
        f.write("Color Mapping:\\n")
        f.write("- Deep Blue: Resting potential (~-85 mV)\\n")
        f.write("- Light Blue/Cyan: Threshold crossing (-70 to -50 mV)\\n")
        f.write("- Green: Early depolarization (-50 to -20 mV)\\n")
        f.write("- Yellow: Mid depolarization (-20 to 0 mV)\\n")
        f.write("- Orange: Late depolarization (0 to +30 mV)\\n")
        f.write("- Red: Peak activation (+30 to +100 mV)\\n\\n")
        f.write("The electrical stimulus is applied at the apex (bottom tip)\\n")
        f.write("of the heart and propagates upward toward the base.\\n\\n")
        f.write("Files:\\n")
        for i, (time_step, time_ms) in enumerate(zip(time_steps, time_values)):
            f.write(f"  heart_timestep_{i:02d}_{time_ms:03.1f}ms.png - Time {time_ms:.1f} ms\\n")
    
    print(f"Summary written to: {summary_file}")

if __name__ == "__main__":
    export_all_timesteps()