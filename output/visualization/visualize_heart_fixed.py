#!/usr/bin/env python3
"""
Fixed interactive visualization that handles X11 threading issues
"""

import vtk
import os
import sys

def fix_x11_threading():
    """Fix X11 threading issues"""
    os.environ['QT_X11_NO_MITSHM'] = '1'
    os.environ['LIBGL_ALWAYS_SOFTWARE'] = '1'

def create_robust_heart_visualization():
    """Create a robust interactive heart visualization with animation"""
    
    # Fix threading issues first
    fix_x11_threading()
    
    vtu_file = "/home/orlando/Thesis/testoutput/Monodomain3dRabbitHeart/vtk_output/results.vtu"
    
    print("Loading cardiac simulation data...")
    
    # Read the VTU file
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_file)
    reader.Update()
    
    data = reader.GetOutput()
    point_data = data.GetPointData()
    
    # Find voltage arrays
    voltage_arrays = []
    for i in range(point_data.GetNumberOfArrays()):
        array = point_data.GetArray(i)
        name = array.GetName()
        if name and ('V_' in name):
            voltage_arrays.append(name)
    
    voltage_arrays.sort()
    print(f"Found {len(voltage_arrays)} voltage time steps")
    
    # Analyze voltage ranges for each time step
    print("\\nVoltage analysis by time step:")
    for i, array_name in enumerate(voltage_arrays[:6]):  # Show first 6
        array = point_data.GetArray(array_name)
        v_range = array.GetRange()
        time_ms = i * 0.2
        print(f"  {array_name} (t={time_ms:.1f}ms): {v_range[0]:.1f} to {v_range[1]:.1f} mV (span: {v_range[1]-v_range[0]:.1f})")
    
    # Use a time step with good variation (V_000006 showed best results)
    if "V_000006" in voltage_arrays:
        active_array = "V_000006"
    elif len(voltage_arrays) > 5:
        active_array = voltage_arrays[5]
    else:
        active_array = voltage_arrays[-1]
    
    print(f"\\nUsing voltage array: {active_array} for visualization")
    
    # Set active scalars
    data.GetPointData().SetActiveScalars(active_array)
    voltage_array = data.GetPointData().GetArray(active_array)
    voltage_range = voltage_array.GetRange()
    
    print(f"Voltage range: {voltage_range[0]:.1f} to {voltage_range[1]:.1f} mV")
    
    # Convert to numpy to check actual variation
    import numpy as np
    numpy_voltages = np.zeros(voltage_array.GetNumberOfTuples())
    for j in range(voltage_array.GetNumberOfTuples()):
        numpy_voltages[j] = voltage_array.GetValue(j)
    
    unique_vals = len(np.unique(numpy_voltages))
    active_cells = np.sum(numpy_voltages > -80)
    print(f"Unique voltage values: {unique_vals:,}")
    print(f"Active cells (>-80mV): {active_cells:,} ({100*active_cells/len(numpy_voltages):.2f}%)")
    
    # Create surface filter for better performance
    surface_filter = vtk.vtkDataSetSurfaceFilter()
    surface_filter.SetInputData(data)
    surface_filter.Update()
    
    # Create optimized color lookup table
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(256)
    
    # Use the full voltage range but optimize for cardiac data
    min_voltage = voltage_range[0]
    max_voltage = voltage_range[1]
    lut.SetRange(min_voltage, max_voltage)
    
    print(f"Color mapping range: {min_voltage:.1f} to {max_voltage:.1f} mV")
    
    # Create VERY HIGH CONTRAST cardiac color scheme
    for i in range(256):
        val = i / 255.0
        voltage = min_voltage + val * (max_voltage - min_voltage)
        
        # EXTREME contrast color mapping - make small voltage changes very visible
        if voltage < -85:  # Deep resting - dark blue
            r, g, b = 0.0, 0.0, 0.3
        elif voltage < -84:  # Slightly less resting - blue
            r, g, b = 0.0, 0.0, 1.0
        elif voltage < -83:  # Light blue
            r, g, b = 0.0, 0.5, 1.0
        elif voltage < -82:  # Cyan
            r, g, b = 0.0, 1.0, 1.0
        elif voltage < -80:  # Light cyan
            r, g, b = 0.0, 1.0, 0.5
        elif voltage < -70:  # Green
            r, g, b = 0.0, 1.0, 0.0
        elif voltage < -50:  # Yellow-green
            r, g, b = 0.5, 1.0, 0.0
        elif voltage < -20:  # Yellow
            r, g, b = 1.0, 1.0, 0.0
        elif voltage < 0:    # Orange
            r, g, b = 1.0, 0.5, 0.0
        elif voltage < 20:   # Red-orange
            r, g, b = 1.0, 0.2, 0.0
        else:               # Bright red for any activation
            r, g, b = 1.0, 0.0, 0.0
        
        lut.SetTableValue(i, r, g, b, 1.0)
    
    lut.Build()
    
    # Create mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(surface_filter.GetOutputPort())
    mapper.SetScalarRange(min_voltage, max_voltage)
    mapper.SetLookupTable(lut)
    mapper.SetScalarModeToUsePointData()
    
    # Create actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    
    # Create renderer with better settings
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(0.0, 0.0, 0.0)  # Black background
    
    # Add lighting for better 3D appearance
    light = vtk.vtkLight()
    light.SetPosition(1, 1, 1)
    light.SetIntensity(0.8)
    renderer.AddLight(light)
    
    # Create render window with better settings
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(1400, 1000)
    render_window.SetWindowName("Rabbit Heart - Interactive Voltage Visualization")
    
    # Disable some problematic features
    render_window.SetMultiSamples(0)  # Disable antialiasing
    
    # Add color bar
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(lut)
    scalar_bar.SetTitle("Voltage (mV)")
    scalar_bar.SetNumberOfLabels(8)
    scalar_bar.SetPosition(0.87, 0.1)
    scalar_bar.SetWidth(0.1)
    scalar_bar.SetHeight(0.8)
    
    # Style the color bar
    scalar_bar.GetLabelTextProperty().SetColor(1, 1, 1)
    scalar_bar.GetTitleTextProperty().SetColor(1, 1, 1)
    scalar_bar.GetLabelTextProperty().SetFontSize(12)
    scalar_bar.GetTitleTextProperty().SetFontSize(14)
    
    renderer.AddActor2D(scalar_bar)
    
    # Add informative text
    info_text = vtk.vtkTextActor()
    info_text.SetInput(f"Rabbit Heart Electrophysiology\\n" +
                      f"Time: {active_array}\\n" +
                      f"Voltage Range: {min_voltage:.0f} to {max_voltage:.0f} mV\\n" +
                      f"Blue=Resting, Red=Active\\n" +
                      f"Controls: Left=Rotate, Right=Zoom, 'q'=Quit")
    info_text.SetPosition(10, 10)
    info_text.GetTextProperty().SetColor(1, 1, 1)
    info_text.GetTextProperty().SetFontSize(14)
    info_text.GetTextProperty().SetBackgroundColor(0, 0, 0)
    info_text.GetTextProperty().SetBackgroundOpacity(0.8)
    renderer.AddActor2D(info_text)
    
    # Set up camera for optimal heart viewing
    camera = renderer.GetActiveCamera()
    bounds = data.GetBounds()
    center = [(bounds[0] + bounds[1]) / 2,
              (bounds[2] + bounds[3]) / 2, 
              (bounds[4] + bounds[5]) / 2]
    
    # Position camera for good heart view (slightly from below to see apex)
    distance = max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4]) * 2
    camera.SetPosition(center[0] + distance*0.3, 
                      center[1] - distance*0.8, 
                      center[2] + distance*0.2)
    camera.SetFocalPoint(center[0], center[1], center[2])
    camera.SetViewUp(0, 0, 1)
    renderer.ResetCamera()
    
    # Create custom interaction style with time step control
    class HeartInteractionStyle(vtk.vtkInteractorStyleTrackballCamera):
        def __init__(self):
            super().__init__()
            self.current_time_step = voltage_arrays.index(active_array)
            self.voltage_arrays = voltage_arrays
            self.data = data
            self.mapper = mapper
            self.info_text = info_text
            
        def OnKeyPress(self):
            key = self.GetInteractor().GetKeySym()
            
            if key == 'n' or key == 'Right':  # Next time step
                self.current_time_step = (self.current_time_step + 1) % len(self.voltage_arrays)
                self.update_time_step()
            elif key == 'p' or key == 'Left':  # Previous time step  
                self.current_time_step = (self.current_time_step - 1) % len(self.voltage_arrays)
                self.update_time_step()
            elif key == 'space':  # Auto-animate
                self.animate()
            
            # Call parent
            super().OnKeyPress()
            
        def update_time_step(self):
            array_name = self.voltage_arrays[self.current_time_step]
            self.data.GetPointData().SetActiveScalars(array_name)
            
            # Update info text
            time_ms = self.current_time_step * 0.2
            v_array = self.data.GetPointData().GetArray(array_name)
            v_range = v_array.GetRange()
            
            self.info_text.SetInput(f"Rabbit Heart Electrophysiology\\n" +
                                  f"Time: {time_ms:.1f}ms ({array_name})\\n" +
                                  f"Voltage: {v_range[0]:.0f} to {v_range[1]:.0f} mV\\n" +
                                  f"Controls: n/p=time, space=animate, q=quit")
            
            # Force render
            self.GetInteractor().GetRenderWindow().Render()
            print(f"Switched to {array_name} (t={time_ms:.1f}ms)")
            
        def animate(self):
            """Auto-animate through all time steps"""
            print("Starting animation...")
            for i in range(len(self.voltage_arrays)):
                self.current_time_step = i
                self.update_time_step()
                import time
                time.sleep(0.5)  # 0.5 second per frame
    
    # Create interactor with custom style
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)
    
    # Set custom interaction style
    style = HeartInteractionStyle()
    interactor.SetInteractorStyle(style)
    
    print("\\n=== Interactive Cardiac Visualization ===")
    print("Features:")
    print("- 3D rabbit heart with HIGH CONTRAST voltage coloring")
    print("- Real cardiac simulation data with time animation")
    print("- Enhanced visualization for small voltage changes")
    print("\\nControls:")
    print("- Mouse drag: Rotate/zoom/pan the heart")
    print("- 'n' or Right Arrow: Next time step")
    print("- 'p' or Left Arrow: Previous time step")  
    print("- Space bar: Auto-animate through all time steps")
    print("- 'r': Reset camera view")
    print("- 'q': Quit visualization")
    print("\\nStarting visualization...")
    
    try:
        # Initialize and start
        interactor.Initialize()
        render_window.Render()
        interactor.Start()
        
    except Exception as e:
        print(f"Interactive mode failed: {e}")
        print("\\nFalling back to image export mode...")
        
        # Export a high-quality image instead
        render_window.SetOffScreenRendering(1)
        render_window.Render()
        
        # Save the image
        writer = vtk.vtkPNGWriter()
        window_to_image = vtk.vtkWindowToImageFilter()
        window_to_image.SetInput(render_window)
        window_to_image.Update()
        
        output_file = "/home/orlando/Thesis/heart_interactive_fallback.png"
        writer.SetFileName(output_file)
        writer.SetInputConnection(window_to_image.GetOutputPort())
        writer.Write()
        
        print(f"Image saved to: {output_file}")
        print("You can view this image to see the heart visualization")

if __name__ == "__main__":
    create_robust_heart_visualization()