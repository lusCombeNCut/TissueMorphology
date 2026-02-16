#!/usr/bin/env python3
"""
Plot cell invasion analysis for different ECM fiber orientations.

This script analyzes cell migration data from the DynamicECMInvasion2d simulations
and plots:
1. Invasion front position over time for each ECM type
2. Invasion speed comparison
3. Cell density profiles at different timepoints
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import xml.etree.ElementTree as ET
from typing import Dict, List, Tuple

# Font size constants for easy customization
FONT_SIZE_TITLE = 22
FONT_SIZE_LABEL = 18
FONT_SIZE_TICK = 18
FONT_SIZE_LEGEND = 18
FONT_SIZE_SUPTITLE = 22

# Try to import seaborn for styling, but continue without it
try:
    import seaborn as sns
    sns.set_style("whitegrid")
except ImportError:
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')

plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 11

def read_vtu_file(filepath: Path) -> Tuple[np.ndarray, np.ndarray]:
    """
    Read cell positions from a VTU file.
    
    Returns:
        positions: Nx2 array of cell positions
        cell_data: dictionary of cell data arrays
    """
    try:
        # Use VTK library if available, otherwise parse XML manually
        try:
            import vtk
            from vtk.util.numpy_support import vtk_to_numpy
            
            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(str(filepath))
            reader.Update()
            
            output = reader.GetOutput()
            points = output.GetPoints()
            
            if points:
                positions = vtk_to_numpy(points.GetData())[:, :2]  # Take x,y only
                return positions, {}
            else:
                return np.array([]), {}
                
        except ImportError:
            # Fallback: try to read ASCII format VTU
            tree = ET.parse(filepath)
            root = tree.getroot()
            
            # Find Points data
            points_data = root.find('.//Points/DataArray')
            if points_data is not None and points_data.text and points_data.get('format') == 'ascii':
                coords = np.fromstring(points_data.text, sep=' ')
                # Reshape to Nx3, then take first 2 columns for 2D
                positions = coords.reshape(-1, 3)[:, :2]
            else:
                # Try to get number of points from Piece element
                piece = root.find('.//Piece')
                if piece is not None:
                    num_points = int(piece.get('NumberOfPoints', 0))
                    if num_points > 0:
                        print(f"  Note: {filepath.name} has {num_points} points in binary format")
                        print(f"        Install python3-vtk to read binary VTU files")
                return np.array([]), {}
            
            return positions, {}
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return np.array([]), {}

def read_pvd_file(pvd_path: Path) -> List[Tuple[float, Path]]:
    """
    Read PVD file and return list of (time, vtu_filepath) tuples.
    """
    tree = ET.parse(pvd_path)
    root = tree.getroot()
    
    timesteps = []
    for dataset in root.findall('.//DataSet'):
        time = float(dataset.get('timestep'))
        filename = dataset.get('file')
        vtu_path = pvd_path.parent / filename
        timesteps.append((time, vtu_path))
    
    return sorted(timesteps, key=lambda x: x[0])

def calculate_invasion_metrics(positions: np.ndarray, domain_height: float = 1000.0) -> Dict[str, float]:
    """
    Calculate invasion metrics from cell positions.
    
    Returns:
        dict with:
        - max_y: maximum y position (invasion front)
        - mean_y: mean y position
        - std_y: standard deviation of y positions
        - num_cells: number of cells
    """
    if len(positions) == 0:
        return {
            'max_y': 0.0,
            'mean_y': 0.0,
            'std_y': 0.0,
            'num_cells': 0,
            'percentile_90': 0.0,
            'percentile_75': 0.0
        }
    
    y_positions = positions[:, 1]
    
    return {
        'max_y': np.max(y_positions),
        'mean_y': np.mean(y_positions),
        'std_y': np.std(y_positions),
        'num_cells': len(positions),
        'percentile_90': np.percentile(y_positions, 90),
        'percentile_75': np.percentile(y_positions, 75)
    }

def analyze_ecm_type(ecm_type: str, base_dir: Path) -> Dict[str, np.ndarray]:
    """
    Analyze invasion for a specific ECM type.
    
    Returns:
        dict with arrays of:
        - times
        - max_y (invasion front)
        - mean_y
        - num_cells
        - percentile_90
    """
    results_dir = base_dir / ecm_type
    
    # Look for results_from_time_0 subdirectory (common in Chaste output)
    if (results_dir / "results_from_time_0").exists():
        results_dir = results_dir / "results_from_time_0"
    
    pvd_file = results_dir / "results.pvd"
    
    if not pvd_file.exists():
        print(f"Warning: {pvd_file} not found")
        return None
    
    timesteps = read_pvd_file(pvd_file)
    
    times = []
    max_y_values = []
    mean_y_values = []
    num_cells_values = []
    percentile_90_values = []
    percentile_75_values = []
    
    for time, vtu_path in timesteps:
        if not vtu_path.exists():
            continue
        
        positions, _ = read_vtu_file(vtu_path)
        metrics = calculate_invasion_metrics(positions)
        
        times.append(time)
        max_y_values.append(metrics['max_y'])
        mean_y_values.append(metrics['mean_y'])
        num_cells_values.append(metrics['num_cells'])
        percentile_90_values.append(metrics['percentile_90'])
        percentile_75_values.append(metrics['percentile_75'])
    
    return {
        'times': np.array(times),
        'max_y': np.array(max_y_values),
        'mean_y': np.array(mean_y_values),
        'num_cells': np.array(num_cells_values),
        'percentile_90': np.array(percentile_90_values),
        'percentile_75': np.array(percentile_75_values)
    }

def calculate_invasion_speed(times: np.ndarray, positions: np.ndarray, 
                            window_hours: float = 3.0) -> np.ndarray:
    """
    Calculate invasion speed using a sliding window.
    
    Args:
        times: time array in minutes
        positions: position array (e.g., max_y)
        window_hours: window size in hours for speed calculation
    
    Returns:
        speeds: array of speeds in μm/hour
    """
    window_minutes = window_hours * 60
    speeds = np.zeros_like(times)
    
    for i in range(len(times)):
        # Find points within window before current time
        mask = (times <= times[i]) & (times >= times[i] - window_minutes)
        if np.sum(mask) >= 2:
            t_window = times[mask]
            p_window = positions[mask]
            
            # Linear fit
            if len(t_window) > 1:
                slope, _ = np.polyfit(t_window, p_window, 1)
                speeds[i] = slope * 60  # Convert from μm/min to μm/hour
    
    return speeds

def plot_invasion_comparison(results_dict: Dict[str, Dict], output_dir: Path):
    """
    Create invasion comparison plots (2 metrics, excluding mixed configuration).
    """
    ecm_types = ['perpendicular', 'parallel', 'random']
    colors = {
        'perpendicular': '#2ecc71',  # green
        'parallel': '#e74c3c',       # red
        'random': '#95a5a6',         # gray
    }
    labels = {
        'perpendicular': 'Perpendicular (Fast)',
        'parallel': 'Parallel (Slow)',
        'random': 'Random',
    }
    
    # Filter to only available results
    available_types = [t for t in ecm_types if t in results_dict and results_dict[t] is not None]
    
    if not available_types:
        print("No results available to plot")
        return
    
    # Create figure with 2 subplots in a row
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    
    # 1. Mean invasion depth (90th percentile)
    ax1 = axes[0]
    for ecm_type in available_types:
        data = results_dict[ecm_type]
        times_hours = data['times'] / 60
        ax1.plot(times_hours, data['percentile_90'], 
                label=labels[ecm_type], 
                color=colors[ecm_type], 
                linewidth=2.5)
    
    ax1.set_xlabel('Time (hours)', fontsize=FONT_SIZE_LABEL, fontweight='bold')
    ax1.set_ylabel('90th Percentile Depth (μm)', fontsize=FONT_SIZE_LABEL, fontweight='bold')
    ax1.set_title('Cell Invasion Front\n(90th Percentile)', fontsize=FONT_SIZE_TITLE, fontweight='bold')
    ax1.tick_params(labelsize=FONT_SIZE_TICK)
    ax1.grid(True, alpha=0.3)
    
    # 2. Invasion speed over time
    ax2 = axes[1]
    for ecm_type in available_types:
        data = results_dict[ecm_type]
        times_hours = data['times'] / 60
        speeds = calculate_invasion_speed(data['times'], data['max_y'])
        # Smooth speeds for visualization
        if len(speeds) > 5:
            from scipy.ndimage import gaussian_filter1d
            speeds_smooth = gaussian_filter1d(speeds, sigma=2)
        else:
            speeds_smooth = speeds
        ax2.plot(times_hours, speeds_smooth, 
                label=labels[ecm_type], 
                color=colors[ecm_type], 
                linewidth=2.5, 
                alpha=0.8)
    
    ax2.set_xlabel('Time (hours)', fontsize=FONT_SIZE_LABEL, fontweight='bold')
    ax2.set_ylabel('Invasion Speed (μm/hour)', fontsize=FONT_SIZE_LABEL, fontweight='bold')
    ax2.set_title('Invasion Speed Over Time\n(3h window)', fontsize=FONT_SIZE_TITLE, fontweight='bold')
    ax2.tick_params(labelsize=FONT_SIZE_TICK)
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    
    # Add single legend at the bottom center
    handles, labels_list = ax1.get_legend_handles_labels()
    fig.legend(handles, labels_list, loc='lower center', bbox_to_anchor=(0.5, 0.02), 
              ncol=3, framealpha=0.9, fontsize=FONT_SIZE_LEGEND)
    
    # plt.suptitle('ECM Fiber Orientation Effects on Cell Invasion', 
    #             fontsize=16, fontweight='bold', y=1.02)
    
    plt.tight_layout(rect=[0, 0.15, 1, 1])
    
    # Save figure
    output_file = output_dir / 'invasion_analysis.pdf'
    plt.savefig(output_file, dpi=300, format='pdf', bbox_inches='tight')
    print(f"Saved: {output_file}")
    
    plt.close()

def create_invasion_summary_table(results_dict: Dict[str, Dict], output_dir: Path):
    """
    Create a summary table of invasion metrics.
    """
    summary_file = output_dir / 'invasion_summary.txt'
    
    with open(summary_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("INVASION ANALYSIS SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        
        for ecm_type, data in results_dict.items():
            if data is None:
                continue
            
            f.write(f"\n{ecm_type.upper()} ECM\n")
            f.write("-" * 40 + "\n")
            
            if len(data['max_y']) > 0:
                final_time = data['times'][-1] / 60  # hours
                final_invasion = data['max_y'][-1]
                final_cells = data['num_cells'][-1]
                
                # Calculate average speed
                avg_speed = final_invasion / final_time if final_time > 0 else 0
                
                f.write(f"Final time: {final_time:.1f} hours\n")
                f.write(f"Final invasion depth: {final_invasion:.2f} μm\n")
                f.write(f"Average invasion speed: {avg_speed:.2f} μm/hour\n")
                f.write(f"Final cell count: {final_cells}\n")
                f.write(f"90th percentile depth: {data['percentile_90'][-1]:.2f} μm\n")
                f.write(f"75th percentile depth: {data['percentile_75'][-1]:.2f} μm\n")
            else:
                f.write("No data available\n")
    
    print(f"Saved: {summary_file}")

def main():
    """
    Main analysis function.
    """
    # Define paths - go up to Thesis directory
    script_dir = Path(__file__).parent.resolve()
    # Assuming script is in Chaste/projects/TissueMorphology
    thesis_dir = script_dir.parent.parent.parent
    base_dir = thesis_dir / "testoutput" / "DynamicECMInvasion2d"
    output_dir = script_dir / "invasion_analysis_output"
    output_dir.mkdir(exist_ok=True)
    
    print("Analyzing ECM invasion simulations...")
    print(f"Script directory: {script_dir}")
    print(f"Thesis directory: {thesis_dir}")
    print(f"Looking for data in: {base_dir}")
    print(f"Base dir exists: {base_dir.exists()}")
    
    # Analyze each ECM type
    results = {}
    ecm_types = ['perpendicular', 'parallel', 'random', 'mixed']
    
    for ecm_type in ecm_types:
        print(f"\nAnalyzing {ecm_type} ECM...")
        results[ecm_type] = analyze_ecm_type(ecm_type, base_dir)
        
        if results[ecm_type] is not None:
            n_timepoints = len(results[ecm_type]['times'])
            print(f"  Found {n_timepoints} timepoints")
            if n_timepoints > 0:
                print(f"  Time range: {results[ecm_type]['times'][0]/60:.1f} - {results[ecm_type]['times'][-1]/60:.1f} hours")
                print(f"  Final invasion: {results[ecm_type]['max_y'][-1]:.2f} μm")
    
    # Create plots
    print("\nGenerating plots...")
    plot_invasion_comparison(results, output_dir)
    
    # Create summary table
    create_invasion_summary_table(results, output_dir)
    
    print("\n✓ Analysis complete!")
    print(f"Output saved to: {output_dir}")

if __name__ == '__main__':
    # Try to import scipy for smoothing, but continue if not available
    try:
        from scipy.ndimage import gaussian_filter1d
    except ImportError:
        print("Warning: scipy not available, speed plots won't be smoothed")
        # Define a dummy function
        def gaussian_filter1d(x, sigma):
            return x
    
    main()
