#!/usr/bin/env python3
"""
Plot averaged cell invasion analysis for different ECM fiber orientations.

This script analyzes cell migration data from multiple runs of DynamicECMInvasion2d
simulations and plots averages with standard deviations:
1. 90th percentile invasion depth over time
2. Invasion speed comparison
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import xml.etree.ElementTree as ET
from typing import Dict, List, Tuple
from scipy import interpolate
from scipy.ndimage import gaussian_filter1d

# Font size constants matching export_ecm_snapshots.py
FONT_SIZE_TITLE = 18
FONT_SIZE_LABEL = 16
FONT_SIZE_TICK = 16
FONT_SIZE_LEGEND = 16
FONT_SIZE_SUPTITLE = 18

# Try to import seaborn for styling, but continue without it
try:
    import seaborn as sns
    sns.set_style("whitegrid")
except ImportError:
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')

plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 11

def read_vtu_file(filepath: Path) -> np.ndarray:
    """Read cell positions from a VTU file."""
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
            return positions
        else:
            return np.array([])
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return np.array([])

def read_pvd_file(pvd_path: Path) -> List[Tuple[float, Path]]:
    """Read PVD file and return list of (time, vtu_filepath) tuples."""
    try:
        tree = ET.parse(pvd_path)
        root = tree.getroot()
        
        timesteps = []
        for dataset in root.findall('.//DataSet'):
            time = float(dataset.get('timestep'))
            filename = dataset.get('file')
            vtu_path = pvd_path.parent / filename
            timesteps.append((time, vtu_path))
        
        return sorted(timesteps, key=lambda x: x[0])
    except Exception as e:
        print(f"  Warning: Could not parse {pvd_path.name}: {e}")
        return []

def analyze_single_run(run_dir: Path) -> Dict[str, np.ndarray]:
    """
    Analyze invasion for a single run.
    
    Returns dict with arrays: times, percentile_90, num_cells
    """
    # Look for results_from_time_0 subdirectory
    if (run_dir / "results_from_time_0").exists():
        results_dir = run_dir / "results_from_time_0"
    else:
        results_dir = run_dir
    
    pvd_file = results_dir / "results.pvd"
    
    if not pvd_file.exists():
        return None
    
    timesteps = read_pvd_file(pvd_file)
    
    if not timesteps:
        return None
    
    times = []
    percentile_90_values = []
    num_cells_values = []
    
    for time, vtu_path in timesteps:
        if not vtu_path.exists():
            continue
        
        positions = read_vtu_file(vtu_path)
        
        if len(positions) == 0:
            continue
        
        y_positions = positions[:, 1]
        
        times.append(time)
        percentile_90_values.append(np.percentile(y_positions, 90))
        num_cells_values.append(len(positions))
    
    return {
        'times': np.array(times),
        'percentile_90': np.array(percentile_90_values),
        'num_cells': np.array(num_cells_values)
    }

def analyze_ecm_type_multiple_runs(ecm_type: str, base_dir: Path) -> Dict[str, np.ndarray]:
    """
    Analyze invasion for multiple runs of a specific ECM type.
    
    Returns dict with:
        - times: common time array
        - percentile_90_mean: average 90th percentile across runs
        - percentile_90_std: standard deviation across runs
        - num_cells_mean: average cell count
        - num_cells_std: std of cell count
    """
    ecm_dir = base_dir / ecm_type
    
    if not ecm_dir.exists():
        print(f"Warning: {ecm_dir} not found")
        return None
    
    # Find all run directories (run_12345, run_13345, etc.)
    run_dirs = sorted([d for d in ecm_dir.iterdir() if d.is_dir() and d.name.startswith('run_')])
    
    if not run_dirs:
        print(f"Warning: No run directories found in {ecm_dir}")
        return None
    
    print(f"Analyzing {ecm_type} ECM...")
    print(f"  Found {len(run_dirs)} runs")
    
    # Analyze each run
    all_runs = []
    for run_dir in run_dirs:
        result = analyze_single_run(run_dir)
        if result is not None:
            all_runs.append(result)
    
    if not all_runs:
        print(f"  No valid data found")
        return None
    
    print(f"  Successfully loaded {len(all_runs)} runs")
    
    # Find common time grid (use the first run's times as reference)
    common_times = all_runs[0]['times']
    print(f"  Time range: {common_times[0]:.1f} - {common_times[-1]:.1f} minutes")
    
    # Interpolate all runs onto common time grid
    percentile_90_all = []
    num_cells_all = []
    
    for run_data in all_runs:
        # Interpolate percentile_90
        f_p90 = interpolate.interp1d(run_data['times'], run_data['percentile_90'], 
                                      kind='linear', bounds_error=False, fill_value='extrapolate')
        percentile_90_all.append(f_p90(common_times))
        
        # Interpolate num_cells
        f_nc = interpolate.interp1d(run_data['times'], run_data['num_cells'], 
                                     kind='linear', bounds_error=False, fill_value='extrapolate')
        num_cells_all.append(f_nc(common_times))
    
    # Convert to arrays and compute statistics
    percentile_90_all = np.array(percentile_90_all)
    num_cells_all = np.array(num_cells_all)
    
    percentile_90_mean = np.mean(percentile_90_all, axis=0)
    percentile_90_std = np.std(percentile_90_all, axis=0)
    
    num_cells_mean = np.mean(num_cells_all, axis=0)
    num_cells_std = np.std(num_cells_all, axis=0)
    
    # Convert times from minutes to hours
    times_hours = common_times / 60.0
    
    print(f"  Final invasion: {percentile_90_mean[-1]:.2f} ± {percentile_90_std[-1]:.2f} μm")
    
    return {
        'times': times_hours,
        'percentile_90_mean': percentile_90_mean,
        'percentile_90_std': percentile_90_std,
        'num_cells_mean': num_cells_mean,
        'num_cells_std': num_cells_std,
        'percentile_90_all': percentile_90_all  # Keep individual runs for speed calc
    }

def calculate_invasion_speed(times: np.ndarray, positions: np.ndarray, 
                            window_hours: float = 3.0) -> np.ndarray:
    """
    Calculate invasion speed using a sliding window.
    
    Args:
        times: time array in hours
        positions: position array (e.g., percentile_90)
        window_hours: window size in hours for speed calculation
    
    Returns:
        speeds: array of speeds in μm/hour
    """
    speeds = np.zeros_like(times)
    
    for i in range(len(times)):
        # Find points within window before current time
        mask = (times <= times[i]) & (times >= times[i] - window_hours)
        if np.sum(mask) >= 2:
            t_window = times[mask]
            p_window = positions[mask]
            
            # Linear fit
            if len(t_window) > 1:
                slope, _ = np.polyfit(t_window, p_window, 1)
                speeds[i] = slope  # Already in μm/hour
    
    return speeds

def plot_invasion_comparison(results_dict: Dict[str, Dict], output_dir: Path):
    """
    Create invasion comparison plots with averages and standard deviations.
    """
    ecm_types = ['perpendicular', 'parallel', 'random']
    colors = {
        'perpendicular': '#2ecc71',  # green
        'parallel': '#e74c3c',       # red
        'random': '#95a5a6',         # gray
    }
    labels = {
        'perpendicular': 'Perpendicular',
        'parallel': 'Parallel',
        'random': 'Random',
    }
    
    # Filter to only available results
    available_types = [t for t in ecm_types if t in results_dict and results_dict[t] is not None]
    
    if not available_types:
        print("No results available to plot")
        return
    
    # Create figure with 2 subplots
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    
    # 1. Mean invasion depth (90th percentile) with std
    ax1 = axes[0]
    for ecm_type in available_types:
        data = results_dict[ecm_type]
        times = data['times']
        mean = data['percentile_90_mean']
        std = data['percentile_90_std']
        
        # Plot mean line
        ax1.plot(times, mean, color=colors[ecm_type], 
                label=labels[ecm_type], linewidth=2.5)
        
        # Plot shaded std region
        ax1.fill_between(times, mean - std, mean + std, 
                        color=colors[ecm_type], alpha=0.2)
    
    ax1.set_xlabel('Time (hours)', fontsize=FONT_SIZE_LABEL, fontweight='bold')
    ax1.set_ylabel('90th Percentile Depth (μm)', fontsize=FONT_SIZE_LABEL, fontweight='bold')
    # ax1.set_title('Cell Invasion Front\n(90th Percentile)', 
    #              fontsize=FONT_SIZE_TITLE, fontweight='bold')
    ax1.set_title('(A)', 
                 fontsize=FONT_SIZE_TITLE, fontweight='bold')
    ax1.tick_params(labelsize=FONT_SIZE_TICK)
    ax1.grid(True, alpha=0.3)
    
    # 2. Invasion speed with std
    ax2 = axes[1]
    for ecm_type in available_types:
        data = results_dict[ecm_type]
        times = data['times']
        
        # Calculate speeds for each run
        all_speeds = []
        for run_p90 in data['percentile_90_all']:
            speed = calculate_invasion_speed(times, run_p90, window_hours=3.0)
            # Smooth with Gaussian filter
            speed_smooth = gaussian_filter1d(speed, sigma=2)
            all_speeds.append(speed_smooth)
        
        all_speeds = np.array(all_speeds)
        speed_mean = np.mean(all_speeds, axis=0)
        speed_std = np.std(all_speeds, axis=0)
        
        # Plot mean speed
        ax2.plot(times, speed_mean, color=colors[ecm_type], 
                label=labels[ecm_type], linewidth=2.5, alpha=0.8)
        
        # Plot shaded std region
        ax2.fill_between(times, speed_mean - speed_std, speed_mean + speed_std,
                        color=colors[ecm_type], alpha=0.2)
    
    ax2.set_xlabel('Time (hours)', fontsize=FONT_SIZE_LABEL, fontweight='bold')
    ax2.set_ylabel('Invasion Speed (μm/hour)', fontsize=FONT_SIZE_LABEL, fontweight='bold')
    ax2.set_title('(B)', 
                 fontsize=FONT_SIZE_TITLE, fontweight='bold')
    # ax2.set_title('Invasion Speed Over Time\n(3h window)', 
    #          fontsize=FONT_SIZE_TITLE, fontweight='bold')

                
    ax2.tick_params(labelsize=FONT_SIZE_TICK)
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    
    # Add single legend at the bottom center
    handles, labels_list = ax1.get_legend_handles_labels()
    fig.legend(handles, labels_list, loc='lower center', bbox_to_anchor=(0.5, 0.02), 
              ncol=3, framealpha=0.9, fontsize=FONT_SIZE_LEGEND)
    
    plt.tight_layout(rect=[0, 0.15, 1, 1])
    
    # Save figure
    output_file = output_dir / 'invasion_analysis_averaged.pdf'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def main():
    """Main execution function."""
    print("Analyzing ECM invasion simulations (averaged over multiple runs)...")
    
    # Find script directory and project paths
    script_dir = Path(__file__).parent.resolve()
    print(f"Script directory: {script_dir}")
    
    # Find Thesis directory (go up until we find it)
    thesis_dir = script_dir
    while thesis_dir.name != 'Thesis' and thesis_dir != thesis_dir.parent:
        thesis_dir = thesis_dir.parent
    
    if thesis_dir.name != 'Thesis':
        thesis_dir = Path.home() / 'Thesis'
    
    print(f"Thesis directory: {thesis_dir}")
    
    # Base directory for simulation outputs
    base_dir = thesis_dir / 'testoutput' / 'DynamicECMInvasion2d'
    print(f"Looking for data in: {base_dir}")
    print(f"Base dir exists: {base_dir.exists()}")
    print()
    
    if not base_dir.exists():
        print(f"Error: Base directory not found: {base_dir}")
        return
    
    # Analyze each ECM type
    results = {}
    for ecm_type in ['perpendicular', 'parallel', 'random']:
        results[ecm_type] = analyze_ecm_type_multiple_runs(ecm_type, base_dir)
        print()
    
    # Create output directory
    output_dir = script_dir / 'invasion_analysis_output'
    output_dir.mkdir(exist_ok=True)
    
    # Generate plots
    print("Generating plots...")
    plot_invasion_comparison(results, output_dir)
    
    # Save summary statistics
    summary_file = output_dir / 'invasion_summary_averaged.txt'
    with open(summary_file, 'w') as f:
        f.write("ECM Invasion Analysis Summary (Averaged over 20 runs)\n")
        f.write("=" * 60 + "\n\n")
        
        for ecm_type in ['perpendicular', 'parallel', 'random']:
            if ecm_type in results and results[ecm_type] is not None:
                data = results[ecm_type]
                f.write(f"{ecm_type.upper()} ECM:\n")
                f.write(f"  Final invasion (90th percentile): {data['percentile_90_mean'][-1]:.2f} ± {data['percentile_90_std'][-1]:.2f} μm\n")
                f.write(f"  Final cell count: {data['num_cells_mean'][-1]:.1f} ± {data['num_cells_std'][-1]:.1f}\n")
                f.write(f"  Time range: {data['times'][0]:.1f} - {data['times'][-1]:.1f} hours\n")
                f.write("\n")
    
    print(f"Saved: {summary_file}")
    print(f"\n✓ Analysis complete!")
    print(f"Output saved to: {output_dir}")

if __name__ == '__main__':
    main()
