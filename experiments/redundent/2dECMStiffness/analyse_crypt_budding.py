#!/usr/bin/env python3
"""
analyse_crypt_budding.py

Post-processing script for crypt budding stiffness sweep simulations.
Reads Chaste VTU output and summary CSVs to:
  1. Count the number of crypt buds formed at each stiffness level
  2. Plot distributions: number of crypts vs ECM stiffness
  3. Generate summary statistics for comparison with experimental data

Crypt detection method (polar coordinates for circular organoids):
  - Load final cell positions from VTU or viznodes
  - Compute centroid and convert (x,y) -> (theta, r)
  - Smooth the r(theta) profile around the ring
  - Detect local maxima in r(theta) = outward radial protrusions = crypt buds

Usage:
  python analyse_crypt_budding.py --data-dir /path/to/testoutput/CryptBudding2d_NodeBased
  python analyse_crypt_budding.py --data-dir /path/to/testoutput/CryptBudding2d_VertexBased

  # Or analyse both and overlay:
  python analyse_crypt_budding.py \\
      --node-dir /path/to/CryptBudding2d_NodeBased \\
      --vertex-dir /path/to/CryptBudding2d_VertexBased
"""

import os
import sys
import glob
import argparse
import numpy as np
import csv
from collections import defaultdict

try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend for HPC
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("WARNING: matplotlib not available. Text-only output.")

try:
    from scipy.signal import find_peaks, savgol_filter
    from scipy.ndimage import gaussian_filter1d
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("WARNING: scipy not available. Using simple peak detection.")


def load_final_positions_from_vtu(vtu_path):
    """
    Parse cell positions from a Chaste VTU (VTK XML) file.
    Returns arrays of x, y coordinates.
    Handles both ASCII inline and binary/appended VTU formats.
    """
    # Try VTK Python library first (handles all formats)
    try:
        import vtk
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(vtu_path)
        reader.Update()
        output = reader.GetOutput()
        points = output.GetPoints()
        n_points = points.GetNumberOfPoints()
        x = np.array([points.GetPoint(i)[0] for i in range(n_points)])
        y = np.array([points.GetPoint(i)[1] for i in range(n_points)])
        return x, y
    except ImportError:
        pass

    # Fallback: parse ASCII inline VTU
    import xml.etree.ElementTree as ET

    tree = ET.parse(vtu_path)
    root = tree.getroot()

    for piece in root.iter('Piece'):
        for points in piece.iter('Points'):
            for data_array in points.iter('DataArray'):
                fmt = data_array.get('format', 'ascii')
                if fmt != 'ascii':
                    raise ValueError(
                        f"VTU uses '{fmt}' format — install the vtk Python "
                        f"package (`pip install vtk`) to read binary VTU files, "
                        f"or ensure results.viznodes is available.")
                text = data_array.text
                if text is None:
                    continue
                values = [float(v) for v in text.strip().split()]
                n_points = len(values) // 3
                x = np.array([values[3*i] for i in range(n_points)])
                y = np.array([values[3*i+1] for i in range(n_points)])
                return x, y

    raise ValueError(f"Could not parse positions from {vtu_path}")


def load_positions_from_viznodes(viznodes_path):
    """
    Load cell positions from a Chaste results.viznodes file.

    Chaste viznodes format: each line is one timestep, starting with the
    timestep index, followed by x y coordinate pairs for every cell:
        0  x0 y0 x1 y1 x2 y2 ...
        1  x0 y0 x1 y1 x2 y2 ...

    We read the LAST line to get the final-time positions.
    """
    last_line = None
    with open(viznodes_path, 'r') as f:
        for line in f:
            stripped = line.strip()
            if stripped:
                last_line = stripped

    if last_line is None:
        raise ValueError(f"Empty viznodes file: {viznodes_path}")

    parts = last_line.split()
    # First element is the timestep index; remaining are x,y pairs
    coords = [float(v) for v in parts[1:]]
    if len(coords) < 2 or len(coords) % 2 != 0:
        raise ValueError(f"Unexpected coordinate count in {viznodes_path}")

    x = np.array(coords[0::2])
    y = np.array(coords[1::2])
    return x, y


def load_final_positions(data_dir):
    """
    Load final cell positions from a Chaste output directory.
    Tries, in order:
      1. results.viznodes  (fast, always ASCII)
      2. Last VTU file     (may need vtk library for binary format)
    """
    # 1. Try viznodes (most reliable — always ASCII in Chaste)
    viznodes = os.path.join(data_dir, 'results.viznodes')
    if os.path.exists(viznodes):
        return load_positions_from_viznodes(viznodes)

    # 2. Fallback: try VTU files
    vtu_files = sorted(glob.glob(os.path.join(data_dir, 'results_*.vtu')))
    if vtu_files:
        return load_final_positions_from_vtu(vtu_files[-1])

    raise FileNotFoundError(f"No position data found in {data_dir}")


def count_crypts_from_positions(x, y, min_prominence=0.5, min_angular_sep_deg=30.0):
    """
    Count crypt buds by detecting outward radial protrusions from a
    circular organoid using polar coordinates.

    Method:
      1. Compute centroid of all cells
      2. Convert (x, y) -> (theta, r) relative to centroid
      3. Sort by theta and compute r(theta) profile
      4. Smooth and detect local maxima in r(theta)
      5. Each significant radial peak = one crypt bud

    Parameters:
      x, y: cell positions
      min_prominence: minimum radial prominence (in cell diameters)
                      for a peak to count as a crypt
      min_angular_sep_deg: minimum angular separation between detected
                           crypts (degrees)

    Returns:
      num_crypts: number of detected crypt buds
      crypt_angles: theta values (radians) of crypt centres
      profile: (theta_sorted, r_smooth) for plotting
    """
    if len(x) < 10:
        return 0, [], (np.array([]), np.array([]))

    # Centroid
    cx, cy = np.mean(x), np.mean(y)

    # Polar coordinates
    dx = x - cx
    dy = y - cy
    theta = np.arctan2(dy, dx)   # range [-pi, pi]
    r = np.sqrt(dx**2 + dy**2)

    if r.max() - r.min() < 1e-6:
        return 0, [], (np.array([]), np.array([]))

    # Sort by angle
    order = np.argsort(theta)
    theta_sorted = theta[order]
    r_sorted = r[order]

    # Bin into angular bins for a smooth profile
    n_bins = max(36, len(x) // 2)  # at least 10-degree resolution
    bin_edges = np.linspace(-np.pi, np.pi, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # Mean radius in each angular bin
    r_profile = np.full(n_bins, np.nan)
    for i in range(n_bins):
        mask = (theta >= bin_edges[i]) & (theta < bin_edges[i+1])
        if mask.any():
            r_profile[i] = np.mean(r[mask])

    # Interpolate NaN bins (sparse regions)
    valid = ~np.isnan(r_profile)
    if valid.sum() < 3:
        return 0, [], (bin_centers, r_profile)

    r_profile[~valid] = np.interp(
        bin_centers[~valid], bin_centers[valid], r_profile[valid])

    # Smooth the profile (circular: wrap edges)
    if HAS_SCIPY:
        # Pad circularly for smooth wraparound
        pad = n_bins // 4
        r_padded = np.concatenate([r_profile[-pad:], r_profile, r_profile[:pad]])
        window = min(max(5, n_bins // 5), len(r_padded))
        if window % 2 == 0:
            window -= 1
        if window >= 3:
            r_smooth_padded = savgol_filter(r_padded, window, 2)
        else:
            r_smooth_padded = r_padded.copy()
        r_smooth = r_smooth_padded[pad:pad + n_bins]
    else:
        # Simple moving average with circular padding
        pad = n_bins // 4
        r_padded = np.concatenate([r_profile[-pad:], r_profile, r_profile[:pad]])
        kernel_size = max(3, n_bins // 10)
        kernel = np.ones(kernel_size) / kernel_size
        r_smooth_padded = np.convolve(r_padded, kernel, mode='same')
        r_smooth = r_smooth_padded[pad:pad + n_bins]

    # Detect local maxima (outward protrusions = crypt buds)
    min_distance = max(2, int(min_angular_sep_deg / (360.0 / n_bins)))

    if HAS_SCIPY:
        # Also pad for circular peak detection
        r_circ = np.concatenate([r_smooth, r_smooth, r_smooth])
        peaks_circ, props = find_peaks(r_circ,
                                       distance=min_distance,
                                       prominence=min_prominence)
        # Keep only peaks in the middle copy
        peaks = peaks_circ[(peaks_circ >= n_bins) & (peaks_circ < 2 * n_bins)] - n_bins
        # Deduplicate (circular wrap)
        peaks = np.unique(peaks)
    else:
        # Simple peak detection
        peaks = []
        for i in range(n_bins):
            prev_i = (i - 1) % n_bins
            next_i = (i + 1) % n_bins
            if r_smooth[i] > r_smooth[prev_i] and r_smooth[i] > r_smooth[next_i]:
                # Check prominence against local baseline
                baseline = min(r_smooth[(i - min_distance) % n_bins],
                               r_smooth[(i + min_distance) % n_bins])
                if r_smooth[i] - baseline > min_prominence:
                    peaks.append(i)
        peaks = np.array(peaks)

        # Enforce minimum angular separation
        if len(peaks) > 1:
            filtered = [peaks[0]]
            for p in peaks[1:]:
                if min(abs(p - filtered[-1]), n_bins - abs(p - filtered[-1])) >= min_distance:
                    filtered.append(p)
            peaks = np.array(filtered)

    crypt_angles = bin_centers[peaks] if len(peaks) > 0 else []
    return len(peaks), crypt_angles, (bin_centers, r_smooth)


def analyse_sweep_directory(data_dir, model_label=""):
    """
    Scan a stiffness sweep output directory and count crypts for each
    (stiffness, replicate) combination.

    Expected structure:
      data_dir/
        stiffness_0.5/run_0/crypt_summary.csv
        stiffness_0.5/run_1/crypt_summary.csv
        ...
        stiffness_50.0/run_9/crypt_summary.csv

    Returns:
      results: dict mapping stiffness → list of crypt counts
    """
    results = defaultdict(list)

    # Find all stiffness directories
    stiffness_dirs = sorted(glob.glob(os.path.join(data_dir, 'stiffness_*')))

    if not stiffness_dirs:
        print(f"No stiffness_* directories found in {data_dir}")
        return results

    for s_dir in stiffness_dirs:
        # Extract stiffness value
        dirname = os.path.basename(s_dir)
        try:
            stiffness = float(dirname.replace('stiffness_', ''))
        except ValueError:
            continue

        # Find all run directories
        run_dirs = sorted(glob.glob(os.path.join(s_dir, 'run_*')))

        for r_dir in run_dirs:
            run_label = os.path.basename(r_dir)

            try:
                # Chaste may place output directly in run_* or inside
                # a results_from_time_* subdirectory.  Check both.
                search_dirs = [r_dir]
                for sub in sorted(glob.glob(os.path.join(r_dir, 'results_from_time_*'))):
                    if os.path.isdir(sub):
                        search_dirs.append(sub)

                x, y = None, None
                for search_dir in search_dirs:
                    try:
                        x, y = load_final_positions(search_dir)
                        break
                    except (FileNotFoundError, ValueError):
                        continue

                if x is None:
                    print(f"  SKIP {dirname}/{run_label}: no output data")
                    continue

                n_crypts, positions, _ = count_crypts_from_positions(x, y)
                results[stiffness].append(n_crypts)
                print(f"  {model_label} {dirname}/{run_label}: {n_crypts} crypts")

            except Exception as e:
                print(f"  ERROR {dirname}/{run_label}: {e}")
                continue

    return results


def plot_crypt_vs_stiffness(results, model_label, output_path, color='steelblue'):
    """
    Plot: number of crypts vs ECM stiffness (box plot + scatter).
    """
    if not HAS_MATPLOTLIB:
        return

    stiffnesses = sorted(results.keys())
    data = [results[s] for s in stiffnesses]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Box plot
    ax = axes[0]
    bp = ax.boxplot(data, positions=range(len(stiffnesses)),
                     widths=0.5, patch_artist=True,
                     boxprops=dict(facecolor=color, alpha=0.4),
                     medianprops=dict(color='black', linewidth=2))
    ax.set_xticks(range(len(stiffnesses)))
    ax.set_xticklabels([f'{s:.1f}' for s in stiffnesses], rotation=45)
    ax.set_xlabel('ECM Stiffness')
    ax.set_ylabel('Number of Crypts')
    ax.set_title(f'{model_label}: Crypts vs ECM Stiffness')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.grid(axis='y', alpha=0.3)

    # Scatter with jitter
    for i, (s, counts) in enumerate(zip(stiffnesses, data)):
        jitter = np.random.normal(0, 0.08, len(counts))
        ax.scatter(np.full(len(counts), i) + jitter, counts,
                   color=color, alpha=0.6, s=30, zorder=5)

    # Summary statistics
    ax2 = axes[1]
    means = [np.mean(d) if d else 0 for d in data]
    stds = [np.std(d) if d else 0 for d in data]
    ax2.errorbar(stiffnesses, means, yerr=stds, fmt='o-',
                  color=color, capsize=5, markersize=8, linewidth=2)
    ax2.set_xlabel('ECM Stiffness')
    ax2.set_ylabel('Number of Crypts (mean ± SD)')
    ax2.set_title(f'{model_label}: Mean Crypts vs Stiffness')
    ax2.set_xscale('log')
    ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {output_path}")
    plt.close()


def plot_comparison(node_results, vertex_results, output_path):
    """
    Overlay node-based and vertex-based results on the same plot.
    """
    if not HAS_MATPLOTLIB:
        return

    fig, ax = plt.subplots(figsize=(8, 5))

    for results, label, color, marker in [
        (node_results, 'Node-based', 'steelblue', 'o'),
        (vertex_results, 'Vertex-based', 'coral', 's')
    ]:
        if not results:
            continue

        stiffnesses = sorted(results.keys())
        means = [np.mean(results[s]) if results[s] else 0 for s in stiffnesses]
        stds = [np.std(results[s]) if results[s] else 0 for s in stiffnesses]

        ax.errorbar(stiffnesses, means, yerr=stds,
                     fmt=f'{marker}-', color=color, capsize=5,
                     markersize=8, linewidth=2, label=label)

    ax.set_xlabel('ECM Stiffness', fontsize=12)
    ax.set_ylabel('Number of Crypts (mean ± SD)', fontsize=12)
    ax.set_title('Crypt Budding vs ECM Stiffness\n(Node-based vs Vertex-based)', fontsize=13)
    ax.set_xscale('log')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend(fontsize=11)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {output_path}")
    plt.close()


def print_summary_table(results, model_label):
    """
    Print a text summary table of results.
    """
    print(f"\n{'='*60}")
    print(f"  {model_label}: Crypt Count Summary")
    print(f"{'='*60}")
    print(f"  {'Stiffness':>10}  {'N':>4}  {'Mean':>6}  {'SD':>6}  {'Min':>4}  {'Max':>4}")
    print(f"  {'-'*10}  {'-'*4}  {'-'*6}  {'-'*6}  {'-'*4}  {'-'*4}")

    for s in sorted(results.keys()):
        counts = results[s]
        if counts:
            print(f"  {s:10.1f}  {len(counts):4d}  {np.mean(counts):6.2f}  "
                  f"{np.std(counts):6.2f}  {min(counts):4d}  {max(counts):4d}")
        else:
            print(f"  {s:10.1f}  {'N/A':>4}  {'N/A':>6}  {'N/A':>6}  {'N/A':>4}  {'N/A':>4}")

    print(f"{'='*60}\n")


def save_results_csv(results, output_path, model_label):
    """
    Save results to a CSV file for further analysis.
    """
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['model', 'stiffness', 'replicate', 'num_crypts'])
        for s in sorted(results.keys()):
            for i, count in enumerate(results[s]):
                writer.writerow([model_label, s, i, count])
    print(f"  Saved CSV: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Analyse crypt budding simulations: count crypts vs ECM stiffness')

    parser.add_argument('--data-dir', type=str, default=None,
                        help='Path to single model output directory')
    parser.add_argument('--node-dir', type=str, default=None,
                        help='Path to node-based model output directory')
    parser.add_argument('--vertex-dir', type=str, default=None,
                        help='Path to vertex-based model output directory')
    parser.add_argument('--output-dir', type=str, default='crypt_budding_analysis',
                        help='Directory for analysis output (default: crypt_budding_analysis)')
    parser.add_argument('--min-depth', type=float, default=0.5,
                        help='Minimum crypt depth for detection (default: 0.5)')

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    node_results = {}
    vertex_results = {}

    # Single directory mode
    if args.data_dir:
        label = 'Node-based' if 'NodeBased' in args.data_dir else 'Vertex-based'
        print(f"\nAnalysing {label} results from: {args.data_dir}")
        results = analyse_sweep_directory(args.data_dir, label)
        print_summary_table(results, label)

        if results:
            save_results_csv(results, os.path.join(args.output_dir, 'crypt_counts.csv'), label)
            if HAS_MATPLOTLIB:
                plot_crypt_vs_stiffness(
                    results, label,
                    os.path.join(args.output_dir, 'crypts_vs_stiffness.png'))
        return

    # Dual directory mode
    if args.node_dir:
        print(f"\nAnalysing node-based results from: {args.node_dir}")
        node_results = analyse_sweep_directory(args.node_dir, 'Node-based')
        print_summary_table(node_results, 'Node-based')
        save_results_csv(node_results,
                          os.path.join(args.output_dir, 'node_crypt_counts.csv'),
                          'Node-based')
        if HAS_MATPLOTLIB:
            plot_crypt_vs_stiffness(
                node_results, 'Node-based',
                os.path.join(args.output_dir, 'node_crypts_vs_stiffness.png'),
                color='steelblue')

    if args.vertex_dir:
        print(f"\nAnalysing vertex-based results from: {args.vertex_dir}")
        vertex_results = analyse_sweep_directory(args.vertex_dir, 'Vertex-based')
        print_summary_table(vertex_results, 'Vertex-based')
        save_results_csv(vertex_results,
                          os.path.join(args.output_dir, 'vertex_crypt_counts.csv'),
                          'Vertex-based')
        if HAS_MATPLOTLIB:
            plot_crypt_vs_stiffness(
                vertex_results, 'Vertex-based',
                os.path.join(args.output_dir, 'vertex_crypts_vs_stiffness.png'),
                color='coral')

    # Comparison plot
    if node_results and vertex_results and HAS_MATPLOTLIB:
        plot_comparison(node_results, vertex_results,
                         os.path.join(args.output_dir, 'comparison_crypts_vs_stiffness.png'))

    if not node_results and not vertex_results:
        print("\nNo data directories specified. Use --data-dir, --node-dir, or --vertex-dir.")
        parser.print_help()


if __name__ == '__main__':
    main()
