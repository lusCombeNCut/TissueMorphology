#!/usr/bin/env python3
"""
analyse_crypt_budding.py

Post-processing for CryptBuddingApp stiffness sweep results.
Supports all four model types: node2d, vertex2d, node3d, vertex3d.

Crypt detection:
  - 2D: polar coordinates (theta, r) from centroid, peak detection in r(theta)
  - 3D: spherical coordinates, peak detection in radial profile by (theta, phi) binning

Output directory structure expected (as produced by CryptBuddingApp):
  <base>/CryptBudding/<model>/stiffness_<X>/run_<N>/

Usage:
  # Analyse a single model type:
  python analyse_crypt_budding.py --base-dir /path/to/testoutput --model node2d

  # Analyse all models found:
  python analyse_crypt_budding.py --base-dir /path/to/testoutput

  # Compare specific models:
  python analyse_crypt_budding.py --base-dir /path/to/testoutput --model node2d vertex2d

  # Custom output folder:
  python analyse_crypt_budding.py --base-dir /path/to/testoutput -o results/
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
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("WARNING: matplotlib not available. Text-only output.")

try:
    from scipy.signal import find_peaks, savgol_filter
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("WARNING: scipy not available. Using simple peak detection.")

# SimpleCryptCount method (optional)
try:
    from simple_crypt_count import (count_crypts_simple_method, DEFAULT_PARAMS,
                                    load_final_outline, plot_debug_analysis,
                                    plot_crypt_outline, load_final_vertex_boundary)
    HAS_SIMPLE_CRYPT_COUNT = True
except ImportError:
    HAS_SIMPLE_CRYPT_COUNT = False


# =====================================================================
# Data loading
# =====================================================================

def load_positions_from_vtu(vtu_path, dim=2):
    """Load cell positions from a Chaste VTU file."""
    try:
        import vtk
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(vtu_path)
        reader.Update()
        output = reader.GetOutput()
        points = output.GetPoints()
        n = points.GetNumberOfPoints()
        coords = np.array([points.GetPoint(i)[:dim] for i in range(n)])
        return coords
    except ImportError:
        pass

    import xml.etree.ElementTree as ET
    tree = ET.parse(vtu_path)
    for piece in tree.getroot().iter('Piece'):
        for points in piece.iter('Points'):
            for da in points.iter('DataArray'):
                if da.get('format', 'ascii') != 'ascii':
                    raise ValueError(f"Binary VTU — install vtk: pip install vtk")
                vals = [float(v) for v in da.text.strip().split()]
                n = len(vals) // 3
                coords = np.array([[vals[3*i + d] for d in range(dim)] for i in range(n)])
                return coords
    raise ValueError(f"Could not parse {vtu_path}")


def load_positions_from_viznodes(path, dim=2):
    """Load final cell positions from results.viznodes."""
    last_line = None
    with open(path) as f:
        for line in f:
            if line.strip():
                last_line = line.strip()
    if not last_line:
        raise ValueError(f"Empty viznodes: {path}")

    parts = last_line.split()
    coords = [float(v) for v in parts[1:]]
    n = len(coords) // dim
    return np.array(coords[:n * dim]).reshape(n, dim)


def load_final_positions(data_dir, dim=2):
    """Load positions, trying viznodes then VTU."""
    viznodes = os.path.join(data_dir, 'results.viznodes')
    if os.path.exists(viznodes):
        return load_positions_from_viznodes(viznodes, dim)

    vtu_files = sorted(glob.glob(os.path.join(data_dir, 'results_*.vtu')))
    if vtu_files:
        return load_positions_from_vtu(vtu_files[-1], dim)

    raise FileNotFoundError(f"No position data in {data_dir}")


def load_summary_csv(csv_path):
    """Load crypt_summary.csv into a dict of lists."""
    data = defaultdict(list)
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            for key, val in row.items():
                try:
                    data[key].append(float(val))
                except ValueError:
                    data[key].append(val)
    return dict(data)


# =====================================================================
# Crypt detection — 2D (polar)
# =====================================================================

def count_crypts_2d(positions, min_prominence=0.5, min_angular_sep_deg=30.0):
    """
    Detect outward radial protrusions in a 2D organoid.
    Returns (num_crypts, crypt_angles, (theta_bins, r_smooth)).
    """
    if len(positions) < 10:
        return 0, [], (np.array([]), np.array([]))

    cx, cy = positions.mean(axis=0)
    dx = positions[:, 0] - cx
    dy = positions[:, 1] - cy
    theta = np.arctan2(dy, dx)
    r = np.hypot(dx, dy)

    if r.ptp() < 1e-6:
        return 0, [], (np.array([]), np.array([]))

    n_bins = max(36, len(positions) // 2)
    edges = np.linspace(-np.pi, np.pi, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    r_profile = np.full(n_bins, np.nan)
    for i in range(n_bins):
        mask = (theta >= edges[i]) & (theta < edges[i + 1])
        if mask.any():
            r_profile[i] = r[mask].mean()

    valid = ~np.isnan(r_profile)
    if valid.sum() < 3:
        return 0, [], (centers, r_profile)
    r_profile[~valid] = np.interp(centers[~valid], centers[valid], r_profile[valid])

    # Circular smoothing
    pad = n_bins // 4
    r_pad = np.concatenate([r_profile[-pad:], r_profile, r_profile[:pad]])
    if HAS_SCIPY:
        win = min(max(5, n_bins // 5), len(r_pad))
        if win % 2 == 0:
            win -= 1
        r_smooth = savgol_filter(r_pad, win, 2)[pad:pad + n_bins] if win >= 3 else r_pad[pad:pad + n_bins]
    else:
        k = max(3, n_bins // 10)
        r_smooth = np.convolve(r_pad, np.ones(k) / k, mode='same')[pad:pad + n_bins]

    # Peak detection
    min_dist = max(2, int(min_angular_sep_deg / (360.0 / n_bins)))

    if HAS_SCIPY:
        r_circ = np.tile(r_smooth, 3)
        peaks, _ = find_peaks(r_circ, distance=min_dist, prominence=min_prominence)
        peaks = peaks[(peaks >= n_bins) & (peaks < 2 * n_bins)] - n_bins
        peaks = np.unique(peaks)
    else:
        peaks = []
        for i in range(n_bins):
            if (r_smooth[i] > r_smooth[(i - 1) % n_bins] and
                    r_smooth[i] > r_smooth[(i + 1) % n_bins]):
                baseline = min(r_smooth[(i - min_dist) % n_bins],
                               r_smooth[(i + min_dist) % n_bins])
                if r_smooth[i] - baseline > min_prominence:
                    peaks.append(i)
        peaks = np.array(peaks)
        if len(peaks) > 1:
            filtered = [peaks[0]]
            for p in peaks[1:]:
                if min(abs(p - filtered[-1]), n_bins - abs(p - filtered[-1])) >= min_dist:
                    filtered.append(p)
            peaks = np.array(filtered)

    angles = centers[peaks] if len(peaks) > 0 else []
    return len(peaks), angles, (centers, r_smooth)


# =====================================================================
# Crypt detection — 3D (spherical)
# =====================================================================

def count_crypts_3d(positions, min_prominence=1.0, min_sep_deg=30.0):
    """
    Detect outward radial protrusions on a 3D organoid sphere.
    Uses spherical-coordinate binning of r(theta, phi).
    Returns (num_crypts, peak_locations, profile_data).
    """
    if len(positions) < 20:
        return 0, [], {}

    centroid = positions.mean(axis=0)
    rel = positions - centroid
    r = np.linalg.norm(rel, axis=1)

    if r.ptp() < 1e-6:
        return 0, [], {}

    # Spherical coords
    theta = np.arctan2(rel[:, 1], rel[:, 0])  # azimuthal [-pi, pi]
    phi = np.arccos(np.clip(rel[:, 2] / r, -1, 1))  # polar [0, pi]

    # Bin into (theta, phi) grid
    n_theta = max(12, int(np.sqrt(len(positions))))
    n_phi = max(6, n_theta // 2)
    theta_edges = np.linspace(-np.pi, np.pi, n_theta + 1)
    phi_edges = np.linspace(0, np.pi, n_phi + 1)

    r_grid = np.full((n_theta, n_phi), np.nan)
    for it in range(n_theta):
        for ip in range(n_phi):
            mask = ((theta >= theta_edges[it]) & (theta < theta_edges[it + 1]) &
                    (phi >= phi_edges[ip]) & (phi < phi_edges[ip + 1]))
            if mask.any():
                r_grid[it, ip] = r[mask].mean()

    # Fill NaNs with nearest-neighbour (simple)
    from scipy.ndimage import generic_filter
    def nanmean_filter(vals):
        v = vals[~np.isnan(vals)]
        return v.mean() if len(v) > 0 else np.nan

    valid_frac = (~np.isnan(r_grid)).mean()
    if valid_frac < 0.2:
        return 0, [], {}

    # Simple iterative NaN filling
    for _ in range(5):
        if not np.isnan(r_grid).any():
            break
        filled = r_grid.copy()
        for it in range(n_theta):
            for ip in range(n_phi):
                if np.isnan(r_grid[it, ip]):
                    neighbours = []
                    for di in [-1, 0, 1]:
                        for dj in [-1, 0, 1]:
                            if di == 0 and dj == 0:
                                continue
                            ni = (it + di) % n_theta
                            nj = ip + dj
                            if 0 <= nj < n_phi and not np.isnan(r_grid[ni, nj]):
                                neighbours.append(r_grid[ni, nj])
                    if neighbours:
                        filled[it, ip] = np.mean(neighbours)
        r_grid = filled

    # Flatten to 1D for peak detection (treat as unwrapped sphere surface)
    # Simple approach: count peaks in each phi band
    total_peaks = 0
    peak_locs = []

    for ip in range(n_phi):
        band = r_grid[:, ip]
        if np.isnan(band).any():
            continue
        mean_r = np.nanmean(band)
        # Skip polar caps (too few cells for reliable detection)
        phi_center = 0.5 * (phi_edges[ip] + phi_edges[ip + 1])
        if phi_center < 0.3 or phi_center > np.pi - 0.3:
            continue

        # Peak detection in this phi band (circular in theta)
        min_dist = max(2, int(min_sep_deg / (360.0 / n_theta)))
        if HAS_SCIPY:
            band_circ = np.tile(band, 3)
            peaks, _ = find_peaks(band_circ, distance=min_dist,
                                  prominence=min_prominence)
            peaks = peaks[(peaks >= n_theta) & (peaks < 2 * n_theta)] - n_theta
            peaks = np.unique(peaks)
        else:
            peaks = []
            for i in range(n_theta):
                if (band[i] > band[(i-1) % n_theta] and
                        band[i] > band[(i+1) % n_theta] and
                        band[i] - mean_r > min_prominence):
                    peaks.append(i)
            peaks = np.array(peaks)

        total_peaks += len(peaks)
        for p in peaks:
            theta_c = 0.5 * (theta_edges[p] + theta_edges[p + 1])
            peak_locs.append((theta_c, phi_center))

    return total_peaks, peak_locs, {'r_grid': r_grid}


# =====================================================================
# Shape metrics
# =====================================================================

def compute_circularity_2d(positions):
    """
    Circularity of a 2D point cloud (convex hull).
    circularity = 4*pi*Area / Perimeter^2
    Returns 1.0 for a perfect circle, <1 for irregular shapes.
    """
    if len(positions) < 3:
        return np.nan
    try:
        from scipy.spatial import ConvexHull
        hull = ConvexHull(positions)
        area = hull.volume  # in 2D, 'volume' gives area
        verts = positions[hull.vertices]
        perimeter = 0.0
        for i in range(len(verts)):
            perimeter += np.linalg.norm(verts[i] - verts[(i + 1) % len(verts)])
        if perimeter < 1e-12:
            return np.nan
        return 4.0 * np.pi * area / (perimeter ** 2)
    except Exception:
        return np.nan


def compute_sphericity_3d(positions):
    """
    Sphericity of a 3D point cloud (convex hull).
    sphericity = (pi^(1/3) * (6V)^(2/3)) / A
    Returns 1.0 for a perfect sphere, <1 for irregular shapes.
    """
    if len(positions) < 4:
        return np.nan
    try:
        from scipy.spatial import ConvexHull
        hull = ConvexHull(positions)
        volume = hull.volume
        area = hull.area
        if area < 1e-12:
            return np.nan
        return (np.pi ** (1.0 / 3.0) * (6.0 * volume) ** (2.0 / 3.0)) / area
    except Exception:
        return np.nan


# =====================================================================
# Analysis
# =====================================================================

def analyse_model(base_dir, model_type, method='polar', simple_params=None, use_outline=False,
                  debug_plots_dir=None):
    """
    Scan output directory for a model and count crypts per (stiffness, run).
    
    Args:
        base_dir: Base output directory
        model_type: Model type (node2d, vertex2d, etc.)
        method: 'polar' (radial peak detection) or 'simple' (SimpleCryptCount)
        simple_params: Parameters dict for SimpleCryptCount method
        use_outline: If True, use VTP outline files for SimpleCryptCount method (2D only)
        debug_plots_dir: If provided, save debug plots for run_1 of each stiffness
    
    Returns: dict { stiffness: [{'num_crypts': int, 'circularity': float,
                                  'sphericity': float}, ...] }
    """
    # Try both with and without 'CryptBudding' subdirectory
    model_dir = os.path.join(base_dir, 'CryptBudding', model_type)
    if not os.path.isdir(model_dir):
        model_dir = os.path.join(base_dir, model_type)
    if not os.path.isdir(model_dir):
        print(f"  Directory not found: {model_dir}")
        return {}

    dim = 3 if '3d' in model_type else 2
    results = defaultdict(list)
    
    # Validate method
    if method == 'simple' and dim == 3:
        print(f"  WARNING: SimpleCryptCount method only supports 2D. Falling back to polar.")
        method = 'polar'
    if method == 'simple' and not HAS_SIMPLE_CRYPT_COUNT:
        print(f"  WARNING: SimpleCryptCount not available. Falling back to polar.")
        method = 'polar'

    stiffness_dirs = sorted(glob.glob(os.path.join(model_dir, 'stiffness_*')))
    if not stiffness_dirs:
        print(f"  No stiffness_* directories in {model_dir}")
        return results

    for s_dir in stiffness_dirs:
        try:
            stiffness = float(os.path.basename(s_dir).replace('stiffness_', ''))
        except ValueError:
            continue

        run_dirs = sorted(glob.glob(os.path.join(s_dir, 'run_*')))
        for r_dir in run_dirs:
            run_label = os.path.basename(r_dir)

            # Search for output in the run dir and any results_from_time_* sub
            search_dirs = [r_dir]
            for sub in sorted(glob.glob(os.path.join(r_dir, 'results_from_time_*'))):
                if os.path.isdir(sub):
                    search_dirs.append(sub)

            positions = None
            outline_data = None
            cell_types = None
            
            for sd in search_dirs:
                # Try to load outline data (for SimpleCryptCount with --use-outline)
                if method == 'simple' and use_outline:
                    try:
                        outline_data, cell_types = load_final_outline(sd)
                        break
                    except (FileNotFoundError, ValueError):
                        pass
                
                # Try to load vertex mesh boundary (for vertex2d model)
                if method == 'simple' and 'vertex' in model_type and dim == 2:
                    try:
                        outline_data, cell_types = load_final_vertex_boundary(sd)
                        break
                    except (FileNotFoundError, ValueError, Exception):
                        pass
                
                # Fall back to cell positions
                try:
                    positions = load_final_positions(sd, dim)
                    break
                except (FileNotFoundError, ValueError):
                    continue

            if positions is None and outline_data is None:
                print(f"  SKIP {model_type}/stiffness_{stiffness}/{run_label}: no data")
                continue

            if dim == 2:
                if method == 'simple':
                    # Use SimpleCryptCount method
                    if outline_data is not None:
                        # Use VTP outline or vertex mesh directly
                        result = count_crypts_simple_method(outline_data, simple_params,
                                                           boundary_is_ordered=True)
                    else:
                        # Use cell positions
                        result = count_crypts_simple_method(positions, simple_params)
                    n_crypts = result.num_crypts
                    circ = result.circularity
                    
                    # Generate debug plot for run_1 of each stiffness
                    if debug_plots_dir and run_label == 'run_1' and HAS_MATPLOTLIB:
                        os.makedirs(debug_plots_dir, exist_ok=True)
                        debug_path = os.path.join(debug_plots_dir, 
                                                  f'{model_type}_stiffness_{stiffness}_crypt_outline.png')
                        title = f'{model_type} - ECM Stiffness = {stiffness}'
                        plot_crypt_outline(result, output_path=debug_path, title=title)
                else:
                    # Use polar radial detection method
                    n_crypts, crypt_angles, profile_data = count_crypts_2d(positions)
                    circ = compute_circularity_2d(positions)
                sph = np.nan
            else:
                n_crypts, _, _ = count_crypts_3d(positions)
                circ = np.nan
                sph = compute_sphericity_3d(positions)

            shape_val = circ if dim == 2 else sph
            shape_name = 'circularity' if dim == 2 else 'sphericity'
            results[stiffness].append({
                'num_crypts': n_crypts,
                'circularity': circ,
                'sphericity': sph,
            })
            print(f"  {model_type}/stiffness_{stiffness}/{run_label}: "
                  f"{n_crypts} crypts, {shape_name}={shape_val:.4f}")

    return dict(results)


# =====================================================================
# Plotting
# =====================================================================

COLORS = {
    'node2d': 'steelblue',
    'vertex2d': 'coral',
    'node3d': 'forestgreen',
    'vertex3d': 'mediumpurple',
}

LABELS = {
    'node2d': '2D Node-based',
    'vertex2d': '2D Vertex-based',
    'node3d': '3D Node-based',
    'vertex3d': '3D Vertex-based',
}


def plot_single_model(results, model_type, output_path):
    """Box + scatter plot for one model."""
    if not HAS_MATPLOTLIB or not results:
        return

    color = COLORS.get(model_type, 'gray')
    label = LABELS.get(model_type, model_type)
    stiffnesses = sorted(results.keys())
    data = [[d['num_crypts'] for d in results[s]] for s in stiffnesses]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Box plot
    ax1.boxplot(data, positions=range(len(stiffnesses)), widths=0.5,
                patch_artist=True,
                boxprops=dict(facecolor=color, alpha=0.4),
                medianprops=dict(color='black', linewidth=2))
    for i, counts in enumerate(data):
        jitter = np.random.normal(0, 0.08, len(counts))
        ax1.scatter(np.full(len(counts), i) + jitter, counts,
                    color=color, alpha=0.6, s=30, zorder=5)
    ax1.set_xticks(range(len(stiffnesses)))
    ax1.set_xticklabels([f'{s:.1f}' for s in stiffnesses], rotation=45)
    ax1.set_xlabel('ECM Stiffness')
    ax1.set_ylabel('Number of Crypts')
    ax1.set_title(f'{label}: Crypts vs ECM Stiffness')
    ax1.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.grid(axis='y', alpha=0.3)

    # Mean ± SD (log-scale x)
    means = [np.mean(d) if d else 0 for d in data]
    stds = [np.std(d) if d else 0 for d in data]
    ax2.errorbar(stiffnesses, means, yerr=stds, fmt='o-',
                 color=color, capsize=5, markersize=8, linewidth=2)
    ax2.set_xlabel('ECM Stiffness')
    ax2.set_ylabel('Number of Crypts (mean ± SD)')
    ax2.set_title(f'{label}: Mean Crypts vs Stiffness')
    ax2.set_xscale('log')
    ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {output_path}")
    plt.close()


def plot_comparison(all_results, output_path):
    """Overlay multiple models on one plot."""
    if not HAS_MATPLOTLIB:
        return

    fig, ax = plt.subplots(figsize=(9, 5))

    markers = ['o', 's', '^', 'D']
    for i, (model_type, results) in enumerate(sorted(all_results.items())):
        if not results:
            continue
        color = COLORS.get(model_type, 'gray')
        label = LABELS.get(model_type, model_type)
        stiffnesses = sorted(results.keys())
        counts = {s: [d['num_crypts'] for d in results[s]] for s in stiffnesses}
        means = [np.mean(counts[s]) for s in stiffnesses]
        stds = [np.std(counts[s]) for s in stiffnesses]
        ax.errorbar(stiffnesses, means, yerr=stds,
                    fmt=f'{markers[i % len(markers)]}-',
                    color=color, capsize=5, markersize=8,
                    linewidth=2, label=label)

    ax.set_xlabel('ECM Stiffness', fontsize=12)
    ax.set_ylabel('Number of Crypts (mean ± SD)', fontsize=12)
    ax.set_title('Crypt Budding vs ECM Stiffness', fontsize=13)
    ax.set_xscale('log')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend(fontsize=11)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {output_path}")
    plt.close()


# =====================================================================
# Summary CSV / table
# =====================================================================

def print_summary_table(results, model_type):
    label = LABELS.get(model_type, model_type)
    dim = 3 if '3d' in model_type else 2
    shape_name = 'Sphericity' if dim == 3 else 'Circularity'
    shape_key = 'sphericity' if dim == 3 else 'circularity'

    print(f"\n{'=' * 72}")
    print(f"  {label}: Crypt Count & {shape_name} Summary")
    print(f"{'=' * 72}")
    print(f"  {'Stiffness':>10}  {'N':>4}  {'Mean':>6}  {'SD':>6}  "
          f"{'Min':>4}  {'Max':>4}  {shape_name:>11}")
    print(f"  {'-'*10}  {'-'*4}  {'-'*6}  {'-'*6}  "
          f"{'-'*4}  {'-'*4}  {'-'*11}")
    for s in sorted(results.keys()):
        entries = results[s]
        if entries:
            c = [d['num_crypts'] for d in entries]
            sv = [d[shape_key] for d in entries]
            sv_valid = [v for v in sv if not np.isnan(v)]
            sv_str = f"{np.mean(sv_valid):.4f}" if sv_valid else "N/A"
            print(f"  {s:10.1f}  {len(c):4d}  {np.mean(c):6.2f}  "
                  f"{np.std(c):6.2f}  {min(c):4d}  {max(c):4d}  "
                  f"{sv_str:>11}")
    print(f"{'=' * 72}\n")


def save_results_csv(all_results, output_path):
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['model', 'stiffness', 'replicate', 'num_crypts',
                         'circularity', 'sphericity'])
        for model_type in sorted(all_results):
            for s in sorted(all_results[model_type]):
                for i, entry in enumerate(all_results[model_type][s]):
                    circ = entry['circularity']
                    sph = entry['sphericity']
                    writer.writerow([
                        model_type, s, i, entry['num_crypts'],
                        f"{circ:.6f}" if not np.isnan(circ) else '',
                        f"{sph:.6f}" if not np.isnan(sph) else '',
                    ])
    print(f"  Saved CSV: {output_path}")


def load_summary_results(base_dir, model_type):
    """
    Load time-series summary data from crypt_summary.csv files.
    Returns dict { stiffness: { run: DataFrame-like dict } }
    """
    # Try both with and without 'CryptBudding' subdirectory
    model_dir = os.path.join(base_dir, 'CryptBudding', model_type)
    if not os.path.isdir(model_dir):
        model_dir = os.path.join(base_dir, model_type)
    results = defaultdict(dict)

    for s_dir in sorted(glob.glob(os.path.join(model_dir, 'stiffness_*'))):
        try:
            stiffness = float(os.path.basename(s_dir).replace('stiffness_', ''))
        except ValueError:
            continue
        for r_dir in sorted(glob.glob(os.path.join(s_dir, 'run_*'))):
            run_num = os.path.basename(r_dir)
            csv_path = os.path.join(r_dir, 'crypt_summary.csv')
            if os.path.exists(csv_path):
                results[stiffness][run_num] = load_summary_csv(csv_path)

    return dict(results)


# =====================================================================
# Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Analyse CryptBuddingApp output: count crypts vs ECM stiffness')

    parser.add_argument('--base-dir', type=str, required=True,
                        help='Base output directory (contains CryptBudding/<model>/...)')
    parser.add_argument('--model', nargs='+', default=None,
                        choices=['node2d', 'vertex2d', 'node3d', 'vertex3d'],
                        help='Model type(s) to analyse (default: all found)')
    parser.add_argument('-o', '--output-dir', type=str,
                        default='crypt_budding_analysis',
                        help='Output directory for plots/CSVs')
    parser.add_argument('--min-prominence', type=float, default=0.5,
                        help='Minimum radial prominence for crypt detection (polar method)')
    
    # SimpleCryptCount method options
    parser.add_argument('--method', type=str, default='polar',
                        choices=['polar', 'simple'],
                        help='Crypt counting method: polar (radial peaks) or simple (SimpleCryptCount)')
    parser.add_argument('--use-outline', action='store_true',
                        help='Use VTP outline files for SimpleCryptCount (more accurate for 2D)')
    parser.add_argument('--fourier-harmonics', type=int, default=25,
                        help='Fourier harmonics for SimpleCryptCount method')
    parser.add_argument('--min-area', type=float, default=0.0666,
                        help='Min normalized crypt area (SimpleCryptCount)')
    parser.add_argument('--max-area', type=float, default=0.2736,
                        help='Max normalized crypt area (SimpleCryptCount)')
    parser.add_argument('--min-arc-length', type=float, default=0.1466,
                        help='Min normalized arc length (SimpleCryptCount)')
    parser.add_argument('--debug-plots', action='store_true',
                        help='Save debug plots showing crypt outlines for run_1 of each stiffness')

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Set up SimpleCryptCount parameters if using that method
    simple_params = None
    use_outline = args.use_outline if hasattr(args, 'use_outline') else False
    
    if args.method == 'simple':
        if not HAS_SIMPLE_CRYPT_COUNT:
            print("ERROR: SimpleCryptCount method not available. ")
            print("       Check that simple_crypt_count.py is in the same directory.")
            sys.exit(1)
        simple_params = DEFAULT_PARAMS.copy()
        simple_params['fourier_harmonics'] = args.fourier_harmonics
        simple_params['min_area'] = args.min_area
        simple_params['max_area'] = args.max_area
        simple_params['min_arc_length'] = args.min_arc_length
        print(f"Using SimpleCryptCount method with parameters:")
        print(f"  Fourier harmonics: {simple_params['fourier_harmonics']}")
        print(f"  Area range: [{simple_params['min_area']:.4f}, {simple_params['max_area']:.4f})")
        print(f"  Min arc length: {simple_params['min_arc_length']:.4f}")
        if use_outline:
            print(f"  Using VTP outline files for boundary data")

    # Determine which models to analyse
    all_models = ['node2d', 'vertex2d', 'node3d', 'vertex3d']
    if args.model:
        models = args.model
    else:
        models = []
        for m in all_models:
            d = os.path.join(args.base_dir, 'CryptBudding', m)
            if not os.path.isdir(d):
                d = os.path.join(args.base_dir, m)
            if os.path.isdir(d):
                models.append(m)
        if not models:
            print(f"No model directories found under {args.base_dir}/")
            sys.exit(1)

    print(f"Models to analyse: {', '.join(models)}")
    print(f"Base directory:    {args.base_dir}")
    print(f"Output directory:  {args.output_dir}")

    all_results = {}

    for model_type in models:
        print(f"\n--- Analysing {LABELS.get(model_type, model_type)} ---")
        
        # Set up debug plots directory if requested (only for SimpleCryptCount method)
        debug_plots_dir = None
        if args.debug_plots and args.method == 'simple':
            debug_plots_dir = os.path.join(args.output_dir, 'debug_plots')
        
        results = analyse_model(args.base_dir, model_type, 
                               method=args.method, simple_params=simple_params,
                               use_outline=use_outline,
                               debug_plots_dir=debug_plots_dir)

        if results:
            all_results[model_type] = results
            print_summary_table(results, model_type)

            if HAS_MATPLOTLIB:
                method_suffix = f'_{args.method}' if args.method != 'polar' else ''
                if use_outline:
                    method_suffix += '_outline'
                plot_single_model(
                    results, model_type,
                    os.path.join(args.output_dir, f'{model_type}_crypts_vs_stiffness{method_suffix}.png'))

    if all_results:
        save_results_csv(all_results,
                         os.path.join(args.output_dir, 'all_crypt_counts.csv'))

        if len(all_results) > 1 and HAS_MATPLOTLIB:
            plot_comparison(all_results,
                           os.path.join(args.output_dir, 'comparison_crypts_vs_stiffness.png'))

    print("\nDone.")


if __name__ == '__main__':
    main()
