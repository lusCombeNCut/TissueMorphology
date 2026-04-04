#!/usr/bin/env python3
"""
Generate comparison analysis plots for 2D models (Node2d and Vertex2d).

Reads:
  - all_crypt_counts.csv (from reanalyse_crypt_counts.py)
  - crypt_summary.csv files from each simulation (time-series data)

Produces:
  - Crypt count vs stiffness (boxplot + scatter)
  - Circularity vs stiffness
  - Mean radius vs stiffness
  - Radius variance vs stiffness
  - Cell count vs stiffness
  - Time-series comparison at selected stiffnesses
  - Combined multi-panel figure

Usage:
  python3 plot_2d_comparison.py
"""

import os
import sys
import csv
import glob
import json
import numpy as np
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# ===========================================================================
# Configuration
# ===========================================================================
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    'results-download', '03-31_08_56_51-merged')
ANALYSIS_DIR = os.path.join(BASE, 'analysis_ecm_stiffness_2026-03-31_08-56-51')
CRYPT_CSV = os.path.join(ANALYSIS_DIR, 'crypt_analysis_output', 'all_crypt_counts.csv')
PLOTS_DIR = os.path.join(ANALYSIS_DIR, 'plots')
os.makedirs(PLOTS_DIR, exist_ok=True)

MODEL_DIRS = {
    'node2d':   os.path.join(BASE, '16535750_node2d_2026-03-31_08-56-51'),
    'vertex2d': os.path.join(BASE, '16535751_vertex2d_2026-03-31_08-56-51'),
}

MODEL_LABELS = {
    'node2d':   '2D Node-based',
    'vertex2d': '2D Vertex-based',
}

MODEL_COLORS = {
    'node2d':   '#1f77b4',
    'vertex2d': '#ff7f0e',
}

STIFFNESSES = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 35.0, 50.0, 70.0, 100.0]
N_RUNS = 10


# ===========================================================================
# Data loading
# ===========================================================================
def load_crypt_counts():
    """Load crypt count data from CSV."""
    data = defaultdict(lambda: defaultdict(list))
    circ_data = defaultdict(lambda: defaultdict(list))

    with open(CRYPT_CSV) as f:
        reader = csv.DictReader(f)
        for row in reader:
            model = row['model']
            stiffness = float(row['stiffness'])
            data[model][stiffness].append(int(row['num_crypts']))
            circ_val = row['circularity']
            if circ_val and circ_val != 'nan':
                circ_data[model][stiffness].append(float(circ_val))

    return data, circ_data


def load_timeseries_data():
    """Load time-series summary data from crypt_summary.csv files."""
    ts_data = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    # ts_data[model][stiffness][run] = list of dicts with time-series

    for model, model_dir in MODEL_DIRS.items():
        if not os.path.isdir(model_dir):
            continue
        for stiffness in STIFFNESSES:
            for run in range(N_RUNS):
                csv_path = os.path.join(model_dir, f's{stiffness}_r{run}',
                                        'results_from_time_0', 'crypt_summary.csv')
                if not os.path.exists(csv_path):
                    continue
                try:
                    with open(csv_path) as f:
                        reader = csv.DictReader(f)
                        rows = []
                        for row in reader:
                            rows.append({
                                'time': float(row['time']),
                                'num_cells': int(row['num_cells']),
                                'mean_r': float(row['mean_r']),
                                'var_r': float(row['var_r']),
                                'r_range': float(row['r_range']),
                                'max_vol': float(row['max_vol']),
                                'min_vol': float(row['min_vol']),
                            })
                        ts_data[model][stiffness][run] = rows
                except Exception:
                    pass

    return ts_data


def get_final_metrics(ts_data):
    """Extract final-timestep metrics from time-series data."""
    final = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    # final[model][stiffness][metric] = list of values

    for model in ts_data:
        for stiffness in ts_data[model]:
            for run in ts_data[model][stiffness]:
                rows = ts_data[model][stiffness][run]
                if not rows:
                    continue
                last = rows[-1]
                for key in ['num_cells', 'mean_r', 'var_r', 'r_range']:
                    final[model][stiffness][key].append(last[key])

    return final


# ===========================================================================
# Plotting helpers
# ===========================================================================
def plot_metric_vs_stiffness(metric_data, ylabel, title, filename,
                              log_x=True, ylim=None):
    """Create a box+scatter plot of a metric vs stiffness for both 2D models."""
    fig, ax = plt.subplots(figsize=(10, 6))

    x_positions = np.arange(len(STIFFNESSES))
    width = 0.35

    for i, model in enumerate(['node2d', 'vertex2d']):
        if model not in metric_data:
            continue

        positions = x_positions + (i - 0.5) * width
        box_data = []
        for s in STIFFNESSES:
            vals = metric_data[model].get(s, [])
            box_data.append(vals if vals else [0])

        bp = ax.boxplot(box_data, positions=positions, widths=width * 0.8,
                        patch_artist=True, showfliers=False,
                        medianprops=dict(color='black', linewidth=1.5))

        color = MODEL_COLORS[model]
        for patch in bp['boxes']:
            patch.set_facecolor(color)
            patch.set_alpha(0.5)

        # Scatter individual points
        for j, s in enumerate(STIFFNESSES):
            vals = metric_data[model].get(s, [])
            if vals:
                jitter = np.random.normal(0, width * 0.1, len(vals))
                ax.scatter([positions[j]] * len(vals) + jitter, vals,
                          c=color, s=15, alpha=0.6, zorder=3, edgecolors='none')

        # Mean line
        means = [np.mean(metric_data[model].get(s, [0])) for s in STIFFNESSES]
        ax.plot(positions, means, '-o', color=color, markersize=5, linewidth=1.5,
                label=MODEL_LABELS[model], zorder=4)

    ax.set_xticks(x_positions)
    ax.set_xticklabels([str(s) for s in STIFFNESSES], rotation=45)
    ax.set_xlabel('ECM Stiffness (CD$^{-2}$)')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    if ylim is not None:
        ax.set_ylim(ylim)

    plt.tight_layout()
    path = os.path.join(PLOTS_DIR, filename)
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")


def plot_timeseries(ts_data, metric, ylabel, title, filename,
                    stiffnesses_to_show=None):
    """Plot time-series of a metric for selected stiffnesses."""
    if stiffnesses_to_show is None:
        stiffnesses_to_show = [0.5, 5.0, 20.0, 50.0, 100.0]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

    for ax, model in zip(axes, ['node2d', 'vertex2d']):
        if model not in ts_data:
            ax.set_title(f'{MODEL_LABELS.get(model, model)} — no data')
            continue

        cmap = plt.cm.viridis
        norm = plt.Normalize(vmin=np.log10(min(stiffnesses_to_show)),
                             vmax=np.log10(max(stiffnesses_to_show)))

        for stiffness in stiffnesses_to_show:
            color = cmap(norm(np.log10(stiffness)))
            for run in ts_data[model].get(stiffness, {}):
                rows = ts_data[model][stiffness][run]
                if not rows:
                    continue
                times = [r['time'] for r in rows]
                vals = [r[metric] for r in rows]
                ax.plot(times, vals, color=color, alpha=0.3, linewidth=0.8)

            # Mean line
            all_runs = ts_data[model].get(stiffness, {})
            if all_runs:
                # Align to common time points
                first_run = list(all_runs.values())[0]
                times = [r['time'] for r in first_run]
                all_vals = []
                for run_rows in all_runs.values():
                    if len(run_rows) == len(times):
                        all_vals.append([r[metric] for r in run_rows])
                if all_vals:
                    mean_vals = np.mean(all_vals, axis=0)
                    ax.plot(times, mean_vals, color=color, linewidth=2.5,
                           label=f's={stiffness}')

        ax.set_xlabel('Time (h)')
        ax.set_ylabel(ylabel)
        ax.set_title(f'{MODEL_LABELS[model]}')
        ax.legend(fontsize=8, loc='upper left')
        ax.grid(True, alpha=0.3)

    fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    path = os.path.join(PLOTS_DIR, filename)
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")


def plot_combined_summary(crypt_data, circ_data, final_metrics):
    """Create a combined multi-panel summary figure."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # Helper: plot metric vs stiffness on a given axis
    def _plot_on_ax(ax, metric_data, ylabel, title):
        for model in ['node2d', 'vertex2d']:
            if model not in metric_data:
                continue
            means = []
            stds = []
            for s in STIFFNESSES:
                vals = metric_data[model].get(s, [0])
                means.append(np.mean(vals))
                stds.append(np.std(vals))
            color = MODEL_COLORS[model]
            ax.errorbar(STIFFNESSES, means, yerr=stds, fmt='-o',
                       color=color, markersize=4, linewidth=1.5, capsize=3,
                       label=MODEL_LABELS[model])
        ax.set_xlabel('ECM Stiffness (CD$^{-2}$)')
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xscale('log')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    # Panel 1: Crypt count
    _plot_on_ax(axes[0, 0], crypt_data, 'Number of Crypts',
                'Crypt Count vs ECM Stiffness')

    # Panel 2: Circularity
    _plot_on_ax(axes[0, 1], circ_data, 'Circularity',
                'Circularity vs ECM Stiffness')

    # Panel 3: Number of cells
    cell_data = defaultdict(dict)
    for model in final_metrics:
        for s in final_metrics[model]:
            cell_data[model][s] = final_metrics[model][s].get('num_cells', [])
    _plot_on_ax(axes[0, 2], cell_data, 'Number of Cells',
                'Cell Count vs ECM Stiffness')

    # Panel 4: Mean radius
    r_data = defaultdict(dict)
    for model in final_metrics:
        for s in final_metrics[model]:
            r_data[model][s] = final_metrics[model][s].get('mean_r', [])
    _plot_on_ax(axes[1, 0], r_data, 'Mean Radius (CD)',
                'Mean Radius vs ECM Stiffness')

    # Panel 5: Radius variance
    var_data = defaultdict(dict)
    for model in final_metrics:
        for s in final_metrics[model]:
            var_data[model][s] = final_metrics[model][s].get('var_r', [])
    _plot_on_ax(axes[1, 1], var_data, 'Radius Variance (CD²)',
                'Radius Variance vs ECM Stiffness')

    # Panel 6: Radius range
    range_data = defaultdict(dict)
    for model in final_metrics:
        for s in final_metrics[model]:
            range_data[model][s] = final_metrics[model][s].get('r_range', [])
    _plot_on_ax(axes[1, 2], range_data, 'Radius Range (CD)',
                'Radius Range vs ECM Stiffness')

    fig.suptitle('ECM Stiffness Sweep — 2D Model Comparison\n(Node2d vs Vertex2d)',
                 fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    path = os.path.join(PLOTS_DIR, 'comparison_2d_summary.png')
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")


def plot_crypt_count_detailed(crypt_data):
    """Dedicated high-quality crypt count comparison plot."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: Linear scale with boxplots
    ax = axes[0]
    x_positions = np.arange(len(STIFFNESSES))
    width = 0.35

    for i, model in enumerate(['node2d', 'vertex2d']):
        if model not in crypt_data:
            continue
        positions = x_positions + (i - 0.5) * width
        box_data = [crypt_data[model].get(s, [0]) for s in STIFFNESSES]

        bp = ax.boxplot(box_data, positions=positions, widths=width * 0.8,
                        patch_artist=True, showfliers=True,
                        medianprops=dict(color='black', linewidth=1.5))

        color = MODEL_COLORS[model]
        for patch in bp['boxes']:
            patch.set_facecolor(color)
            patch.set_alpha(0.5)

        means = [np.mean(crypt_data[model].get(s, [0])) for s in STIFFNESSES]
        ax.plot(positions, means, 'D', color=color, markersize=5, zorder=5,
                label=f'{MODEL_LABELS[model]} (mean)')

    ax.set_xticks(x_positions)
    ax.set_xticklabels([str(s) for s in STIFFNESSES], rotation=45)
    ax.set_xlabel('ECM Stiffness (CD$^{-2}$)')
    ax.set_ylabel('Number of Crypts')
    ax.set_title('Crypt Count — Boxplot')
    ax.legend(fontsize=8)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.grid(True, alpha=0.3, axis='y')

    # Right: Log-scale mean ± std
    ax = axes[1]
    for model in ['node2d', 'vertex2d']:
        if model not in crypt_data:
            continue
        means = []
        stds = []
        for s in STIFFNESSES:
            vals = crypt_data[model].get(s, [0])
            means.append(np.mean(vals))
            stds.append(np.std(vals))
        color = MODEL_COLORS[model]
        ax.errorbar(STIFFNESSES, means, yerr=stds, fmt='-o',
                   color=color, markersize=5, linewidth=2, capsize=4,
                   label=MODEL_LABELS[model])

    ax.set_xscale('log')
    ax.set_xlabel('ECM Stiffness (CD$^{-2}$)')
    ax.set_ylabel('Number of Crypts (mean ± std)')
    ax.set_title('Crypt Count — Log Scale')
    ax.legend()
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.grid(True, alpha=0.3)

    fig.suptitle('Crypt Count vs ECM Stiffness — 2D Models',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    path = os.path.join(PLOTS_DIR, 'comparison_2d_crypt_count.png')
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")


def plot_crypt_count_and_circularity(crypt_data, circ_data):
    """Combined crypt count + circularity on shared x-axis."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    for model in ['node2d', 'vertex2d']:
        color = MODEL_COLORS[model]
        label = MODEL_LABELS[model]

        if model in crypt_data:
            means = [np.mean(crypt_data[model].get(s, [0])) for s in STIFFNESSES]
            stds = [np.std(crypt_data[model].get(s, [0])) for s in STIFFNESSES]
            ax1.errorbar(STIFFNESSES, means, yerr=stds, fmt='-o', color=color,
                        markersize=5, linewidth=2, capsize=4, label=label)

        if model in circ_data:
            means = [np.mean(circ_data[model].get(s, [0.5])) for s in STIFFNESSES]
            stds = [np.std(circ_data[model].get(s, [0.5])) for s in STIFFNESSES]
            ax2.errorbar(STIFFNESSES, means, yerr=stds, fmt='-o', color=color,
                        markersize=5, linewidth=2, capsize=4, label=label)

    ax1.set_ylabel('Number of Crypts')
    ax1.set_title('Crypt Count vs ECM Stiffness')
    ax1.legend()
    ax1.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.grid(True, alpha=0.3)

    ax2.set_xlabel('ECM Stiffness (CD$^{-2}$)')
    ax2.set_ylabel('Circularity')
    ax2.set_title('Circularity vs ECM Stiffness')
    ax2.set_xscale('log')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    fig.suptitle('Crypt Morphology — 2D Model Comparison',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    path = os.path.join(PLOTS_DIR, 'comparison_2d_crypt_and_circularity.png')
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")


# ===========================================================================
# Main
# ===========================================================================
def main():
    print("=" * 60)
    print("Generating 2D comparison plots")
    print("=" * 60)

    # Load data
    print("\nLoading crypt count data...")
    crypt_data, circ_data = load_crypt_counts()

    print("Loading time-series data...")
    ts_data = load_timeseries_data()
    final_metrics = get_final_metrics(ts_data)

    # Generate individual metric plots
    print("\nGenerating plots...")

    plot_crypt_count_detailed(crypt_data)
    plot_crypt_count_and_circularity(crypt_data, circ_data)

    # Metrics from final timestep
    for metric, ylabel, title, fname in [
        ('num_cells', 'Number of Cells', 'Cell Count vs ECM Stiffness',
         'comparison_2d_num_cells.png'),
        ('mean_r', 'Mean Radius (CD)', 'Mean Radius vs ECM Stiffness',
         'comparison_2d_mean_r.png'),
        ('var_r', 'Radius Variance (CD²)', 'Radius Variance vs ECM Stiffness',
         'comparison_2d_var_r.png'),
        ('r_range', 'Radius Range (CD)', 'Radius Range vs ECM Stiffness',
         'comparison_2d_r_range.png'),
    ]:
        metric_dict = defaultdict(dict)
        for model in final_metrics:
            for s in final_metrics[model]:
                metric_dict[model][s] = final_metrics[model][s].get(metric, [])
        plot_metric_vs_stiffness(metric_dict, ylabel, title, fname)

    # Time-series plots
    for metric, ylabel, title, fname in [
        ('num_cells', 'Number of Cells', 'Cell Growth Over Time',
         'comparison_2d_timeseries_cells.png'),
        ('mean_r', 'Mean Radius (CD)', 'Mean Radius Over Time',
         'comparison_2d_timeseries_mean_r.png'),
        ('var_r', 'Radius Variance (CD²)', 'Radius Variance Over Time',
         'comparison_2d_timeseries_var_r.png'),
    ]:
        plot_timeseries(ts_data, metric, ylabel, title, fname)

    # Combined summary figure
    plot_combined_summary(crypt_data, circ_data, final_metrics)

    print(f"\nAll plots saved to: {PLOTS_DIR}")
    print("Done!")


if __name__ == '__main__':
    main()
