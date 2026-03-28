#!/usr/bin/env python3
"""
analyse_and_plot.py — Lumen Pressure Sweep

Generates summary plots for lumen pressure sweep simulations.
Investigates how hydrostatic lumen pressure affects crypt budding,
organoid growth, and morphological instability.

Usage:
  python analyse_and_plot.py --data-dir /path/to/sim_output/<RUN_TAG>
  python analyse_and_plot.py --data-dir /path/to/merged --model node2d --full
"""

import os
import sys
import argparse
import glob

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_utils import *

try:
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), 'ECMStiffnessSweep'))
    from simple_crypt_count import count_crypts_simple_method, load_final_outline
    HAS_CRYPT_COUNT = True
except ImportError:
    HAS_CRYPT_COUNT = False

PARAM_NAME = 'Lumen Pressure'
PARAM_UNIT = ''
PARAM_JSON_PATH = 'forces.LumenPressureForce.lumenPressure'


def load_pressure_sweep(base_dir, model_filter=None):
    """Load data organised by lumen pressure rather than stiffness."""
    csvs = find_summary_csvs(base_dir)
    if not csvs:
        print(f"  No crypt_summary.csv found under {base_dir}")
        return {}

    sweep = defaultdict(lambda: defaultdict(list))

    for csv_path, params_path in csvs:
        params = load_params_json(params_path) if params_path else None
        model = identify_model_type(csv_path, params)
        if model_filter and model != model_filter:
            continue

        pressure = None
        if params:
            pressure = extract_param_from_json(params, PARAM_JSON_PATH)

        # Fall back to directory pattern p<val>_r<run>
        import re
        if pressure is None:
            m = re.search(r'/p([\d.]+)_r(\d+)/', csv_path)
            if m:
                pressure = float(m.group(1))

        if pressure is None:
            continue

        run_num = 0
        if params:
            run_num = params.get('simulation', {}).get('runNumber', 0)

        data = load_crypt_summary(csv_path)
        sweep[model][pressure].append((run_num, data))

    return dict(sweep)


def plot_pressure_specific(sweep_data, model_type, plots_dir):
    """Plots specific to lumen pressure analysis."""
    if not HAS_MATPLOTLIB:
        return

    setup_style()
    param_vals = sorted(sweep_data.keys())
    colours = get_param_colours(param_vals)

    # Growth rate: plot the derivative of mean_r
    fig, ax = plt.subplots(figsize=(10, 5))
    for pval in param_vals:
        runs = sweep_data[pval]
        times, mean, std, n = aggregate_timeseries(runs, 'mean_r')
        if len(times) > 2:
            dt = np.diff(times)
            growth_rate = np.diff(mean) / dt
            mid_times = 0.5 * (times[:-1] + times[1:])
            # Smooth
            if len(growth_rate) > 5:
                kernel = np.ones(5) / 5
                growth_rate = np.convolve(growth_rate, kernel, mode='same')
            ax.plot(mid_times, growth_rate, color=colours[pval],
                    label=f'P={pval}')
    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('Radial Growth Rate (CD/h)')
    ax.set_title(f'{model_type}: Radial Growth Rate vs Lumen Pressure')
    ax.legend(fontsize=8, ncol=2)
    ax.axhline(0, color='gray', linewidth=0.5)
    plt.tight_layout()
    path = os.path.join(plots_dir, f'{model_type}_growth_rate.png')
    plt.savefig(path)
    plt.close()
    print(f"  Saved: {path}")

    # Ratio of r_range to mean_r (normalised budding amplitude)
    fig, ax = plt.subplots(figsize=(10, 5))
    for pval in param_vals:
        runs = sweep_data[pval]
        # Compute ratio per run then aggregate
        ratio_runs = []
        for rn, data in runs:
            if 'r_range' in data and 'mean_r' in data:
                ratio = data['r_range'] / np.maximum(data['mean_r'], 1e-6)
                ratio_data = dict(data)
                ratio_data['ratio'] = ratio
                ratio_runs.append((rn, ratio_data))
        if ratio_runs:
            times, mean, std, n = aggregate_timeseries(ratio_runs, 'ratio')
            if len(times) > 0:
                ax.plot(times, mean, color=colours[pval], label=f'P={pval}')
                ax.fill_between(times, mean - std, mean + std,
                                color=colours[pval], alpha=0.15)
    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('Normalised Budding Amplitude (r_range / mean_r)')
    ax.set_title(f'{model_type}: Normalised Budding vs Lumen Pressure')
    ax.legend(fontsize=8, ncol=2)
    plt.tight_layout()
    path = os.path.join(plots_dir, f'{model_type}_normalised_budding.png')
    plt.savefig(path)
    plt.close()
    print(f"  Saved: {path}")


def analyse_model(base_dir, model_type, plots_dir, full=False):
    """Run analysis for one model type."""
    print(f"\n--- Analysing {model_type} (Lumen Pressure Sweep) ---")

    sweep = load_pressure_sweep(base_dir, model_filter=model_type)
    if model_type not in sweep:
        print(f"  No data found for {model_type}")
        return

    sweep_data = sweep[model_type]
    print(f"  Found {sum(len(v) for v in sweep_data.values())} runs across "
          f"{len(sweep_data)} pressure values")

    generate_standard_plots(sweep_data, PARAM_NAME, PARAM_UNIT, model_type,
                            plots_dir, logx=False)
    plot_pressure_specific(sweep_data, model_type, plots_dir)


def main():
    parser = argparse.ArgumentParser(
        description='Analyse lumen pressure sweep and generate plots')
    parser.add_argument('--data-dir', required=True)
    parser.add_argument('--model', default='all',
                        choices=['node2d', 'vertex2d', 'node3d', 'vertex3d', 'all'])
    parser.add_argument('--output-dir', '-o', default=None)
    parser.add_argument('--full', action='store_true')

    args = parser.parse_args()
    plots_dir = args.output_dir or os.path.join(args.data_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    models = ['node2d', 'vertex2d', 'node3d', 'vertex3d'] \
             if args.model == 'all' else [args.model]

    for model in models:
        analyse_model(args.data_dir, model, plots_dir, full=args.full)

    if len(models) > 1:
        print("\n--- Cross-model comparison ---")
        all_sweep = load_pressure_sweep(args.data_dir)
        if len(all_sweep) > 1:
            for col, ylabel, title in [
                ('num_cells', 'Final Cell Count', 'Cell Count vs Lumen Pressure'),
                ('mean_r', 'Final Mean Radius', 'Organoid Size vs Lumen Pressure'),
                ('r_range', 'Final Radial Range', 'Budding Amplitude vs Lumen Pressure'),
            ]:
                path = os.path.join(plots_dir, f'comparison_{col}_vs_pressure.png')
                plot_multi_model_comparison(all_sweep, col, PARAM_NAME, PARAM_UNIT,
                                           ylabel, title, path, logx=False)

    print(f"\nAll plots saved to: {plots_dir}")


if __name__ == '__main__':
    main()
