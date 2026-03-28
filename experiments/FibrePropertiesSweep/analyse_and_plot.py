#!/usr/bin/env python3
"""
analyse_and_plot.py — Fibre Properties Sweep

Analyses how ECM fibre remodeling rate affects crypt budding morphology.
Higher remodeling rates allow fibres to realign more quickly with cell
traction forces, potentially facilitating or hindering bud formation.

Usage:
  python analyse_and_plot.py --data-dir /path/to/sim_output/<RUN_TAG>
  python analyse_and_plot.py --data-dir /path/to/merged --model node2d
"""

import os
import sys
import re
import argparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_utils import *

PARAM_NAME = 'Fibre Remodel Rate'
PARAM_UNIT = 'h\u207b\u00b9'
PARAM_JSON_PATH = 'forces.GhostNodeECM.fibreRemodelingRate'


def load_fibre_sweep(base_dir, model_filter=None):
    """Load data organised by fibre remodeling rate."""
    csvs = find_summary_csvs(base_dir)
    if not csvs:
        return {}

    sweep = defaultdict(lambda: defaultdict(list))

    for csv_path, params_path in csvs:
        params = load_params_json(params_path) if params_path else None
        model = identify_model_type(csv_path, params)
        if model_filter and model != model_filter:
            continue

        rate = None
        if params:
            rate = extract_param_from_json(params, PARAM_JSON_PATH)
        if rate is None:
            m = re.search(r'/f([\d.]+)_r(\d+)/', csv_path)
            if m:
                rate = float(m.group(1))
        if rate is None:
            continue

        run_num = params.get('simulation', {}).get('runNumber', 0) if params else 0
        data = load_crypt_summary(csv_path)
        sweep[model][rate].append((run_num, data))

    return dict(sweep)


def plot_fibre_specific(sweep_data, model_type, plots_dir):
    """Fibre-specific analysis plots."""
    if not HAS_MATPLOTLIB:
        return

    setup_style()
    param_vals = sorted(sweep_data.keys())
    colours = get_param_colours(param_vals)

    # Variance onset time: time at which var_r first exceeds 2x initial value
    fig, ax = plt.subplots(figsize=(8, 5))
    onset_times = []
    for pval in param_vals:
        runs = sweep_data[pval]
        times_list = []
        for _, data in runs:
            if 'var_r' in data and 'time' in data and len(data['var_r']) > 5:
                t = data['time']
                v = data['var_r']
                init_var = np.mean(v[:3])
                threshold = max(init_var * 3, 0.01)
                exceed_idx = np.where(v > threshold)[0]
                if len(exceed_idx) > 0:
                    times_list.append(t[exceed_idx[0]])
        if times_list:
            onset_times.append((pval, np.mean(times_list), np.std(times_list)))

    if onset_times:
        pvals = [x[0] for x in onset_times]
        means = [x[1] for x in onset_times]
        stds = [x[2] for x in onset_times]
        ax.errorbar(pvals, means, yerr=stds, fmt='o-', capsize=5,
                    color='steelblue', markersize=8, linewidth=2)
        ax.set_xlabel(f'{PARAM_NAME} ({PARAM_UNIT})')
        ax.set_ylabel('Budding Onset Time (hours)')
        ax.set_title(f'{model_type}: Budding Onset vs Fibre Remodeling Rate')
        plt.tight_layout()
        path = os.path.join(plots_dir, f'{model_type}_budding_onset.png')
        plt.savefig(path)
        print(f"  Saved: {path}")
    plt.close()

    # ECM anisotropy effect: var_r at final time vs fibre rate
    fig, ax = plt.subplots(figsize=(8, 5))
    for pval in param_vals:
        runs = sweep_data[pval]
        final_vars = [d['var_r'][-1] for _, d in runs
                      if 'var_r' in d and len(d['var_r']) > 0]
        final_ranges = [d['r_range'][-1] for _, d in runs
                        if 'r_range' in d and len(d['r_range']) > 0]
        if final_vars and final_ranges:
            ax.scatter([np.mean(final_vars)], [np.mean(final_ranges)],
                       color=colours[pval], s=100, zorder=5,
                       label=f'Rate={pval}')
            ax.errorbar([np.mean(final_vars)], [np.mean(final_ranges)],
                        xerr=[np.std(final_vars)], yerr=[np.std(final_ranges)],
                        color=colours[pval], capsize=3, linewidth=1)

    ax.set_xlabel('Final Radial Variance')
    ax.set_ylabel('Final Radial Range (CD)')
    ax.set_title(f'{model_type}: Morphological Deformation vs Fibre Remodeling')
    ax.legend(fontsize=8)
    plt.tight_layout()
    path = os.path.join(plots_dir, f'{model_type}_deformation_scatter.png')
    plt.savefig(path)
    plt.close()
    print(f"  Saved: {path}")


def analyse_model(base_dir, model_type, plots_dir, full=False):
    print(f"\n--- Analysing {model_type} (Fibre Properties Sweep) ---")
    sweep = load_fibre_sweep(base_dir, model_filter=model_type)
    if model_type not in sweep:
        print(f"  No data found for {model_type}")
        return
    sweep_data = sweep[model_type]
    print(f"  Found {sum(len(v) for v in sweep_data.values())} runs across "
          f"{len(sweep_data)} fibre rate values")

    generate_standard_plots(sweep_data, PARAM_NAME, PARAM_UNIT, model_type,
                            plots_dir, logx=False)
    plot_fibre_specific(sweep_data, model_type, plots_dir)


def main():
    parser = argparse.ArgumentParser(
        description='Analyse fibre properties sweep')
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
        all_sweep = load_fibre_sweep(args.data_dir)
        if len(all_sweep) > 1:
            for col, ylabel, title in [
                ('r_range', 'Final Radial Range', 'Budding Amplitude vs Fibre Rate'),
                ('var_r', 'Final Radial Variance', 'Radial Variance vs Fibre Rate'),
            ]:
                path = os.path.join(plots_dir, f'comparison_{col}_vs_fibre_rate.png')
                plot_multi_model_comparison(all_sweep, col, PARAM_NAME,
                    PARAM_UNIT, ylabel, title, path)

    print(f"\nAll plots saved to: {plots_dir}")


if __name__ == '__main__':
    main()
