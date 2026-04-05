#!/usr/bin/env python3
"""
analyse_and_plot.py — Viscoelastic Properties Sweep

Analyses how ECM viscoelastic relaxation time (tau) affects crypt budding.
Shorter tau = more fluid-like ECM (rapid stress relaxation).
Longer tau = more elastic ECM (sustained resistance to deformation).

Includes ghost node VTP analysis for viscoelastic strain tracking.

Usage:
  python analyse_and_plot.py --data-dir /path/to/sim_output/<RUN_TAG>
  python analyse_and_plot.py --data-dir /path/to/merged --model node2d --full
"""

import os
import sys
import re
import argparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_utils import *
from analysis_utils import _save_fig

PARAM_NAME = 'Relaxation Time'
PARAM_UNIT = 'hours'
PARAM_JSON_PATH = 'forces.GhostNodeECM.ViscoelasticECM.ghostRelaxationTime'


def load_visco_sweep(base_dir, model_filter=None):
    """Load data organised by relaxation time."""
    csvs = find_summary_csvs(base_dir)
    if not csvs:
        return {}

    sweep = defaultdict(lambda: defaultdict(list))

    for csv_path, params_path in csvs:
        params = load_params_json(params_path) if params_path else None
        model = identify_model_type(csv_path, params)
        if model_filter and model != model_filter:
            continue

        tau = None
        if params:
            tau = extract_param_from_json(params, PARAM_JSON_PATH)
        if tau is None:
            m = re.search(r'/tau([\d.]+)_r(\d+)/', csv_path)
            if m:
                tau = float(m.group(1))
        if tau is None:
            continue

        run_num = params.get('simulation', {}).get('runNumber', 0) if params else 0
        data = load_crypt_summary(csv_path)
        sweep[model][tau].append((run_num, data))

    return dict(sweep)


def plot_visco_specific(sweep_data, model_type, plots_dir, base_dir):
    """Viscoelastic-specific analysis plots."""
    if not HAS_MATPLOTLIB:
        return

    setup_style()
    param_vals = sorted(sweep_data.keys())
    colours = get_param_colours(param_vals)

    # Phase diagram: final var_r vs tau (log scale)
    fig, ax = plt.subplots(figsize=figsize())
    for pval in param_vals:
        runs = sweep_data[pval]
        final_vars = [d['var_r'][-1] for _, d in runs
                      if 'var_r' in d and len(d['var_r']) > 0]
        if final_vars:
            ax.scatter([pval] * len(final_vars), final_vars,
                       color=colours[pval], alpha=0.5, s=40, zorder=5)
            ax.plot(pval, np.mean(final_vars), 'D', color=colours[pval],
                    markersize=10, zorder=6, markeredgecolor='black',
                    markeredgewidth=0.5)

    ax.set_xlabel(f'{PARAM_NAME} ({PARAM_UNIT})')
    ax.set_ylabel('Final Radial Variance')
    ax.set_title(f'{model_type}: Morphological Instability vs Relaxation Time')
    ax.set_xscale('log')
    plt.tight_layout()
    path = os.path.join(plots_dir, f'{model_type}_instability_vs_tau.svg')
    _save_fig(path)
    plt.close()
    print(f"  Saved: {path}")

    # Ghost node strain analysis (from VTP files)
    csvs = find_summary_csvs(base_dir)
    first_run_dirs = {}
    for csv_path, params_path in csvs:
        params = load_params_json(params_path) if params_path else None
        model = identify_model_type(csv_path, params)
        if model != model_type:
            continue
        tau = extract_param_from_json(params, PARAM_JSON_PATH) if params else None
        if tau is not None and tau not in first_run_dirs:
            first_run_dirs[tau] = os.path.dirname(csv_path)

    # Ghost node plots: strain, density, node count
    fig, axes = plt.subplots(2, 2, figsize=figsize(2, 2))
    has_data = False

    for pval in param_vals:
        if pval not in first_run_dirs:
            continue
        data_dir = first_run_dirs[pval]
        ghost_ts = get_ghost_node_timeseries(data_dir)
        if 'times' not in ghost_ts or len(ghost_ts['times']) == 0:
            continue
        has_data = True
        t = ghost_ts['times']
        c = colours[pval]
        lbl = f'{pval} {PARAM_UNIT}'

        axes[0, 0].plot(t, ghost_ts['n_nodes'], color=c, label=lbl)
        axes[0, 1].plot(t, ghost_ts['mean_density'], color=c, label=lbl)

        if 'mean_strain' in ghost_ts:
            axes[1, 0].plot(t, ghost_ts['mean_strain'], color=c, label=lbl)

        if 'mean_displacement' in ghost_ts:
            axes[1, 1].plot(t, ghost_ts['mean_displacement'], color=c, label=lbl)

    if has_data:
        axes[0, 0].set_ylabel('Active Ghost Nodes')
        axes[0, 0].set_title('ECM Node Count Over Time')

        axes[0, 1].set_ylabel('Mean ECM Density')
        axes[0, 1].set_title('ECM Degradation Over Time')

        axes[1, 0].set_xlabel('Timestep')
        axes[1, 0].set_ylabel('Mean Rest Length Strain')
        axes[1, 0].set_title('Viscoelastic Strain Relaxation')
        axes[1, 0].axhline(0, color='gray', linewidth=0.5)

        axes[1, 1].set_xlabel('Timestep')
        axes[1, 1].set_ylabel('Mean Node Displacement (CD)')
        axes[1, 1].set_title('ECM Node Displacement')

        for ax in axes.flat:
            ax.grid(alpha=0.3)

        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center',
                   bbox_to_anchor=(0.5, -0.05), ncol=3, title=PARAM_NAME)
        plt.tight_layout()
        path = os.path.join(plots_dir, f'{model_type}_ghost_node_viscoelastic.svg')
        _save_fig(path)
        print(f"  Saved: {path}")
    plt.close()


def analyse_model(base_dir, model_type, plots_dir, full=False):
    print(f"\n--- Analysing {model_type} (Viscoelastic Sweep) ---")
    sweep = load_visco_sweep(base_dir, model_filter=model_type)
    if model_type not in sweep:
        print(f"  No data found for {model_type}")
        return
    sweep_data = sweep[model_type]
    print(f"  Found {sum(len(v) for v in sweep_data.values())} runs across "
          f"{len(sweep_data)} tau values")

    generate_standard_plots(sweep_data, PARAM_NAME, PARAM_UNIT, model_type,
                            plots_dir, logx=True)
    plot_visco_specific(sweep_data, model_type, plots_dir, base_dir)


def main():
    parser = argparse.ArgumentParser(
        description='Analyse viscoelastic properties sweep')
    parser.add_argument('--data-dir', required=True)
    parser.add_argument('--model', default='all',
                        choices=['node2d', 'vertex2d', 'node3d', 'vertex3d', 'all'])
    parser.add_argument('--output-dir', '-o', default=None)
    parser.add_argument('--full', action='store_true')
    add_common_args(parser)
    args = parser.parse_args()
    apply_common_args(args)

    plots_dir = args.output_dir or os.path.join(args.data_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    models = ['node2d', 'vertex2d', 'node3d', 'vertex3d'] \
             if args.model == 'all' else [args.model]
    for model in models:
        analyse_model(args.data_dir, model, plots_dir, full=args.full)

    if len(models) > 1:
        all_sweep = load_visco_sweep(args.data_dir)
        if len(all_sweep) > 1:
            for col, ylabel, title in [
                ('r_range', 'Final Radial Range', 'Budding Amplitude vs Relaxation Time'),
                ('var_r', 'Final Radial Variance', 'Variance vs Relaxation Time'),
                ('num_cells', 'Final Cell Count', 'Cell Count vs Relaxation Time'),
            ]:
                path = os.path.join(plots_dir, f'comparison_{col}_vs_tau.svg')
                plot_multi_model_comparison(all_sweep, col, PARAM_NAME,
                    PARAM_UNIT, ylabel, title, path, logx=True)

    print(f"\nAll plots saved to: {plots_dir}")


if __name__ == '__main__':
    main()
