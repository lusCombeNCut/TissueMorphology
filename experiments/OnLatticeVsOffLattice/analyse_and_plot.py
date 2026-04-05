#!/usr/bin/env python3
"""
analyse_and_plot.py — On-Lattice vs Off-Lattice ECM Comparison

Compares grid-based (on-lattice) ECMConfinementForce with particle-based
(off-lattice) GhostNodeEcmForce across matched stiffness levels.

Key questions:
  - Do both ECM representations produce similar crypt budding behaviour?
  - How does ECM discretisation affect morphological outcomes?
  - Are there stiffness regimes where the models diverge?

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
from analysis_utils import _save_fig

PARAM_NAME = 'ECM Stiffness'
PARAM_UNIT = ''


def load_lattice_comparison(base_dir, model_filter=None):
    """
    Load data organised by ECM type and stiffness.
    Returns: dict[model] -> dict[ecm_type] -> dict[stiffness] -> list of runs
    """
    csvs = find_summary_csvs(base_dir)
    if not csvs:
        return {}

    # result[model][ecm_type][stiffness] = [(run_num, data), ...]
    result = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    for csv_path, params_path in csvs:
        params = load_params_json(params_path) if params_path else None
        model = identify_model_type(csv_path, params)
        if model_filter and model != model_filter:
            continue

        # Determine ECM type
        ecm_type = None
        if params:
            ecm_type = params.get('simulation', {}).get('ecmType')
        if ecm_type is None:
            # Infer from directory name
            if 'on-lattice' in csv_path:
                ecm_type = 'on-lattice'
            elif 'off-lattice' in csv_path:
                ecm_type = 'off-lattice'
        if ecm_type is None:
            # Infer from feature flags
            if params:
                ghost = params.get('features', {}).get('enableGhostNodeECM', False)
                ecm_type = 'off-lattice' if ghost else 'on-lattice'
        if ecm_type is None:
            continue

        # Get stiffness
        stiffness = None
        if params:
            stiffness = extract_param_from_json(params,
                            'forces.ECMConfinementForce.ecmStiffness')
        if stiffness is None:
            m = re.search(r's([\d.]+)_r(\d+)', csv_path)
            if m:
                stiffness = float(m.group(1))
        if stiffness is None:
            continue

        run_num = params.get('simulation', {}).get('runNumber', 0) if params else 0
        data = load_crypt_summary(csv_path)
        result[model][ecm_type][stiffness].append((run_num, data))

    return dict(result)


def plot_comparison(model_data, model_type, plots_dir):
    """Generate comparison plots between on-lattice and off-lattice ECM."""
    if not HAS_MATPLOTLIB:
        return

    setup_style()

    ecm_styles = {
        'on-lattice':  ('steelblue', 'o', '-',  'Grid-based (on-lattice)'),
        'off-lattice': ('coral',     's', '--', 'Particle-based (off-lattice)'),
    }

    # 1. Mean ± SD summary for each metric
    for col, ylabel, title_suffix in [
        ('num_cells', 'Final Cell Count', 'Cell Count'),
        ('mean_r', 'Final Mean Radius (CD)', 'Organoid Size'),
        ('var_r', 'Final Radial Variance', 'Radial Variance'),
        ('r_range', 'Final Radial Range (CD)', 'Budding Amplitude'),
    ]:
        fig, ax = plt.subplots(figsize=figsize())

        for ecm_type, (colour, marker, ls, label) in ecm_styles.items():
            if ecm_type not in model_data:
                continue
            sweep = model_data[ecm_type]
            stiffnesses = sorted(sweep.keys())
            means = []
            stds = []
            for s in stiffnesses:
                vals = [d[col][-1] for _, d in sweep[s]
                        if col in d and len(d[col]) > 0]
                means.append(np.mean(vals) if vals else np.nan)
                stds.append(np.std(vals) if vals else np.nan)

            ax.errorbar(stiffnesses, means, yerr=stds,
                        fmt=f'{marker}{ls}', color=colour, capsize=5,
                        markersize=8, linewidth=2, label=label)

        ax.set_xlabel(f'{PARAM_NAME}')
        ax.set_ylabel(ylabel)
        ax.set_title(f'{model_type}: {title_suffix} — On-Lattice vs Off-Lattice')
        ax.set_xscale('log')
        _add_legend(ax, title='ECM Type', ncol=2)
        plt.tight_layout()
        path = os.path.join(plots_dir, f'{model_type}_{col}_lattice_comparison.svg')
        _save_fig(path)
        plt.close()
        print(f"  Saved: {path}")

    # 2. Time-series comparison at selected stiffness values
    all_stiffnesses = set()
    for ecm_type_data in model_data.values():
        all_stiffnesses.update(ecm_type_data.keys())
    all_stiffnesses = sorted(all_stiffnesses)

    # Pick 3 representative stiffness values
    if len(all_stiffnesses) >= 3:
        selected = [all_stiffnesses[0],
                    all_stiffnesses[len(all_stiffnesses) // 2],
                    all_stiffnesses[-1]]
    else:
        selected = all_stiffnesses

    for col, ylabel in [('mean_r', 'Mean Radius (CD)'),
                        ('r_range', 'Radial Range (CD)')]:
        fig, axes = plt.subplots(1, len(selected), figsize=figsize(len(selected), 1),
                                 squeeze=False)
        for idx, stiff in enumerate(selected):
            ax = axes[0, idx]
            for ecm_type, (colour, marker, ls, label) in ecm_styles.items():
                if ecm_type not in model_data:
                    continue
                if stiff not in model_data[ecm_type]:
                    continue
                runs = model_data[ecm_type][stiff]
                times, mean, std, n = aggregate_timeseries(runs, col)
                if len(times) > 0:
                    ax.plot(times, mean, color=colour, linestyle=ls, label=label)
                    ax.fill_between(times, mean - std, mean + std,
                                    color=colour, alpha=0.15)
            ax.set_xlabel('Time (hours)')
            ax.set_ylabel(ylabel)
            ax.set_title(f'Stiffness = {stiff}')
            if idx == 0:
                _add_legend(ax, title='ECM Type', ncol=1)

        plt.suptitle(f'{model_type}: {ylabel} Over Time — Lattice Comparison',
                     y=1.02)
        plt.tight_layout()
        path = os.path.join(plots_dir,
                            f'{model_type}_{col}_timeseries_lattice_comparison.svg')
        _save_fig(path)
        plt.close()
        print(f"  Saved: {path}")

    # 3. Difference plot: off-lattice minus on-lattice
    if 'on-lattice' in model_data and 'off-lattice' in model_data:
        fig, ax = plt.subplots(figsize=figsize())
        common_stiff = sorted(set(model_data['on-lattice'].keys()) &
                              set(model_data['off-lattice'].keys()))
        for col, colour, label in [
            ('var_r', 'steelblue', 'Radial Variance'),
            ('r_range', 'coral', 'Radial Range'),
        ]:
            diffs = []
            for s in common_stiff:
                on_vals = [d[col][-1] for _, d in model_data['on-lattice'][s]
                           if col in d and len(d[col]) > 0]
                off_vals = [d[col][-1] for _, d in model_data['off-lattice'][s]
                            if col in d and len(d[col]) > 0]
                if on_vals and off_vals:
                    diffs.append(np.mean(off_vals) - np.mean(on_vals))
                else:
                    diffs.append(np.nan)

            ax.plot(common_stiff, diffs, 'o-', color=colour, label=label,
                    markersize=8, linewidth=2)

        ax.axhline(0, color='gray', linewidth=0.5)
        ax.set_xlabel(PARAM_NAME)
        ax.set_ylabel('Off-Lattice minus On-Lattice')
        ax.set_title(f'{model_type}: ECM Representation Difference')
        ax.set_xscale('log')
        _add_legend(ax, title='Metric', ncol=2)
        plt.tight_layout()
        path = os.path.join(plots_dir, f'{model_type}_lattice_difference.svg')
        _save_fig(path)
        plt.close()
        print(f"  Saved: {path}")


def analyse_model(base_dir, model_type, plots_dir):
    print(f"\n--- Analysing {model_type} (On-Lattice vs Off-Lattice) ---")
    all_data = load_lattice_comparison(base_dir, model_filter=model_type)
    if model_type not in all_data:
        print(f"  No data found for {model_type}")
        return
    model_data = all_data[model_type]
    for ecm_type, sweep in model_data.items():
        n_runs = sum(len(v) for v in sweep.values())
        print(f"  {ecm_type}: {n_runs} runs across {len(sweep)} stiffness values")

    plot_comparison(model_data, model_type, plots_dir)

    # Also generate per-ECM-type standard plots
    for ecm_type, sweep in model_data.items():
        label = f'{model_type}_{ecm_type.replace("-", "_")}'
        generate_standard_plots(sweep, PARAM_NAME, PARAM_UNIT, label,
                                plots_dir, logx=True)


def main():
    parser = argparse.ArgumentParser(
        description='Analyse on-lattice vs off-lattice ECM comparison')
    parser.add_argument('--data-dir', required=True)
    parser.add_argument('--model', default='all',
                        choices=['node2d', 'vertex2d', 'node3d', 'vertex3d', 'all'])
    parser.add_argument('--output-dir', '-o', default=None)
    add_common_args(parser)
    args = parser.parse_args()
    apply_common_args(args)

    plots_dir = args.output_dir or os.path.join(args.data_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    models = ['node2d', 'vertex2d', 'node3d', 'vertex3d'] \
             if args.model == 'all' else [args.model]
    for model in models:
        analyse_model(args.data_dir, model, plots_dir)

    print(f"\nAll plots saved to: {plots_dir}")


if __name__ == '__main__':
    main()
