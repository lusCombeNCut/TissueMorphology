#!/usr/bin/env python3
"""
analyse_and_plot.py — ECM Stiffness Sweep

Generates comprehensive summary plots from ECM stiffness sweep simulation
output. Reads crypt_summary.csv for time-series data, VTU files for
cell-level properties, and ghost node VTP files for ECM statistics.

Produces a plots/ subfolder with all figures.

Usage:
  python analyse_and_plot.py --data-dir /path/to/sim_output/<RUN_TAG>
  python analyse_and_plot.py --data-dir /path/to/merged/CryptBudding --model node2d
  python analyse_and_plot.py --data-dir /path/to/merged --model all --full
"""

import os
import sys
import argparse

# Add parent directory for shared utilities
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_utils import *

# Import crypt counting if available
try:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from simple_crypt_count import count_crypts_simple_method, load_final_outline
    HAS_CRYPT_COUNT = True
except ImportError:
    HAS_CRYPT_COUNT = False
    print("Note: simple_crypt_count not available — skipping crypt counting.")


PARAM_NAME = 'ECM Stiffness'
PARAM_UNIT = 'CD\u207b\u00b2'  # per cell diameter squared
PARAM_JSON_PATH = 'forces.ECMConfinementForce.ecmStiffness'


def run_crypt_counting(base_dir, model_type, plots_dir):
    """
    Run crypt counting on final-timestep output for 2D models.
    Returns dict[stiffness] -> list of (n_crypts, circularity).
    """
    if not HAS_CRYPT_COUNT:
        return {}
    if model_type not in ('node2d', 'vertex2d'):
        print(f"  Crypt counting not supported for {model_type} (3D)")
        return {}

    csvs = find_summary_csvs(base_dir)
    crypt_results = defaultdict(list)

    for csv_path, params_path in csvs:
        params = load_params_json(params_path) if params_path else None
        model = identify_model_type(csv_path, params)
        if model != model_type:
            continue

        stiffness = extract_param_from_json(params, PARAM_JSON_PATH) if params else None
        if stiffness is None:
            stiffness, _ = extract_param_from_dir(csv_path, 'stiffness')
        if stiffness is None:
            continue

        # Find the results_from_time_* directory
        data_dir = os.path.dirname(csv_path)

        try:
            if model_type == 'node2d':
                boundary, cell_types = load_final_outline(data_dir)
                result = count_crypts_simple_method(boundary,
                             boundary_is_ordered=True)
            else:
                # vertex2d: load from VTU
                from simple_crypt_count import load_boundary_from_vtu
                vtu_files = sorted(glob.glob(os.path.join(data_dir, 'results_*.vtu')))
                if not vtu_files:
                    continue
                boundary, cell_types = load_boundary_from_vtu(vtu_files[-1])
                result = count_crypts_simple_method(boundary,
                             boundary_is_ordered=True)

            crypt_results[stiffness].append((result.num_crypts, result.circularity))
        except Exception as e:
            print(f"  Crypt count error (stiffness={stiffness}): {e}")
            continue

    return dict(crypt_results)


def plot_crypt_counts(crypt_results, plots_dir, model_type):
    """Plot crypt count and circularity vs stiffness."""
    if not HAS_MATPLOTLIB or not crypt_results:
        return

    setup_style()
    stiffnesses = sorted(crypt_results.keys())

    # Crypt count box plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Box plot
    ax = axes[0]
    crypt_data = [[nc for nc, _ in crypt_results[s]] for s in stiffnesses]
    colours = get_param_colours(stiffnesses)

    bp = ax.boxplot(crypt_data, positions=range(len(stiffnesses)),
                    widths=0.5, patch_artist=True,
                    medianprops=dict(color='black', linewidth=2))
    for i, patch in enumerate(bp['boxes']):
        patch.set_facecolor(colours[stiffnesses[i]])
        patch.set_alpha(0.4)
    for i, (s, data) in enumerate(zip(stiffnesses, crypt_data)):
        if data:
            jitter = np.random.normal(0, 0.08, len(data))
            ax.scatter(np.full(len(data), i) + jitter, data,
                       color=colours[s], alpha=0.6, s=30, zorder=5)
    ax.set_xticks(range(len(stiffnesses)))
    ax.set_xticklabels([f'{s:.1f}' for s in stiffnesses], rotation=45)
    ax.set_xlabel(f'{PARAM_NAME} ({PARAM_UNIT})')
    ax.set_ylabel('Number of Crypts')
    ax.set_title(f'{model_type}: Crypt Count vs {PARAM_NAME}')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    # Circularity
    ax2 = axes[1]
    circ_data = [[c for _, c in crypt_results[s]] for s in stiffnesses]
    means = [np.mean(d) if d else 0 for d in circ_data]
    stds = [np.std(d) if d else 0 for d in circ_data]
    ax2.errorbar(stiffnesses, means, yerr=stds, fmt='o-', capsize=5,
                 color='steelblue', markersize=8, linewidth=2)
    ax2.set_xlabel(f'{PARAM_NAME} ({PARAM_UNIT})')
    ax2.set_ylabel('Circularity')
    ax2.set_title(f'{model_type}: Circularity vs {PARAM_NAME}')
    ax2.set_xscale('log')
    ax2.set_ylim(0, 1.05)

    plt.tight_layout()
    path = os.path.join(plots_dir, f'{model_type}_crypt_count_and_circularity.png')
    plt.savefig(path)
    plt.close()
    print(f"  Saved: {path}")


def analyse_model(base_dir, model_type, plots_dir, full=False):
    """Run full analysis for one model type."""
    print(f"\n--- Analysing {model_type} ---")

    sweep = load_sweep_data(base_dir, 'stiffness',
                            param_json_path=PARAM_JSON_PATH,
                            model_filter=model_type)
    if model_type not in sweep:
        print(f"  No data found for {model_type}")
        return

    sweep_data = sweep[model_type]
    print(f"  Found {sum(len(v) for v in sweep_data.values())} runs across "
          f"{len(sweep_data)} stiffness values")

    # Standard time-series and summary plots
    generate_standard_plots(sweep_data, PARAM_NAME, PARAM_UNIT, model_type,
                            plots_dir, logx=True)

    # Crypt counting (2D only)
    crypt_results = run_crypt_counting(base_dir, model_type, plots_dir)
    if crypt_results:
        plot_crypt_counts(crypt_results, plots_dir, model_type)

    # Full analysis: parse VTU/VTP files (slower)
    if full:
        print(f"  Running full VTU/VTP analysis for {model_type}...")
        analyse_vtu_data(base_dir, model_type, sweep_data, plots_dir)


def analyse_vtu_data(base_dir, model_type, sweep_data, plots_dir):
    """
    Extended analysis from VTU and VTP files. Produces contact inhibition,
    cell type ratios, lumen force, and ghost node plots.
    """
    if not HAS_MATPLOTLIB:
        return

    setup_style()
    param_vals = sorted(sweep_data.keys())
    colours = get_param_colours(param_vals)

    # For each parameter value, pick first replicate's data directory
    # to get VTU/VTP files
    first_run_dirs = {}
    csvs = find_summary_csvs(base_dir)
    for csv_path, params_path in csvs:
        params = load_params_json(params_path) if params_path else None
        model = identify_model_type(csv_path, params)
        if model != model_type:
            continue
        stiffness = extract_param_from_json(params, PARAM_JSON_PATH) if params else None
        if stiffness is None:
            continue
        if stiffness not in first_run_dirs:
            first_run_dirs[stiffness] = os.path.dirname(csv_path)

    # Contact inhibition
    fig, ax = plt.subplots(figsize=(10, 5))
    has_ci_data = False
    for pval in param_vals:
        if pval not in first_run_dirs:
            continue
        data_dir = first_run_dirs[pval]
        times, fracs = compute_contact_inhibition_timeseries(data_dir)
        if len(times) > 0:
            has_ci_data = True
            ax.plot(times, fracs, color=colours[pval],
                    label=f'Stiffness={pval}')
    if has_ci_data:
        ax.set_xlabel('Timestep')
        ax.set_ylabel('Fraction Contact Inhibited')
        ax.set_title(f'{model_type}: Contact Inhibition Over Time')
        ax.legend(fontsize=8, ncol=2)
        plt.tight_layout()
        path = os.path.join(plots_dir, f'{model_type}_contact_inhibition.png')
        plt.savefig(path)
        print(f"  Saved: {path}")
    plt.close()

    # Cell type ratios (pick a representative stiffness)
    mid_idx = len(param_vals) // 2
    mid_stiffness = param_vals[mid_idx] if param_vals else None
    if mid_stiffness and mid_stiffness in first_run_dirs:
        data_dir = first_run_dirs[mid_stiffness]
        ct_data = compute_cell_type_timeseries(data_dir)
        if 'times' in ct_data and len(ct_data['times']) > 0:
            fig, ax = plt.subplots(figsize=(10, 5))
            type_colours = {'Stem': 'green', 'Transit': 'blue',
                            'Paneth': 'red', 'Enterocyte': 'purple'}
            for name, colour in type_colours.items():
                if name in ct_data:
                    ax.plot(ct_data['times'], ct_data[name], color=colour,
                            label=name)
            ax.set_xlabel('Timestep')
            ax.set_ylabel('Cell Type Fraction')
            ax.set_title(f'{model_type}: Cell Type Ratios (stiffness={mid_stiffness})')
            ax.legend()
            ax.set_ylim(0, 1)
            plt.tight_layout()
            path = os.path.join(plots_dir, f'{model_type}_cell_type_ratios.png')
            plt.savefig(path)
            plt.close()
            print(f"  Saved: {path}")

    # Lumen force (node models primarily)
    fig, ax = plt.subplots(figsize=(10, 5))
    has_lumen = False
    for pval in param_vals:
        if pval not in first_run_dirs:
            continue
        data_dir = first_run_dirs[pval]
        times, forces = compute_lumen_force_timeseries(data_dir)
        if len(times) > 0:
            has_lumen = True
            ax.plot(times, forces, color=colours[pval],
                    label=f'Stiffness={pval}')
    if has_lumen:
        ax.set_xlabel('Timestep')
        ax.set_ylabel('Mean Lumen Force Magnitude')
        ax.set_title(f'{model_type}: Mean Lumen Force Over Time')
        ax.legend(fontsize=8, ncol=2)
        plt.tight_layout()
        path = os.path.join(plots_dir, f'{model_type}_lumen_force.png')
        plt.savefig(path)
        print(f"  Saved: {path}")
    plt.close()

    # Ghost node ECM statistics
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    has_ghost = False
    for pval in param_vals:
        if pval not in first_run_dirs:
            continue
        data_dir = first_run_dirs[pval]
        ghost_ts = get_ghost_node_timeseries(data_dir)
        if 'times' not in ghost_ts or len(ghost_ts['times']) == 0:
            continue
        has_ghost = True
        t = ghost_ts['times']
        c = colours[pval]
        lbl = f'S={pval}'

        axes[0].plot(t, ghost_ts['n_nodes'], color=c, label=lbl)
        axes[1].plot(t, ghost_ts['mean_density'], color=c, label=lbl)
        if 'mean_strain' in ghost_ts:
            axes[2].plot(t, ghost_ts['mean_strain'], color=c, label=lbl)

    if has_ghost:
        axes[0].set_xlabel('Timestep')
        axes[0].set_ylabel('Active Ghost Nodes')
        axes[0].set_title('ECM Node Count')
        axes[0].legend(fontsize=7, ncol=2)
        axes[1].set_xlabel('Timestep')
        axes[1].set_ylabel('Mean ECM Density')
        axes[1].set_title('ECM Density')
        axes[1].legend(fontsize=7, ncol=2)
        axes[2].set_xlabel('Timestep')
        axes[2].set_ylabel('Mean Rest Length Strain')
        axes[2].set_title('Viscoelastic Strain')
        axes[2].legend(fontsize=7, ncol=2)
        plt.tight_layout()
        path = os.path.join(plots_dir, f'{model_type}_ghost_node_ecm.png')
        plt.savefig(path)
        print(f"  Saved: {path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Analyse ECM stiffness sweep and generate plots')
    parser.add_argument('--data-dir', required=True,
                        help='Root directory of sweep output')
    parser.add_argument('--model', default='all',
                        choices=['node2d', 'vertex2d', 'node3d', 'vertex3d', 'all'],
                        help='Model type to analyse (default: all)')
    parser.add_argument('--output-dir', '-o', default=None,
                        help='Output directory for plots (default: plots/ in data-dir)')
    parser.add_argument('--full', action='store_true',
                        help='Run full analysis including VTU/VTP parsing (slower)')

    args = parser.parse_args()

    plots_dir = args.output_dir or os.path.join(args.data_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    models = ['node2d', 'vertex2d', 'node3d', 'vertex3d'] \
             if args.model == 'all' else [args.model]

    for model in models:
        analyse_model(args.data_dir, model, plots_dir, full=args.full)

    # Cross-model comparison plot
    if len(models) > 1:
        print("\n--- Cross-model comparison ---")
        all_sweep = load_sweep_data(args.data_dir, 'stiffness',
                                    param_json_path=PARAM_JSON_PATH)
        if len(all_sweep) > 1:
            for col, ylabel, title in [
                ('num_cells', 'Final Cell Count', 'Cell Count vs ECM Stiffness'),
                ('mean_r', 'Final Mean Radius (CD)', 'Organoid Size vs ECM Stiffness'),
                ('var_r', 'Final Radial Variance', 'Radial Variance vs ECM Stiffness'),
                ('r_range', 'Final Radial Range (CD)', 'Budding Amplitude vs ECM Stiffness'),
            ]:
                path = os.path.join(plots_dir, f'comparison_{col}_vs_stiffness.png')
                plot_multi_model_comparison(all_sweep, col, PARAM_NAME, PARAM_UNIT,
                                           ylabel, title, path, logx=True)

    print(f"\nAll plots saved to: {plots_dir}")


if __name__ == '__main__':
    main()
