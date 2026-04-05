#!/usr/bin/env python3
import os
import sys
import argparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_utils import *
from analysis_utils import _save_fig

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from simple_crypt_count import (count_crypts_simple_method, load_final_outline,
                                load_outline_at_time, load_final_vertex_boundary,
                                load_boundary_from_vtu)

PARAM_NAME = 'ECM Stiffness'
PARAM_UNIT = 'Pa'
PARAM_JSON_PATH = 'forces.GhostNodeECM.ViscoelasticECM.ecmShearModulusPa'
TARGET_CIRCULARITY_TIME_H = 96.0


def run_circularity_at_day4(base_dir, model_type):
    """Load one VTU/VTP per run at day 4 and compute circularity = 4piA/P^2."""
    if model_type not in ('node2d', 'vertex2d'):
        return {}

    csvs = find_summary_csvs(base_dir)
    results = defaultdict(list)
    n = 0
    total = sum(1 for _, pp in csvs
                if identify_model_type('', load_params_json(pp) if pp else None) == model_type)

    for csv_path, params_path in csvs:
        params = load_params_json(params_path) if params_path else None
        if identify_model_type(csv_path, params) != model_type:
            continue

        stiffness = extract_param_from_json(params, PARAM_JSON_PATH) if params else None
        if stiffness is None:
            continue

        dt = params.get('simulation', {}).get('dt', 0.001) if params else 0.001
        data_dir = os.path.dirname(csv_path)
        n += 1

        if model_type == 'node2d':
            boundary, _ = load_outline_at_time(data_dir, TARGET_CIRCULARITY_TIME_H, dt)
        else:
            boundary, _ = load_final_vertex_boundary(
                data_dir, target_time=TARGET_CIRCULARITY_TIME_H, dt=dt)

        result = count_crypts_simple_method(boundary, boundary_is_ordered=True)
        results[stiffness].append(result.circularity)
        print(f"    {n}/{total}  G={stiffness:.0f}Pa  circ={result.circularity:.4f}")

    return dict(results)


def plot_circularity_summary(circ_results, plots_dir, model_type):
    """Boxplot and mean+/-SD of circularity at day 4 vs stiffness."""
    if not circ_results:
        return

    setup_style()
    stiffnesses = sorted(circ_results.keys())
    prefix = model_type
    day = TARGET_CIRCULARITY_TIME_H / 24.0
    colours = get_param_colours(stiffnesses)

    fig, axes = plt.subplots(1, 2, figsize=figsize(2, 1))

    ax = axes[0]
    data = [circ_results[s] for s in stiffnesses]
    bp = ax.boxplot(data, positions=range(len(stiffnesses)), widths=0.5,
                    patch_artist=True, medianprops=dict(color='black', linewidth=2))
    for i, patch in enumerate(bp['boxes']):
        patch.set_facecolor(colours[stiffnesses[i]])
        patch.set_alpha(0.4)
    for i, (s, d) in enumerate(zip(stiffnesses, data)):
        jitter = np.random.normal(0, 0.08, len(d))
        ax.scatter(np.full(len(d), i) + jitter, d,
                   color=colours[s], alpha=0.7, s=30, zorder=5)
    ax.set_xticks(range(len(stiffnesses)))
    ax.set_xticklabels([f'{s:.0f}' for s in stiffnesses], rotation=45)
    ax.set_xlabel(f'{PARAM_NAME} ({PARAM_UNIT})')
    ax.set_ylabel('Circularity (4\u03c0A/P\u00b2)')
    ax.set_ylim(0, 1.05)
    ax.grid(True, alpha=0.3)

    ax2 = axes[1]
    means = [np.mean(d) if d else np.nan for d in data]
    stds  = [np.std(d)  if d else 0      for d in data]
    ax2.errorbar(stiffnesses, means, yerr=stds, fmt='o-', capsize=5,
                 color='steelblue', linewidth=2, markersize=6)
    ax2.set_xlabel(f'{PARAM_NAME} ({PARAM_UNIT})')
    ax2.set_ylabel('Circularity (mean \u00b1 SD)')
    ax2.set_xscale('log')
    ax2.set_ylim(0, 1.05)
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    path = os.path.join(plots_dir, f'{prefix}_circularity_day{day:.0f}.png')
    fig.savefig(path, bbox_inches='tight')
    svg_path = os.path.splitext(path)[0] + '.svg'
    fig.savefig(svg_path, format='svg', bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {path}")

    csv_path = os.path.join(plots_dir, f'{prefix}_circularity_day{day:.0f}.csv')
    with open(csv_path, 'w') as f:
        f.write('stiffness_Pa,circularity\n')
        for s in stiffnesses:
            for c in circ_results[s]:
                f.write(f'{s},{c:.6f}\n')
    print(f"  Saved CSV: {csv_path}")


def run_crypt_counting(base_dir, model_type, plots_dir):
    """Run crypt counting on final-timestep output for 2D models."""
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

        data_dir = os.path.dirname(csv_path)

        if model_type == 'node2d':
            boundary, _ = load_final_outline(data_dir)
        else:
            vtu_files = sorted(glob.glob(os.path.join(data_dir, 'results_*.vtu')))
            if not vtu_files:
                continue
            boundary, _ = load_boundary_from_vtu(vtu_files[-1])

        result = count_crypts_simple_method(boundary, boundary_is_ordered=True)
        crypt_results[stiffness].append((result.num_crypts, result.circularity))

    return dict(crypt_results)


def plot_crypt_counts(crypt_results, plots_dir, model_type):
    """Plot crypt count and circularity vs stiffness."""
    if not crypt_results:
        return

    setup_style()
    stiffnesses = sorted(crypt_results.keys())

    fig, axes = plt.subplots(1, 2, figsize=figsize(2, 1))

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
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    ax2 = axes[1]
    circ_data = [[c for _, c in crypt_results[s]] for s in stiffnesses]
    means = [np.mean(d) if d else 0 for d in circ_data]
    stds = [np.std(d) if d else 0 for d in circ_data]
    ax2.errorbar(stiffnesses, means, yerr=stds, fmt='o-', capsize=5,
                 color='steelblue', markersize=8, linewidth=2)
    ax2.set_xlabel(f'{PARAM_NAME} ({PARAM_UNIT})')
    ax2.set_ylabel('Circularity')
    ax2.set_xscale('log')
    ax2.set_ylim(0, 1.05)

    plt.tight_layout()
    path = os.path.join(plots_dir, f'{model_type}_crypt_count_and_circularity.png')
    _save_fig(path)
    plt.close()
    print(f"  Saved: {path}")


def analyse_model(base_dir, model_type, plots_dir, full=False, crypt_count=True):
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

    generate_standard_plots(sweep_data, PARAM_NAME, PARAM_UNIT, model_type,
                            plots_dir, logx=True, crypt_count=crypt_count)

    if model_type in ('node2d', 'vertex2d'):
        print(f"  Computing circularity at day 4 (t={TARGET_CIRCULARITY_TIME_H:.0f}h) ...")
        circ_results = run_circularity_at_day4(base_dir, model_type)
        if circ_results:
            plot_circularity_summary(circ_results, plots_dir, model_type)

    if crypt_count:
        crypt_results = run_crypt_counting(base_dir, model_type, plots_dir)
        if crypt_results:
            plot_crypt_counts(crypt_results, plots_dir, model_type)
    else:
        print(f"  Skipping full crypt count (--no-crypt-count)")

    if full:
        print(f"  Running full VTU/VTP analysis for {model_type}...")
        analyse_vtu_data(base_dir, model_type, sweep_data, plots_dir)


def analyse_vtu_data(base_dir, model_type, sweep_data, plots_dir):
    """Extended analysis from VTU and VTP files."""
    setup_style()
    param_vals = sorted(sweep_data.keys())
    colours = get_param_colours(param_vals)

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

    fig, ax = plt.subplots(figsize=figsize())
    has_ci_data = False
    for pval in param_vals:
        if pval not in first_run_dirs:
            continue
        data_dir = first_run_dirs[pval]
        times, fracs = compute_contact_inhibition_timeseries(data_dir)
        if len(times) > 0:
            has_ci_data = True
            ax.plot(times, fracs, color=colours[pval], label=f'Stiffness={pval}')
    if has_ci_data:
        ax.set_xlabel('Timestep')
        ax.set_ylabel('Fraction Contact Inhibited')
        ax.legend(ncol=2)
        plt.tight_layout()
        path = os.path.join(plots_dir, f'{model_type}_contact_inhibition.png')
        _save_fig(path)
        print(f"  Saved: {path}")
    plt.close()

    mid_idx = len(param_vals) // 2
    mid_stiffness = param_vals[mid_idx] if param_vals else None
    if mid_stiffness and mid_stiffness in first_run_dirs:
        data_dir = first_run_dirs[mid_stiffness]
        ct_data = compute_cell_type_timeseries(data_dir)
        if 'times' in ct_data and len(ct_data['times']) > 0:
            fig, ax = plt.subplots(figsize=figsize())
            type_colours = {'Stem': 'green', 'Transit': 'blue',
                            'Paneth': 'red', 'Enterocyte': 'purple'}
            for name, colour in type_colours.items():
                if name in ct_data:
                    ax.plot(ct_data['times'], ct_data[name], color=colour, label=name)
            ax.set_xlabel('Timestep')
            ax.set_ylabel('Cell Type Fraction')
            ax.legend()
            ax.set_ylim(0, 1)
            plt.tight_layout()
            path = os.path.join(plots_dir, f'{model_type}_cell_type_ratios.png')
            _save_fig(path)
            plt.close()
            print(f"  Saved: {path}")

    fig, ax = plt.subplots(figsize=figsize())
    has_lumen = False
    for pval in param_vals:
        if pval not in first_run_dirs:
            continue
        data_dir = first_run_dirs[pval]
        times, forces = compute_lumen_force_timeseries(data_dir)
        if len(times) > 0:
            has_lumen = True
            ax.plot(times, forces, color=colours[pval], label=f'Stiffness={pval}')
    if has_lumen:
        ax.set_xlabel('Timestep')
        ax.set_ylabel('Mean Lumen Force Magnitude')
        ax.legend(ncol=2)
        plt.tight_layout()
        path = os.path.join(plots_dir, f'{model_type}_lumen_force.png')
        _save_fig(path)
        print(f"  Saved: {path}")
    plt.close()

    fig, axes = plt.subplots(1, 3, figsize=figsize(3, 1))
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
        axes[0].set_xlabel('Timestep'); axes[0].set_ylabel('Active Ghost Nodes')
        axes[0].legend(fontsize='small', ncol=2)
        axes[1].set_xlabel('Timestep'); axes[1].set_ylabel('Mean ECM Density')
        axes[1].legend(fontsize='small', ncol=2)
        axes[2].set_xlabel('Timestep'); axes[2].set_ylabel('Mean Rest Length Strain')
        axes[2].legend(fontsize='small', ncol=2)
        plt.tight_layout()
        path = os.path.join(plots_dir, f'{model_type}_ghost_node_ecm.png')
        _save_fig(path)
        print(f"  Saved: {path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Analyse ECM stiffness sweep and generate plots')
    parser.add_argument('--data-dir', required=True)
    parser.add_argument('--model', default='all',
                        choices=['node2d', 'vertex2d', 'node3d', 'vertex3d', 'all'])
    parser.add_argument('--output-dir', '-o', default=None)
    parser.add_argument('--full', action='store_true')
    parser.add_argument('--no-crypt-count', action='store_true')
    add_common_args(parser)

    args = parser.parse_args()
    apply_common_args(args)

    plots_dir = args.output_dir or os.path.join(args.data_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    models = (['node2d', 'vertex2d', 'node3d', 'vertex3d']
              if args.model == 'all' else [args.model])

    for model in models:
        analyse_model(args.data_dir, model, plots_dir,
                      full=args.full,
                      crypt_count=not args.no_crypt_count)

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
