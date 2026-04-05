#!/usr/bin/env python3
"""
analyse_and_plot.py — Lumen Pressure Sweep

Generates summary plots for lumen pressure sweep simulations.
Investigates how hydrostatic lumen pressure affects crypt budding,
organoid growth, and morphological instability.

All axes are converted to physical (SI) units:
  - Lumen pressure: dimensionless → Pa
  - Radial distances: CD → µm
  - Time: hours (already physical)

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
from analysis_utils import _save_fig

try:
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), 'ECMStiffnessSweep'))
    from simple_crypt_count import count_crypts_simple_method, load_final_outline
    HAS_CRYPT_COUNT = True
except ImportError:
    HAS_CRYPT_COUNT = False

PARAM_NAME = 'Lumen Pressure'
PARAM_UNIT = 'Pa'
PARAM_JSON_PATH = 'forces.LumenPressureForce.lumenPressure'

# ── Unit conversion constants ────────────────────────────────────────────
# Default physical parameters (can be overridden from params.json)
_DEFAULT_CD_UM = 8.0        # cell diameter in µm
_DEFAULT_F_REF_NN = 5.0     # reference force in nN
_DEFAULT_V_REF_UM_H = 10.0  # reference velocity in µm/h
_T0 = 3600.0                # seconds per hour


def _compute_conversions(cd_um=_DEFAULT_CD_UM, f_ref_nN=_DEFAULT_F_REF_NN,
                         v_ref_um_h=_DEFAULT_V_REF_UM_H):
    """Compute dimensionless → physical conversion factors.

    Returns (L0_um, pressure_to_Pa) where:
      L0_um:          multiply CD values by this to get µm
      pressure_to_Pa: multiply dimensionless pressure by this to get Pa
    """
    L0_m = cd_um * 1e-6
    f_ref_N = f_ref_nN * 1e-9
    v_ref_ms = v_ref_um_h * 1e-6 / _T0
    eta_phys = f_ref_N / v_ref_ms           # N·s/m
    F0 = eta_phys * L0_m / _T0              # reference force in N
    pressure_to_Pa = F0 / (L0_m ** 2)       # Pa per dimensionless pressure unit
    return cd_um, pressure_to_Pa


def _get_conversions_from_params(base_dir):
    """Try to read unit conversion parameters from a params.json in the data."""
    csvs = find_summary_csvs(base_dir)
    for _, params_path in csvs:
        if params_path:
            params = load_params_json(params_path)
            sim = params.get('simulation', params.get('units', {}))
            cd = sim.get('cellDiameterMicrometers', _DEFAULT_CD_UM)
            f_ref = sim.get('referenceForcenN', _DEFAULT_F_REF_NN)
            v_ref = sim.get('referenceVelocityUmPerH', _DEFAULT_V_REF_UM_H)
            return _compute_conversions(cd, f_ref, v_ref)
    return _compute_conversions()


def _convert_sweep_data(sweep_data, L0_um, pressure_to_Pa):
    """Convert sweep data in-place: pressure keys → Pa, distance columns → µm.

    Returns a new dict with converted pressure keys. Data arrays are scaled.
    """
    distance_cols = ('mean_r', 'var_r', 'r_range')
    # var_r is in CD² so needs L0² scaling
    converted = {}
    for pressure_dim, runs in sweep_data.items():
        pressure_pa = round(pressure_dim * pressure_to_Pa, 2)
        new_runs = []
        for run_num, data in runs:
            new_data = dict(data)
            for col in distance_cols:
                if col in new_data and isinstance(new_data[col], np.ndarray):
                    if col == 'var_r':
                        new_data[col] = new_data[col] * (L0_um ** 2)
                    else:
                        new_data[col] = new_data[col] * L0_um
            # Scale volume columns (area in 2D, CD² → µm²)
            for col in ('max_vol', 'min_vol'):
                if col in new_data and isinstance(new_data[col], np.ndarray):
                    new_data[col] = new_data[col] * (L0_um ** 2)
            new_runs.append((run_num, new_data))
        converted[pressure_pa] = new_runs
    return converted


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


def plot_pressure_specific(sweep_data, model_type, plots_dir, L0_um):
    """Plots specific to lumen pressure analysis (in physical units)."""
    if not HAS_MATPLOTLIB:
        return

    setup_style()
    param_vals = sorted(sweep_data.keys())
    colours = get_param_colours(param_vals)

    # Growth rate: derivative of mean_r (already in µm)
    fig, ax = plt.subplots(figsize=figsize())
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
                    label=f'P={pval} Pa')
    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('Radial Growth Rate (\u00b5m/h)')
    ax.set_title(f'{model_type}: Radial Growth Rate vs Lumen Pressure')
    ax.legend(ncol=2)
    ax.axhline(0, color='gray', linewidth=0.5)
    plt.tight_layout()
    path = os.path.join(plots_dir, f'{model_type}_growth_rate.png')
    _save_fig(path)
    plt.close()
    print(f"  Saved: {path}")

    # Ratio of r_range to mean_r (normalised budding amplitude, dimensionless)
    fig, ax = plt.subplots(figsize=figsize())
    for pval in param_vals:
        runs = sweep_data[pval]
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
                ax.plot(times, mean, color=colours[pval], label=f'P={pval} Pa')
                ax.fill_between(times, mean - std, mean + std,
                                color=colours[pval], alpha=0.15)
    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('Normalised Budding Amplitude (r_range / mean_r)')
    ax.set_title(f'{model_type}: Normalised Budding vs Lumen Pressure')
    ax.legend(ncol=2)
    plt.tight_layout()
    path = os.path.join(plots_dir, f'{model_type}_normalised_budding.png')
    _save_fig(path)
    plt.close()
    print(f"  Saved: {path}")


def analyse_model(base_dir, model_type, plots_dir, L0_um, pressure_to_Pa,
                  full=False):
    """Run analysis for one model type."""
    print(f"\n--- Analysing {model_type} (Lumen Pressure Sweep) ---")
    print(f"  Unit conversion: 1 CD = {L0_um} \u00b5m, "
          f"1 dim.pressure = {pressure_to_Pa:.2f} Pa")

    sweep = load_pressure_sweep(base_dir, model_filter=model_type)
    if model_type not in sweep:
        print(f"  No data found for {model_type}")
        return

    sweep_data = sweep[model_type]
    print(f"  Found {sum(len(v) for v in sweep_data.values())} runs across "
          f"{len(sweep_data)} pressure values")

    # Convert to physical units
    sweep_data = _convert_sweep_data(sweep_data, L0_um, pressure_to_Pa)

    # Override standard plot labels for physical units
    generate_standard_plots_physical(sweep_data, model_type, plots_dir, L0_um)
    plot_pressure_specific(sweep_data, model_type, plots_dir, L0_um)


def generate_standard_plots_physical(sweep_data, model_type, plots_dir, L0_um):
    """Generate standard plots with physical unit labels."""
    setup_style()
    prefix = model_type

    ts_plots = [
        ('num_cells', 'Number of Cells', 'Cell Proliferation Over Time'),
        ('mean_r', 'Mean Radial Distance (\u00b5m)',
         'Organoid Growth Over Time'),
        ('var_r', 'Radial Variance (\u00b5m\u00b2)',
         'Radial Variance Over Time'),
        ('r_range', 'Radial Range (\u00b5m)',
         'Budding Amplitude Over Time'),
        ('max_vol', 'Max Cell Area (\u00b5m\u00b2)',
         'Maximum Cell Area Over Time'),
        ('min_vol', 'Min Cell Area (\u00b5m\u00b2)',
         'Minimum Cell Area Over Time'),
    ]

    for col, ylabel, title in ts_plots:
        path = os.path.join(plots_dir, f'{prefix}_{col}_timeseries.png')
        plot_timeseries_by_param(sweep_data, col, PARAM_NAME, PARAM_UNIT,
                                ylabel, title, path, model_type)

    summary_plots = [
        ('num_cells', 'Final Cell Count',
         'Final Cell Count vs ' + PARAM_NAME),
        ('mean_r', 'Final Mean Radius (\u00b5m)',
         'Final Organoid Size vs ' + PARAM_NAME),
        ('var_r', 'Final Radial Variance (\u00b5m\u00b2)',
         'Final Radial Variance vs ' + PARAM_NAME),
        ('r_range', 'Final Radial Range (\u00b5m)',
         'Final Budding Amplitude vs ' + PARAM_NAME),
    ]

    for col, ylabel, title in summary_plots:
        path = os.path.join(plots_dir, f'{prefix}_{col}_boxplot.png')
        plot_final_value_boxplot(sweep_data, col, PARAM_NAME, PARAM_UNIT,
                                ylabel, title, path, model_type)
        path = os.path.join(plots_dir, f'{prefix}_{col}_vs_{PARAM_NAME}.png')
        plot_mean_vs_param(sweep_data, col, PARAM_NAME, PARAM_UNIT,
                           ylabel, title, path, model_type, logx=False)

    # Crypt counting (circularity is dimensionless, no conversion needed)
    if _crypt_count_mod and model_type in ('node2d', 'vertex2d'):
        print(f"  Computing crypt count timeseries for {model_type}...")
        enrich_sweep_with_crypt_counts(sweep_data, model_type)

        has_cc = any(
            'num_crypts' in d
            for runs in sweep_data.values()
            for _, d in runs
        )
        if has_cc:
            for col, ylabel, title in [
                ('circularity', 'Circularity', 'Circularity Over Time'),
            ]:
                path = os.path.join(plots_dir, f'{prefix}_{col}_timeseries.png')
                plot_timeseries_by_param(sweep_data, col, PARAM_NAME,
                                        PARAM_UNIT, ylabel, title, path,
                                        model_type, time_col='_cc_time')

            for col, ylabel, title in [
                ('circularity', 'Final Circularity',
                 'Final Circularity vs ' + PARAM_NAME),
            ]:
                path = os.path.join(plots_dir,
                                    f'{prefix}_{col}_boxplot.png')
                plot_final_value_boxplot(sweep_data, col, PARAM_NAME,
                                        PARAM_UNIT, ylabel, title, path,
                                        model_type)
                path = os.path.join(
                    plots_dir,
                    f'{prefix}_{col}_vs_{PARAM_NAME}.png')
                plot_mean_vs_param(sweep_data, col, PARAM_NAME, PARAM_UNIT,
                                   ylabel, title, path, model_type,
                                   logx=False)

            for col, ylabel, title in [
                ('num_crypts', 'Number of Crypts', 'Crypt Count Over Time'),
            ]:
                path = os.path.join(plots_dir,
                                    f'{prefix}_{col}_timeseries.png')
                plot_timeseries_by_param(sweep_data, col, PARAM_NAME,
                                        PARAM_UNIT, ylabel, title, path,
                                        model_type, time_col='_cc_time')

            for col, ylabel, title in [
                ('num_crypts', 'Final Crypt Count',
                 'Final Crypt Count vs ' + PARAM_NAME),
            ]:
                path = os.path.join(plots_dir,
                                    f'{prefix}_{col}_boxplot.png')
                plot_final_value_boxplot(sweep_data, col, PARAM_NAME,
                                        PARAM_UNIT, ylabel, title, path,
                                        model_type)
                path = os.path.join(
                    plots_dir,
                    f'{prefix}_{col}_vs_{PARAM_NAME}.png')
                plot_mean_vs_param(sweep_data, col, PARAM_NAME, PARAM_UNIT,
                                   ylabel, title, path, model_type,
                                   logx=False)

            cc_csv = os.path.join(plots_dir,
                                  f'{prefix}_crypt_count_timeseries.csv')
            save_crypt_count_csv(sweep_data, PARAM_NAME, cc_csv, model_type)

    print_summary_table(sweep_data, PARAM_NAME, model_type)
    csv_path = os.path.join(plots_dir, f'{prefix}_summary.csv')
    save_summary_csv(sweep_data, PARAM_NAME, csv_path, model_type)


def main():
    parser = argparse.ArgumentParser(
        description='Analyse lumen pressure sweep and generate plots')
    parser.add_argument('--data-dir', required=True)
    parser.add_argument('--model', default='all',
                        choices=['node2d', 'vertex2d', 'node3d', 'vertex3d',
                                 'all'])
    parser.add_argument('--output-dir', '-o', default=None)
    parser.add_argument('--full', action='store_true')
    add_common_args(parser)

    args = parser.parse_args()
    apply_common_args(args)
    plots_dir = args.output_dir or os.path.join(args.data_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    # Determine unit conversions (CLI overrides > params.json > defaults)
    cd = args.cell_diameter or _DEFAULT_CD_UM
    fr = args.ref_force or _DEFAULT_F_REF_NN
    vr = args.ref_velocity or _DEFAULT_V_REF_UM_H
    if args.cell_diameter or args.ref_force or args.ref_velocity:
        L0_um, pressure_to_Pa = _compute_conversions(cd, fr, vr)
    else:
        L0_um, pressure_to_Pa = _get_conversions_from_params(args.data_dir)

    models = ['node2d', 'vertex2d', 'node3d', 'vertex3d'] \
             if args.model == 'all' else [args.model]

    for model in models:
        analyse_model(args.data_dir, model, plots_dir, L0_um, pressure_to_Pa,
                      full=args.full)

    if len(models) > 1:
        print("\n--- Cross-model comparison ---")
        all_sweep = load_pressure_sweep(args.data_dir)
        if len(all_sweep) > 1:
            # Convert all models
            converted_sweep = {}
            for model, sweep_data in all_sweep.items():
                converted_sweep[model] = _convert_sweep_data(
                    sweep_data, L0_um, pressure_to_Pa)

            for col, ylabel, title in [
                ('num_cells', 'Final Cell Count',
                 'Cell Count vs Lumen Pressure'),
                ('mean_r', 'Final Mean Radius (\u00b5m)',
                 'Organoid Size vs Lumen Pressure'),
                ('r_range', 'Final Radial Range (\u00b5m)',
                 'Budding Amplitude vs Lumen Pressure'),
            ]:
                path = os.path.join(plots_dir,
                                    f'comparison_{col}_vs_pressure.png')
                plot_multi_model_comparison(converted_sweep, col, PARAM_NAME,
                                           PARAM_UNIT, ylabel, title, path,
                                           logx=False)

    print(f"\nAll plots saved to: {plots_dir}")


if __name__ == '__main__':
    main()
