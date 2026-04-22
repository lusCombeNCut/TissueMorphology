#!/usr/bin/env python3
"""
lumen_correlations.py — Correlation statistics for the Lumen Pressure Sweep

For each replicate's final-timepoint value computes Pearson r, Spearman ρ,
Kendall τ, and a linear regression (slope, intercept, R²) for:
  • final mean radial distance (µm) vs lumen pressure (Pa)
  • final cell count          vs lumen pressure (Pa)

Results are printed to stdout and saved to lumen_correlation_stats.csv in the
requested output directory.

Usage:
  python3 lumen_correlations.py \
      --data-dir /path/to/16568222_vertex2d_2026-04-05_00-23-06 \
      --output-dir /path/to/lumen_plots
"""

import os
import sys
import argparse
import csv

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import importlib.util

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))

from analysis_utils import find_summary_csvs, load_params_json, load_crypt_summary

# Load the LumenPressureSweep analyse_and_plot.py explicitly by path
_spec = importlib.util.spec_from_file_location(
    'lumen_analyse_and_plot',
    os.path.join(_HERE, 'analyse_and_plot.py'))
_aap = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_aap)

load_pressure_sweep = _aap.load_pressure_sweep
_convert_sweep_data = _aap._convert_sweep_data
_get_conversions_from_params = _aap._get_conversions_from_params
_compute_conversions = _aap._compute_conversions
_DEFAULT_CD_UM = _aap._DEFAULT_CD_UM
_DEFAULT_F_REF_NN = _aap._DEFAULT_F_REF_NN
_DEFAULT_V_REF_UM_H = _aap._DEFAULT_V_REF_UM_H


def collect_final_values(sweep_data):
    """Return arrays (pressures, mean_r_finals, cell_count_finals) for every replicate."""
    pressures, mean_r_vals, cell_vals = [], [], []
    for pressure_pa, runs in sorted(sweep_data.items()):
        for _run_num, data in runs:
            if 'mean_r' in data and len(data['mean_r']) > 0:
                mean_r_vals.append(data['mean_r'][-1])
                cell_vals.append(data['num_cells'][-1])
                pressures.append(pressure_pa)
    return (np.array(pressures, dtype=float),
            np.array(mean_r_vals, dtype=float),
            np.array(cell_vals, dtype=float))


def correlation_report(x, y, x_label, y_label):
    """Compute and print Pearson r, Spearman ρ, Kendall τ, and linear regression."""
    n = len(x)

    pearson_r, pearson_p = stats.pearsonr(x, y)
    spearman_r, spearman_p = stats.spearmanr(x, y)
    kendall_tau, kendall_p = stats.kendalltau(x, y)
    slope, intercept, lin_r, lin_p, lin_se = stats.linregress(x, y)
    r_squared = lin_r ** 2

    print(f"\n{'='*65}")
    print(f"  {y_label} vs {x_label}  (n = {n} replicates)")
    print(f"{'='*65}")
    print(f"  Pearson r          = {pearson_r:+.4f}   p = {pearson_p:.4e}")
    print(f"  Spearman ρ         = {spearman_r:+.4f}   p = {spearman_p:.4e}")
    print(f"  Kendall τ          = {kendall_tau:+.4f}   p = {kendall_p:.4e}")
    print(f"  Linear regression:")
    print(f"    slope            = {slope:.4e}")
    print(f"    intercept        = {intercept:.4f}")
    print(f"    R²               = {r_squared:.4f}")
    print(f"    regression p     = {lin_p:.4e}")
    print(f"    std error        = {lin_se:.4e}")
    print(f"{'='*65}")

    return {
        'y_variable': y_label,
        'x_variable': x_label,
        'n': n,
        'pearson_r': pearson_r,
        'pearson_p': pearson_p,
        'spearman_rho': spearman_r,
        'spearman_p': spearman_p,
        'kendall_tau': kendall_tau,
        'kendall_p': kendall_p,
        'slope': slope,
        'intercept': intercept,
        'r_squared': r_squared,
        'regression_p': lin_p,
        'std_error': lin_se,
    }


def boxplot_with_fit(sweep_data, column, x_label, y_label, title,
                     output_path, stats_row):
    """Redraw a final-value boxplot with a dashed best-fit line and Pearson annotation."""
    plt.rcParams.update({
        'font.family': 'serif',
        'mathtext.fontset': 'dejavuserif',
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'figure.constrained_layout.use': True,
        'font.size': 22,
        'axes.titlesize': 26,
        'axes.labelsize': 24,
        'legend.fontsize': 20,
        'xtick.labelsize': 20,
        'ytick.labelsize': 20,
        'axes.grid': True,
        'grid.alpha': 0.3,
    })

    param_vals = sorted(sweep_data.keys())
    positions = list(range(len(param_vals)))

    # Collect final values per pressure
    final_values = []
    for pval in param_vals:
        vals = [data[column][-1]
                for _, data in sweep_data[pval]
                if column in data and len(data[column]) > 0]
        final_values.append(vals)

    fig, ax = plt.subplots(figsize=(9, 7))
    bp = ax.boxplot(final_values, positions=positions, widths=0.6,
                    patch_artist=True,
                    medianprops=dict(color='black', linewidth=2))
    for patch in bp['boxes']:
        patch.set_facecolor('steelblue')
        patch.set_alpha(0.4)

    # Refit regression on integer positions so the line is straight on the
    # categorical x-axis (pressure spacing is uneven, so using Pa values
    # directly would produce a curved appearance).
    all_pos, all_y = [], []
    for i, pval in enumerate(param_vals):
        for _, data in sweep_data[pval]:
            if column in data and len(data[column]) > 0:
                all_pos.append(i)
                all_y.append(data[column][-1])
    from scipy.stats import linregress as _lr
    fit = _lr(all_pos, all_y)
    x_fit = np.linspace(positions[0], positions[-1], 300)
    y_fit = fit.slope * x_fit + fit.intercept
    ax.plot(x_fit, y_fit, color='firebrick', linewidth=2,
            linestyle='--', label='Linear fit', zorder=3)

    # Exponential fit (y = a * exp(b * Pa)) on Pa values — same x-space as
    # the linear fit so R² values are directly comparable.
    # For plotting, evaluate the fit at a dense set of Pa values and map to
    # plot positions using np.interp on the exact knots (exact at box centres,
    # smoothly interpolated between them).
    exp_row = None
    if column == 'mean_r':
        all_pa = np.array([pval for pval in param_vals
                           for _, data in sweep_data[pval]
                           if column in data and len(data[column]) > 0])
        all_y_arr = np.array(all_y)
        def _exp_func(xv, a, b):
            return a * np.exp(b * xv)
        try:
            popt, _ = curve_fit(_exp_func, all_pa, all_y_arr,
                                p0=[all_y_arr.min(), 1e-3], maxfev=10000)
            a_exp, b_exp = popt
            # Dense Pa values from min to max, mapped to plot positions via
            # linear interp between knots — exact at each box, smooth between.
            pa_knots = np.array(param_vals, dtype=float)
            pa_dense = np.linspace(pa_knots[0], pa_knots[-1], 1000)
            pos_dense = np.interp(pa_dense, pa_knots, positions)
            y_exp_dense = _exp_func(pa_dense, a_exp, b_exp)
            ax.plot(pos_dense, y_exp_dense, color='darkgreen', linewidth=2,
                    linestyle=':', label='Exponential fit', zorder=3)
            # R² and F-test p-value on Pa values (same basis as linear fit)
            y_pred = _exp_func(all_pa, a_exp, b_exp)
            ss_res = np.sum((all_y_arr - y_pred) ** 2)
            ss_tot = np.sum((all_y_arr - all_y_arr.mean()) ** 2)
            r_sq = 1.0 - ss_res / ss_tot
            n_obs, n_params = len(all_y_arr), 2  # a, b
            f_stat = (r_sq / n_params) / ((1 - r_sq) / (n_obs - n_params - 1))
            exp_p = stats.f.sf(f_stat, n_params, n_obs - n_params - 1)
            exp_row = {'a': a_exp, 'b': b_exp, 'r_sq': r_sq, 'p': exp_p}
        except RuntimeError:
            pass

    # Top-left: R² and p (linear)
    r_sq_lin = stats_row['r_squared']
    p = stats_row['pearson_p']
    p_str = f'{p:.2e}'.replace('e-0', 'e-').replace('e+0', 'e+')
    lin_lines = f'Linear:  $R^2 = {r_sq_lin:.3f}$, $p = {p_str}$'
    if exp_row:
        exp_p_str = f'{exp_row["p"]:.2e}'.replace('e-0', 'e-').replace('e+0', 'e+')
        exp_lines = f'\nExponential:  $R^2 = {exp_row["r_sq"]:.3f}$, $p = {exp_p_str}$'
    else:
        exp_lines = ''
    ax.text(0.05, 0.95, lin_lines + exp_lines,
            transform=ax.transAxes,
            verticalalignment='top',
            fontsize=18,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='gray', alpha=0.8))

    # Bottom-right: equations
    s = stats_row['slope']
    b = stats_row['intercept']
    sign = '+' if b >= 0 else '-'
    eq_lin = f'$y = {s:.3f}x {sign} {abs(b):.2f}$'
    if exp_row:
        eq_exp = f'\n$y = {exp_row["a"]:.2f} \\, e^{{{exp_row["b"]:.4f} x}}$'
    else:
        eq_exp = ''
    ax.text(0.97, 0.05, eq_lin + eq_exp,
            transform=ax.transAxes,
            verticalalignment='bottom',
            horizontalalignment='right',
            fontsize=18,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='gray', alpha=0.8))

    ax.set_xticks(positions)
    ax.set_xticklabels([str(v) for v in param_vals], rotation=45)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    svg_path = os.path.splitext(output_path)[0] + '.svg'
    plt.savefig(svg_path, format='svg')
    plt.close()
    print(f"  Saved: {svg_path}")


def scatter_with_fit(x, y, x_label, y_label, title, output_path, stats_row,
                     column=''):
    """Scatter plot of individual replicates with linear best-fit line and Pearson annotation."""
    plt.rcParams.update({
        'font.family': 'serif',
        'mathtext.fontset': 'dejavuserif',
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'figure.constrained_layout.use': True,
        'font.size': 22,
        'axes.titlesize': 26,
        'axes.labelsize': 24,
        'legend.fontsize': 20,
        'xtick.labelsize': 20,
        'ytick.labelsize': 20,
        'axes.grid': True,
        'grid.alpha': 0.3,
    })

    fig, ax = plt.subplots(figsize=(7, 7))

    ax.scatter(x, y, color='steelblue', alpha=0.7, zorder=3, s=40)

    # Linear best-fit line
    x_fit = np.linspace(x.min(), x.max(), 300)
    y_fit = stats_row['slope'] * x_fit + stats_row['intercept']
    ax.plot(x_fit, y_fit, color='firebrick', linewidth=2,
            linestyle='--', label='Linear fit')

    # Exponential fit for mean_r
    exp_row = None
    if column == 'mean_r':
        def _exp_func(xv, a, b):
            return a * np.exp(b * xv)
        try:
            popt, _ = curve_fit(_exp_func, x, y,
                                p0=[y.min(), 1e-3], maxfev=10000)
            a_exp, b_exp = popt
            y_exp = _exp_func(x_fit, a_exp, b_exp)
            ax.plot(x_fit, y_exp, color='darkgreen', linewidth=2,
                    linestyle=':', label='Exponential fit')
            # R² and F-test p-value for the nonlinear fit
            y_pred = _exp_func(x, a_exp, b_exp)
            ss_res = np.sum((y - y_pred) ** 2)
            ss_tot = np.sum((y - y.mean()) ** 2)
            r_sq = 1.0 - ss_res / ss_tot
            n_obs, n_params = len(y), 2  # a, b
            f_stat = (r_sq / n_params) / ((1 - r_sq) / (n_obs - n_params - 1))
            exp_p = stats.f.sf(f_stat, n_params, n_obs - n_params - 1)
            exp_row = {'a': a_exp, 'b': b_exp, 'r_sq': r_sq, 'p': exp_p}
        except RuntimeError:
            pass

    # Top-left: R² and p values
    r_sq_lin = stats_row['r_squared']
    p = stats_row['pearson_p']
    p_str = f'{p:.2e}'.replace('e-0', 'e-').replace('e+0', 'e+')
    lin_lines = f'Linear:  $R^2 = {r_sq_lin:.3f}$, $p = {p_str}$'
    if exp_row:
        exp_p_str = f'{exp_row["p"]:.2e}'.replace('e-0', 'e-').replace('e+0', 'e+')
        lin_lines += f'\nExponential:  $R^2 = {exp_row["r_sq"]:.3f}$, $p = {exp_p_str}$'
    ax.text(0.05, 0.95, lin_lines,
            transform=ax.transAxes,
            verticalalignment='top',
            fontsize=18,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='gray', alpha=0.8))

    # Bottom-right: equations
    s = stats_row['slope']
    b = stats_row['intercept']
    sign = '+' if b >= 0 else '-'
    eq_str = f'$y = {s:.3f}x {sign} {abs(b):.2f}$'
    if exp_row:
        eq_str += f'\n$y = {exp_row["a"]:.2f} \\, e^{{{exp_row["b"]:.4f} x}}$'
    ax.text(0.97, 0.05, eq_str,
            transform=ax.transAxes,
            verticalalignment='bottom',
            horizontalalignment='right',
            fontsize=18,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='gray', alpha=0.8))

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    svg_path = os.path.splitext(output_path)[0] + '.svg'
    plt.savefig(svg_path, format='svg')
    plt.close()
    print(f"  Saved: {svg_path}")


def save_stats_csv(rows, output_path):
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n  Saved correlation stats: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Correlation statistics for lumen pressure sweep')
    parser.add_argument('--data-dir', required=True)
    parser.add_argument('--model', default='vertex2d',
                        choices=['node2d', 'vertex2d', 'node3d', 'vertex3d'])
    parser.add_argument('--output-dir', '-o', default=None)
    parser.add_argument('--cell-diameter', type=float, default=None)
    parser.add_argument('--ref-force', type=float, default=None)
    parser.add_argument('--ref-velocity', type=float, default=None)
    args = parser.parse_args()

    plots_dir = args.output_dir or os.path.join(args.data_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    cd = args.cell_diameter or _DEFAULT_CD_UM
    fr = args.ref_force or _DEFAULT_F_REF_NN
    vr = args.ref_velocity or _DEFAULT_V_REF_UM_H
    if args.cell_diameter or args.ref_force or args.ref_velocity:
        L0_um, pressure_to_Pa = _compute_conversions(cd, fr, vr)
    else:
        L0_um, pressure_to_Pa = _get_conversions_from_params(args.data_dir)

    print(f"\nUnit conversion: 1 CD = {L0_um} µm, 1 dim. = {pressure_to_Pa:.2f} Pa")

    sweep_raw = load_pressure_sweep(args.data_dir, model_filter=args.model)
    if args.model not in sweep_raw:
        print(f"No data found for {args.model}")
        sys.exit(1)

    sweep_data = _convert_sweep_data(sweep_raw[args.model], L0_um, pressure_to_Pa)

    pressures, mean_r_vals, cell_vals = collect_final_values(sweep_data)

    rows = []
    rows.append(correlation_report(
        pressures, mean_r_vals,
        'Lumen Pressure (Pa)', 'Final Mean Radius (µm)'))
    rows.append(correlation_report(
        pressures, cell_vals,
        'Lumen Pressure (Pa)', 'Final Cell Count'))

    csv_path = os.path.join(plots_dir, 'lumen_correlation_stats.csv')
    save_stats_csv(rows, csv_path)

    scatter_with_fit(
        pressures, mean_r_vals,
        'Lumen Pressure (Pa)', 'Final Mean Radial Distance (µm)',
        'Final Mean Radial Distance vs Lumen Pressure',
        os.path.join(plots_dir, 'vertex2d_mean_r_vs_pressure_scatter.svg'),
        rows[0], column='mean_r')

    scatter_with_fit(
        pressures, cell_vals,
        'Lumen Pressure (Pa)', 'Final Cell Count',
        'Final Cell Count vs Lumen Pressure',
        os.path.join(plots_dir, 'vertex2d_cell_count_vs_pressure_scatter.svg'),
        rows[1], column='num_cells')

    boxplot_with_fit(
        sweep_data, 'mean_r',
        'Lumen Pressure (Pa)', 'Final Mean Radial Distance (µm)',
        'Final Mean Radial Distance vs Lumen Pressure',
        os.path.join(plots_dir, 'vertex2d_mean_r_boxplot.svg'),
        rows[0])

    boxplot_with_fit(
        sweep_data, 'num_cells',
        'Lumen Pressure (Pa)', 'Final Cell Count',
        'Final Cell Count vs Lumen Pressure',
        os.path.join(plots_dir, 'vertex2d_num_cells_boxplot.svg'),
        rows[1])


if __name__ == '__main__':
    main()
