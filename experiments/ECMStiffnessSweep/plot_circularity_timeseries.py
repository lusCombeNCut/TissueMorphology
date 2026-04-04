#!/usr/bin/env python3
"""Compute circularity timeseries for all runs and plot mean +/- SD over time."""

import os
import sys
import glob
import json
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from simple_crypt_count import (
    load_outline_from_vtp, load_boundary_from_vtu, compute_circularity
)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_utils import (
    find_summary_csvs, load_params_json, identify_model_type,
    extract_param_from_json, get_param_colours
)

PARAM_JSON_PATH = 'forces.GhostNodeECM.ViscoelasticECM.ecmShearModulusPa'

plt.rcParams.update({
    'font.family': 'serif',
    'mathtext.fontset': 'dejavuserif',
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'font.size': 11,
    'axes.labelsize': 13,
    'legend.fontsize': 9,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'lines.linewidth': 1.8,
    'axes.grid': True,
    'grid.alpha': 0.3,
})


def compute_circularity_timeseries(data_dir, model_type, dt):
    """Compute circularity at every saved timestep for one run.

    Returns arrays (times_hours, circularities).
    """
    if model_type == 'node2d':
        files = sorted(glob.glob(os.path.join(data_dir, 'outline_*.vtp')))
        def get_step(f):
            return int(os.path.basename(f).replace('outline_', '').replace('.vtp', ''))
    else:
        files = sorted([f for f in glob.glob(os.path.join(data_dir, 'results_*.vtu'))
                        if 'ecm_grid' not in f])
        def get_step(f):
            return int(os.path.basename(f).replace('results_', '').replace('.vtu', ''))

    files.sort(key=get_step)
    if not files:
        return np.array([]), np.array([])

    times = []
    circs = []
    for f in files:
        step = get_step(f)
        t_hours = step * dt

        if model_type == 'node2d':
            boundary, _ = load_outline_from_vtp(f)
        else:
            boundary, _ = load_boundary_from_vtu(f)

        c = compute_circularity(boundary)
        times.append(t_hours)
        circs.append(c)

    return np.array(times), np.array(circs)


def gather_all_runs(base_dir, model_type):
    """Gather circularity timeseries for all runs, grouped by stiffness.

    Returns dict[stiffness] -> list of (times, circs) arrays.
    """
    csvs = find_summary_csvs(base_dir)
    results = {}
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

        dt = params.get('simulation', {}).get('dt', 0.001)
        data_dir = os.path.dirname(csv_path)
        n += 1
        print(f"  [{n}/{total}] G={stiffness:.0f} Pa — {os.path.basename(os.path.dirname(data_dir))}")

        times, circs = compute_circularity_timeseries(data_dir, model_type, dt)
        if len(times) > 0:
            results.setdefault(stiffness, []).append((times, circs))

    return results


def plot_circularity_over_time(results, model_type, output_path):
    """Plot mean circularity +/- SD over time, one line per stiffness."""
    stiffnesses = sorted(results.keys())
    colours = get_param_colours(stiffnesses)

    fig, ax = plt.subplots(figsize=(8, 5))

    for s in stiffnesses:
        runs = results[s]
        ref_times = max(runs, key=lambda r: len(r[0]))[0]

        interp_circs = np.full((len(runs), len(ref_times)), np.nan)
        for i, (t, c) in enumerate(runs):
            if len(t) > 1:
                interp_circs[i, :] = np.interp(ref_times, t, c)
            elif len(t) == 1:
                interp_circs[i, :] = c[0]

        mean = np.nanmean(interp_circs, axis=0)
        std = np.nanstd(interp_circs, axis=0)

        ax.plot(ref_times, mean, color=colours[s], label=f'{s:.0f} Pa')
        ax.fill_between(ref_times, mean - std, mean + std,
                        color=colours[s], alpha=0.15)

    ax.set_xlabel('Time (hours)')
    ax.set_ylabel(r'Circularity ($4\pi A / P^2$)')
    ax.set_ylim(0, 1.05)
    ax.legend(title='ECM Shear\nModulus', ncol=2, loc='best')

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {output_path}")


def main():
    base = '/home/orlando/Thesis/Chaste/projects/TissueMorphology/experiments/ECMStiffnessSweep/2d_sim_04042026'
    out_dir = os.path.join(base, 'plots_circularity_timeseries')
    os.makedirs(out_dir, exist_ok=True)

    node2d_dir = os.path.join(base, '16564369_node2d_2026-04-04_12-50-05')
    vertex2d_dir = os.path.join(base, '16564370_vertex2d_2026-04-04_12-50-05')

    for model_type, data_dir in [('node2d', node2d_dir), ('vertex2d', vertex2d_dir)]:
        print(f"\n=== {model_type} ===")
        results = gather_all_runs(data_dir, model_type)
        if not results:
            print(f"  No data found for {model_type}")
            continue
        out_path = os.path.join(out_dir, f'{model_type}_circularity_timeseries.png')
        plot_circularity_over_time(results, model_type, out_path)

        csv_path = os.path.join(out_dir, f'{model_type}_circularity_timeseries.csv')
        with open(csv_path, 'w') as f:
            f.write('model,stiffness_Pa,time_hours,mean_circularity,std_circularity,n_runs\n')
            for s in sorted(results.keys()):
                runs = results[s]
                ref_times = max(runs, key=lambda r: len(r[0]))[0]
                interp = np.full((len(runs), len(ref_times)), np.nan)
                for i, (t, c) in enumerate(runs):
                    if len(t) > 1:
                        interp[i, :] = np.interp(ref_times, t, c)
                for j, t in enumerate(ref_times):
                    vals = interp[:, j]
                    valid = vals[~np.isnan(vals)]
                    f.write(f'{model_type},{s},{t:.2f},{np.mean(valid):.6f},{np.std(valid):.6f},{len(valid)}\n')
        print(f"Saved CSV: {csv_path}")


if __name__ == '__main__':
    main()
