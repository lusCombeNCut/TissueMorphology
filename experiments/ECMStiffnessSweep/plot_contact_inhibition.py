#!/usr/bin/env python3
"""
Parse contactinhibition.dat files from node2d and vertex2d simulations,
compute proportion of contact-inhibited cells at each timepoint, write a CSV,
and produce comparison plots.

Note: node2d runs at G >= 700 Pa crashed (ECM instability) and have no timeseries data.
"""

import os
import sys
import re
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import csv

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from analysis_utils import setup_style, figsize

NODE2D_DIR   = '/home/orlando/Thesis/sim_results/2d_sim_04042026/16564369_node2d_2026-04-04_12-50-05'
VERTEX2D_DIR = '/home/orlando/Thesis/sim_results/2d_sim_04042026/16564370_vertex2d_2026-04-04_12-50-05'
OUT_DIR      = '/home/orlando/Thesis/sim_results/2d_sim_04042026/plots_vertex2d'

# ---------------------------------------------------------------------------
# Font sizes
# ---------------------------------------------------------------------------
LEGEND_FONTSIZE       = 8
LEGEND_TITLE_FONTSIZE = 8

RUN_RE = re.compile(r'^g(\d+)_r(\d+)$')


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_contactinhibition_dat(path):
    """
    Read a contactinhibition.dat file.

    Each line: <time_hours>  <0|1> <0|1> ... (one value per cell)

    Returns list of (time_h, proportion_inhibited) tuples.
    Skips lines with only t=0 data if the file has only 1 line (no timeseries).
    """
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            vals = line.split()
            if len(vals) < 2:
                continue
            time_h = float(vals[0])
            cell_vals = [float(v) for v in vals[1:]]
            if not cell_vals:
                continue
            prop = sum(v > 0.5 for v in cell_vals) / len(cell_vals)
            rows.append((time_h, prop, len(cell_vals)))
    return rows


def collect_all_runs(base_dir, model_name):
    """
    Walk all g{stiffness}_r{rep} directories, parse contactinhibition.dat.
    Returns list of dicts with keys: model, stiffness, replicate, time_h, prop_inhibited, n_cells
    """
    records = []
    run_dirs = sorted(os.listdir(base_dir))
    for run_name in run_dirs:
        m = RUN_RE.match(run_name)
        if not m:
            continue
        stiffness = int(m.group(1))
        replicate = int(m.group(2))

        ci_path = None
        crypt_dir = os.path.join(base_dir, run_name, 'CryptBudding')
        if not os.path.isdir(crypt_dir):
            continue
        for hash_dir in os.listdir(crypt_dir):
            candidate = os.path.join(crypt_dir, hash_dir,
                                     'results_from_time_0', 'contactinhibition.dat')
            if os.path.isfile(candidate):
                ci_path = candidate
                break

        if ci_path is None:
            continue

        rows = parse_contactinhibition_dat(ci_path)
        if not rows:
            continue
        # Skip runs that only have t=0 (crashed)
        if len(rows) == 1:
            continue

        for time_h, prop, n_cells in rows:
            records.append({
                'model':          model_name,
                'stiffness':      stiffness,
                'replicate':      replicate,
                'time_h':         time_h,
                'prop_inhibited': prop,
                'n_cells':        n_cells,
            })

    return records


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------

def write_csv(records, out_path):
    fields = ['model', 'stiffness', 'replicate', 'time_h', 'prop_inhibited', 'n_cells']
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(records)
    print(f'  Saved CSV: {out_path}')


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def final_props(records, model_name):
    """Return {stiffness: [proportion_inhibited, ...]} for the latest timepoint per run."""
    from collections import defaultdict

    # group by stiffness/replicate, keep final timepoint
    grouped = defaultdict(list)
    run_data = defaultdict(list)
    for r in records:
        if r['model'] != model_name:
            continue
        key = (r['stiffness'], r['replicate'])
        run_data[key].append((r['time_h'], r['prop_inhibited']))

    result = defaultdict(list)
    for (stiff, rep), ts in run_data.items():
        _, final_prop = max(ts, key=lambda x: x[0])
        result[stiff].append(final_prop)
    return dict(result)


def timeseries_mean(records, model_name):
    """
    Return {stiffness: {time_h: [prop, ...]}} for mean+/-SD timeseries.
    """
    from collections import defaultdict
    data = defaultdict(lambda: defaultdict(list))
    for r in records:
        if r['model'] != model_name:
            continue
        data[r['stiffness']][r['time_h']].append(r['prop_inhibited'])
    return {s: dict(t) for s, t in data.items()}


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

COLOURS = {'node2d': '#2e86ab', 'vertex2d': '#e84855'}


def plot_final_comparison(records, out_dir):
    """Scatter of final proportion inhibited vs stiffness, side-by-side panels."""
    setup_style()
    np.random.seed(42)

    fig, axes = plt.subplots(1, 2, figsize=figsize(2, 1), sharey=True)

    models = [('node2d', axes[0]), ('vertex2d', axes[1])]
    for model_name, ax in models:
        colour = COLOURS[model_name]
        fp = final_props(records, model_name)
        if not fp:
            ax.set_title(f'{model_name}\n(no data)', fontsize=13)
            continue

        stiffnesses = sorted(fp.keys())
        for s in stiffnesses:
            vals = np.array(fp[s])
            jitter = np.random.uniform(-0.03, 0.03, size=len(vals))
            xs = s * np.exp(jitter)
            ax.scatter(xs, vals * 100, color=colour, alpha=0.65, s=45,
                       linewidths=0, zorder=3)
            ax.scatter([s], [vals.mean() * 100], color=colour,
                       s=90, marker='D', edgecolors='black',
                       linewidths=0.8, zorder=4)

        ax.set_xscale('log')
        ax.set_xticks(stiffnesses)
        ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
        ax.tick_params(axis='x', rotation=45)
        ax.set_xlabel('ECM Stiffness (Pa)')
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f%%'))
        ax.grid(True, which='both', alpha=0.25)
        ax.set_title(model_name)

    axes[0].set_ylabel('Proportion Contact Inhibited')

    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='grey',
               alpha=0.7, markersize=7, label='Individual replicate'),
        Line2D([0], [0], marker='D', color='w', markerfacecolor='grey',
               markeredgecolor='black', markersize=8, label='Mean'),
    ]
    axes[1].legend(handles=handles, loc='best')

    # fig.suptitle('Final Proportion of Contact-Inhibited Cells vs ECM Stiffness', y=1.01)
    fig.tight_layout()
    path = os.path.join(out_dir, 'contact_inhibition_final_scatter.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {path}')


def plot_timeseries_comparison(records, out_dir):
    """Mean ± SD timeseries of proportion inhibited, one panel per model, all stiffnesses overlaid."""
    setup_style()

    # collect available stiffnesses per model
    node2d_ts  = timeseries_mean(records, 'node2d')
    vertex2d_ts = timeseries_mean(records, 'vertex2d')

    all_stiffnesses = sorted(
        set(list(node2d_ts.keys()) + list(vertex2d_ts.keys()))
    )

    # colour map across all stiffnesses
    cmap = plt.cm.viridis
    norm = plt.Normalize(vmin=min(all_stiffnesses), vmax=max(all_stiffnesses))

    fig, axes = plt.subplots(1, 2, figsize=figsize(2, 1), sharey=True)

    for (model_name, ts_data, ax) in [
        ('node2d',   node2d_ts,   axes[0]),
        ('vertex2d', vertex2d_ts, axes[1]),
    ]:
        stiffnesses = sorted(ts_data.keys())
        for s in stiffnesses:
            colour = cmap(norm(s))
            times = sorted(ts_data[s].keys())
            means = [np.mean(ts_data[s][t]) * 100 for t in times]
            sds   = [np.std(ts_data[s][t])  * 100 for t in times]
            ax.plot(times, means, color=colour, linewidth=1.5, label=f'{s} Pa')
            ax.fill_between(times,
                            [m - sd for m, sd in zip(means, sds)],
                            [m + sd for m, sd in zip(means, sds)],
                            color=colour, alpha=0.15)

        ax.set_xlabel('Time (hours)')
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f%%'))
        ax.grid(True, alpha=0.25)
        ax.set_title(model_name)
        if stiffnesses:
            ax.legend(fontsize=LEGEND_FONTSIZE, ncol=2, loc='lower right',
                      title='ECM stiffness', title_fontsize=LEGEND_TITLE_FONTSIZE)

    axes[0].set_ylabel('Proportion Contact Inhibited')
    # fig.suptitle('Contact Inhibition Over Time', y=1.01)
    fig.tight_layout()
    path = os.path.join(out_dir, 'contact_inhibition_timeseries.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {path}')


def plot_pooled_timeseries(records, out_dir):
    """
    Single panel: for each model, pool all stiffnesses together and plot
    mean ± SD proportion inhibited over time as a line + shaded area.
    """
    setup_style()
    from collections import defaultdict

    fig, ax = plt.subplots(figsize=figsize(1.3, 1))

    for model_name in ('node2d', 'vertex2d'):
        colour = COLOURS[model_name]
        # aggregate all runs at each timepoint regardless of stiffness
        time_data = defaultdict(list)
        for r in records:
            if r['model'] != model_name:
                continue
            time_data[r['time_h']].append(r['prop_inhibited'])

        if not time_data:
            continue

        times = sorted(time_data.keys())
        means = np.array([np.mean(time_data[t]) * 100 for t in times])
        sds   = np.array([np.std(time_data[t])  * 100 for t in times])

        ax.plot(times, means, color=colour, linewidth=2.0, label=model_name)
        ax.fill_between(times, means - sds, means + sds,
                        color=colour, alpha=0.2)

    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('Proportion Contact Inhibited')
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f%%'))
    # ax.set_title('Contact Inhibition Over Time (all stiffnesses pooled)')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.25)
    ax.set_ylim(0, 100)
    fig.tight_layout()
    path = os.path.join(out_dir, 'contact_inhibition_pooled_timeseries.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {path}')


def plot_combined_final(records, out_dir):
    """Single panel comparing node2d and vertex2d final proportions."""
    setup_style()
    np.random.seed(42)

    fig, ax = plt.subplots(figsize=figsize(1.3, 1))

    for model_name in ('node2d', 'vertex2d'):
        colour = COLOURS[model_name]
        fp = final_props(records, model_name)
        if not fp:
            continue
        stiffnesses = sorted(fp.keys())
        first = True
        for s in stiffnesses:
            vals = np.array(fp[s])
            jitter = np.random.uniform(-0.03, 0.03, size=len(vals))
            xs = s * np.exp(jitter)
            ax.scatter(xs, vals * 100, color=colour, alpha=0.55, s=35,
                       linewidths=0, zorder=3,
                       label=model_name if first else '_nolegend_')
            first = False

    ax.set_xscale('log')
    all_stiffs = sorted(set(
        [r['stiffness'] for r in records]
    ))
    ax.set_xticks(all_stiffs)
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.tick_params(axis='x', rotation=45)
    ax.set_xlabel('ECM Stiffness (Pa)')
    ax.set_ylabel('Final Proportion Contact Inhibited')
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f%%'))
    ax.grid(True, which='both', alpha=0.25)
    ax.legend(loc='lower right')
    # ax.set_title('Final Contact Inhibition vs ECM Stiffness')
    fig.tight_layout()
    path = os.path.join(out_dir, 'contact_inhibition_final_combined.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {path}')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    os.makedirs(OUT_DIR, exist_ok=True)

    print('Parsing node2d contact inhibition data...')
    node2d_records = collect_all_runs(NODE2D_DIR, 'node2d')
    print(f'  {len(node2d_records)} timepoint records from node2d')

    print('Parsing vertex2d contact inhibition data...')
    vertex2d_records = collect_all_runs(VERTEX2D_DIR, 'vertex2d')
    print(f'  {len(vertex2d_records)} timepoint records from vertex2d')

    all_records = node2d_records + vertex2d_records

    csv_path = os.path.join(OUT_DIR, 'contact_inhibition_timeseries.csv')
    write_csv(all_records, csv_path)

    print('\nPlotting...')
    plot_final_comparison(all_records, OUT_DIR)
    plot_timeseries_comparison(all_records, OUT_DIR)
    plot_combined_final(all_records, OUT_DIR)
    plot_pooled_timeseries(all_records, OUT_DIR)

    # Summary stats
    print('\nSummary — final proportion contact inhibited:')
    for model_name in ('node2d', 'vertex2d'):
        fp = final_props(all_records, model_name)
        print(f'  {model_name}:')
        for s in sorted(fp.keys()):
            vals = np.array(fp[s])
            print(f'    G={s:5d} Pa  n={len(vals)}  mean={vals.mean()*100:.1f}%  sd={vals.std()*100:.1f}%')
