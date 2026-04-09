#!/usr/bin/env python3
"""
Plot mean ± SD circularity timeseries pooled across ALL ECM stiffnesses,
for node2d and vertex2d on the same axes.

Grand statistics are computed from per-stiffness (mean, std, n) using the
law of total variance:
    grand_mean = mean of group means  (equal group sizes)
    grand_var  = mean(var_i + mean_i^2) - grand_mean^2
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_utils import setup_style

# --- Paths ---
NODE2D_CSV   = '/home/orlando/Thesis/sim_results/2d_sim_04042026/plots_circularity_timeseries/node2d_circularity_timeseries.csv'
VERTEX2D_CSV = '/home/orlando/Thesis/sim_results/2d_sim_04042026/plots_circularity_timeseries/vertex2d_circularity_timeseries.csv'
OUT_DIR      = '/home/orlando/Thesis/sim_results/2d_sim_04042026/plots_circularity_timeseries'

# --- Font sizes ---
FS_LABEL   = 14
FS_TITLE   = 13
FS_TICK    = 12
FS_LEGEND  = 12


def pool_timeseries(csv_path):
    """
    Return (times, grand_mean, grand_std) arrays pooled across all stiffnesses.
    Uses the law of total variance on the per-stiffness (mean, std, n) rows.
    """
    df = pd.read_csv(csv_path)
    times = sorted(df['time_hours'].unique())
    grand_means, grand_stds = [], []

    for t in times:
        rows = df[df['time_hours'] == t]
        means = rows['mean_circularity'].values
        stds  = rows['std_circularity'].values
        ns    = rows['n_runs'].values

        grand_mean = np.average(means, weights=ns)
        # law of total variance: E[Var] + Var[E]
        grand_var  = np.average(stds**2 + means**2, weights=ns) - grand_mean**2
        grand_std  = np.sqrt(max(grand_var, 0.0))

        grand_means.append(grand_mean)
        grand_stds.append(grand_std)

    return np.array(times), np.array(grand_means), np.array(grand_stds)


def plot_pooled(out_dir):
    setup_style()
    os.makedirs(out_dir, exist_ok=True)

    datasets = [
        ('node2d',   NODE2D_CSV,   '#2e86ab'),
        ('vertex2d', VERTEX2D_CSV, '#e84855'),
    ]

    fig, ax = plt.subplots(figsize=(5, 4))

    for label, csv_path, colour in datasets:
        times, means, stds = pool_timeseries(csv_path)
        ax.plot(times, means, color=colour, linewidth=2, label=label)
        ax.fill_between(times, means - stds, means + stds,
                        color=colour, alpha=0.20)

    ax.set_xlabel('Time (hours)', fontsize=FS_LABEL)
    ax.set_ylabel('Circularity', fontsize=FS_LABEL)
    ax.set_ylim(0, 1.05)
    ax.tick_params(labelsize=FS_TICK)
    ax.legend(fontsize=FS_LEGEND)
    # ax.set_title('Circularity Timeseries (pooled across all stiffnesses)', fontsize=FS_TITLE)

    fig.tight_layout()
    out_path = os.path.join(out_dir, 'circularity_pooled_timeseries.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {out_path}')


if __name__ == '__main__':
    plot_pooled(OUT_DIR)
