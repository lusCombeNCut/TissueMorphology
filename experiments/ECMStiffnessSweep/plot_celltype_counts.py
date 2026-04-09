#!/usr/bin/env python3
"""
Cell type count plots for ECM stiffness sweep.

Reads from celltype_counts_all_runs.csv (produced by extract_celltype_counts.py).
Columns: model, stiffness_Pa, replicate, time_hours, stem, ta, differentiated, total

Produces:
  1. celltype_timeseries_vertex2d.png  — vertex2d Stem/TA/Diff timeseries, pooled mean ± SD
  2. celltype_timeseries_node2d.png    — node2d Stem/TA/Diff timeseries for runs with full data
  3. celltype_timeseries_combined.png  — side-by-side comparison of both models
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_utils import setup_style

CSV_PATH = '/home/orlando/Thesis/sim_results/2d_sim_04042026/plots_vertex2d/celltype_counts_all_runs.csv'
OUT_DIR  = '/home/orlando/Thesis/sim_results/2d_sim_04042026/plots_vertex2d'

# --- Font sizes ---
FS_LABEL     = 12
FS_TITLE     = 13
FS_SUPTITLE  = 14
FS_TICK      = 10
FS_LEGEND    = 11

COLOURS = {
    'stem':          '#2e86ab',
    'ta':            '#f4a261',
    'differentiated':'#e84855',
    'total':         '#4c4c4c',
}
LABELS = {
    'stem':          'Stem',
    'ta':            'TA',
    'differentiated':'Differentiated',
    'total':         'Total',
}

CELL_TYPES = ['stem', 'ta', 'differentiated']


# ─── helpers ─────────────────────────────────────────────────────────────────

def pool_timeseries(df, value_col):
    """
    Pool per-run timeseries into grand mean ± SD at each timepoint.
    Returns DataFrame with columns: time_hours, mean, std, n_runs.
    """
    grp = df.groupby('time_hours')[value_col]
    return pd.DataFrame({
        'time_hours': grp.mean().index,
        'mean':       grp.mean().values,
        'std':        grp.std(ddof=1).fillna(0).values,
        'n_runs':     grp.count().values,
    })


def draw_cell_type_panel(ax, df, title):
    """Draw Stem/TA/Diff mean ± SD timeseries for a given DataFrame onto ax."""
    for ct in CELL_TYPES:
        ts = pool_timeseries(df, ct)
        t  = ts['time_hours'].values
        m  = ts['mean'].values
        s  = ts['std'].values
        ax.plot(t, m, color=COLOURS[ct], linewidth=2, label=LABELS[ct])
        ax.fill_between(t, m - s, m + s, color=COLOURS[ct], alpha=0.18)

    ax.set_xlabel('Time (hours)', fontsize=FS_LABEL)
    ax.set_ylabel('Cell Count', fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    ax.set_title(title, fontsize=FS_TITLE)
    ax.legend(fontsize=FS_LEGEND)


# ─── individual plots ─────────────────────────────────────────────────────────

def plot_vertex2d(df, out_dir):
    setup_style()
    fig, ax = plt.subplots(figsize=(9, 5))
    draw_cell_type_panel(
        ax, df,
        'vertex2d: Cell Type Timeseries\n(all stiffnesses pooled, mean \u00b1 SD)'
    )
    fig.tight_layout()
    path = os.path.join(out_dir, 'celltype_timeseries_vertex2d.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {path}')


def plot_node2d(df_full, out_dir):
    """
    df_full: only rows where node2d has a full timeseries (>1 timepoint per run).
    """
    setup_style()
    stiffnesses = sorted(df_full['stiffness_Pa'].unique())
    stiff_str   = ', '.join(str(int(s)) for s in stiffnesses)
    n_runs      = df_full.groupby(['stiffness_Pa', 'replicate']).ngroups

    fig, ax = plt.subplots(figsize=(9, 5))
    draw_cell_type_panel(
        ax, df_full,
        f'node2d: Cell Type Timeseries\n(stiffnesses {stiff_str} Pa, {n_runs} runs pooled, mean \u00b1 SD)'
    )
    fig.tight_layout()
    path = os.path.join(out_dir, 'celltype_timeseries_node2d.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {path}')


def plot_combined(df_v2d, df_n2d_full, out_dir):
    """Side-by-side panel: vertex2d (all stiffnesses) | node2d (full-timeseries runs)."""
    setup_style()
    stiffnesses = sorted(df_n2d_full['stiffness_Pa'].unique())
    stiff_str   = ', '.join(str(int(s)) for s in stiffnesses)

    fig, axes = plt.subplots(1, 2, figsize=(16, 5), sharey=False)

    draw_cell_type_panel(axes[0], df_v2d,
                         'vertex2d\n(all stiffnesses pooled)')
    draw_cell_type_panel(axes[1], df_n2d_full,
                         f'node2d\n(stiffnesses {stiff_str} Pa)')

    # only show legend on left panel
    axes[1].get_legend().remove()

    fig.suptitle('Cell Type Counts vs Time (mean \u00b1 SD)', fontsize=FS_SUPTITLE)
    fig.tight_layout()
    path = os.path.join(out_dir, 'celltype_timeseries_combined.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {path}')


# ─── main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    os.makedirs(OUT_DIR, exist_ok=True)

    df = pd.read_csv(CSV_PATH)

    # vertex2d: all 50 runs, full timeseries
    df_v2d = df[df['model'] == 'vertex2d'].copy()

    # node2d: keep only runs that have more than one timestep
    n2d_ts_counts = (df[df['model'] == 'node2d']
                     .groupby(['stiffness_Pa', 'replicate'])['time_hours']
                     .count())
    full_runs = n2d_ts_counts[n2d_ts_counts > 1].reset_index()[['stiffness_Pa', 'replicate']]
    df_n2d_full = df[df['model'] == 'node2d'].merge(full_runs, on=['stiffness_Pa', 'replicate'])

    print(f'vertex2d runs: {df_v2d.groupby(["stiffness_Pa","replicate"]).ngroups}')
    print(f'node2d full-timeseries runs: {df_n2d_full.groupby(["stiffness_Pa","replicate"]).ngroups}')

    print('\nGenerating plots...')
    plot_vertex2d(df_v2d, OUT_DIR)
    plot_node2d(df_n2d_full, OUT_DIR)
    plot_combined(df_v2d, df_n2d_full, OUT_DIR)

    print('\nDone.')
