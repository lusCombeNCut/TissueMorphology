#!/usr/bin/env python3
"""
Scatter plot of final crypt count vs ECM stiffness for node2d and vertex2d models.
Each point represents one replicate run (with horizontal jitter to separate overlaps).
"""

import os
import sys
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from analysis_utils import setup_style
from simple_crypt_count import count_crypts_simple_method, load_final_outline, load_boundary_from_vtu

NODE2D_DIR  = '/home/orlando/Thesis/sim_results/2d_sim_04042026/16564369_node2d_2026-04-04_12-50-05'
VERTEX2D_DIR = '/home/orlando/Thesis/sim_results/2d_sim_04042026/16564370_vertex2d_2026-04-04_12-50-05'
VERTEX2D_SUMMARY = '/home/orlando/Thesis/sim_results/2d_sim_04042026/plots_vertex2d/vertex2d_summary.csv'
OUT_DIR = '/home/orlando/Thesis/sim_results/2d_sim_04042026/plots_vertex2d'

RUN_RE = re.compile(r'^g(\d+)_r(\d+)$')

# --- Font sizes ---
FS_LABEL  = 12   # axis labels
FS_TITLE  = 13   # subplot titles
FS_SUPTITLE = 14 # figure suptitle
FS_TICK   = 12   # tick labels
FS_LEGEND = 12   # legend text


def collect_node2d_crypt_counts(base_dir):
    """Run crypt counting on final node2d outlines. Returns list of (stiffness, num_crypts)."""
    records = []
    run_dirs = sorted(os.listdir(base_dir))
    total = sum(1 for d in run_dirs if RUN_RE.match(d))
    n = 0
    for run_name in run_dirs:
        m = RUN_RE.match(run_name)
        if not m:
            continue
        stiffness = int(m.group(1))
        run_path = os.path.join(base_dir, run_name, 'CryptBudding')
        if not os.path.isdir(run_path):
            continue
        # find hash subdir
        hash_dirs = [d for d in os.listdir(run_path) if os.path.isdir(os.path.join(run_path, d))]
        if not hash_dirs:
            continue
        data_dir = os.path.join(run_path, hash_dirs[0], 'results_from_time_0')
        boundary, _ = load_final_outline(data_dir)
        if boundary is None or len(boundary) < 3:
            print(f"  Skipping {run_name}: no boundary loaded")
            continue
        result = count_crypts_simple_method(boundary, boundary_is_ordered=True)
        records.append({'stiffness': stiffness, 'num_crypts': result.num_crypts})
        n += 1
        print(f"  Cell-centre {n}/{total}  G={stiffness}Pa  crypts={result.num_crypts}")
    return pd.DataFrame(records)


def collect_vertex2d_crypt_counts(summary_csv):
    """Load pre-computed final crypt counts from vertex2d summary CSV."""
    df = pd.read_csv(summary_csv)
    return df[['ECM Stiffness', 'final_num_crypts']].rename(
        columns={'ECM Stiffness': 'stiffness', 'final_num_crypts': 'num_crypts'})


def jitter(x, amount=0.03):
    """Add small log-scale jitter around x for scatter visibility."""
    return x * np.exp(np.random.uniform(-amount, amount, size=len(x)))


def plot_scatter(node2d_df, vertex2d_df, out_dir):
    setup_style()
    os.makedirs(out_dir, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(7, 3), sharey=True)
    datasets = [
        ('Cell-centre 2D',   node2d_df,   axes[0], '#2e86ab'),
        ('Vertex 2D', vertex2d_df, axes[1], '#e84855'),
    ]

    np.random.seed(42)
    for model_name, df, ax, colour in datasets:
        stiffnesses = sorted(df['stiffness'].unique())
        for s in stiffnesses:
            vals = df.loc[df['stiffness'] == s, 'num_crypts'].values.astype(float)
            xs = jitter(np.full(len(vals), s))
            ax.scatter(xs, vals, color=colour, alpha=0.65, s=40, linewidths=0,
                       zorder=3)
            # mean marker
            ax.scatter([s], [vals.mean()], color=colour, s=80, marker='D',
                       edgecolors='black', linewidths=0.8, zorder=4, alpha=0.9)

        ax.set_xlabel('ECM Stiffness (Pa)', fontsize=FS_LABEL)
        ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        ax.grid(True, which='both', alpha=0.25)
        ax.set_title(model_name, fontsize=FS_TITLE)
        # x-ticks at actual stiffness values
        ax.set_xticks(stiffnesses)
        ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
        ax.tick_params(axis='x', rotation=45, labelsize=FS_TICK)
        ax.tick_params(axis='y', labelsize=FS_TICK)

    axes[0].set_ylabel('Final Crypt Count', fontsize=FS_LABEL)

    # fig.suptitle('Final Crypt Count vs ECM Stiffness', fontsize=FS_SUPTITLE, y=1.01)
    # legend
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='grey',
               markersize=7, alpha=0.7, label='Individual replicate'),
        Line2D([0], [0], marker='D', color='w', markerfacecolor='grey',
               markeredgecolor='black', markersize=8, label='Mean'),
    ]
    axes[0].legend(handles=handles, loc='upper right', fontsize=FS_LEGEND)

    fig.tight_layout()
    out_path = os.path.join(out_dir, 'crypt_count_scatter_by_stiffness.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'\nSaved: {out_path}')

    # combined single-panel version
    fig2, ax2 = plt.subplots(figsize=(5, 3))
    np.random.seed(42)
    colours = {'Cell-centre 2D': '#2e86ab', 'Vertex 2D': '#e84855'}
    for model_name, df, _, colour in datasets:
        stiffnesses = sorted(df['stiffness'].unique())
        first = True
        for s in stiffnesses:
            vals = df.loc[df['stiffness'] == s, 'num_crypts'].values.astype(float)
            xs = jitter(np.full(len(vals), s))
            ax2.scatter(xs, vals, color=colour, alpha=0.55, s=35, linewidths=0,
                        zorder=3, label=model_name if first else '_nolegend_')
            first = False

    ax2.set_xlabel('ECM Stiffness (Pa)', fontsize=FS_LABEL)
    ax2.set_ylabel('Final Crypt Count', fontsize=FS_LABEL)
    ax2.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax2.grid(True, which='both', alpha=0.25)
    all_stiffnesses = sorted(set(node2d_df['stiffness'].tolist() + vertex2d_df['stiffness'].tolist()))
    ax2.set_xticks(all_stiffnesses)
    ax2.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax2.tick_params(axis='x', rotation=45, labelsize=FS_TICK)
    ax2.tick_params(axis='y', labelsize=FS_TICK)
    ax2.legend(fontsize=FS_LEGEND)
    ax2.set_title('Final Crypt Count vs ECM Stiffness', fontsize=FS_TITLE)
    fig2.tight_layout()
    out_path2 = os.path.join(out_dir, 'crypt_count_scatter_combined.png')
    fig2.savefig(out_path2, dpi=150, bbox_inches='tight')
    plt.close(fig2)
    print(f'Saved: {out_path2}')


if __name__ == '__main__':
    print('Collecting node2d crypt counts...')
    node2d_df = collect_node2d_crypt_counts(NODE2D_DIR)

    print('\nLoading vertex2d crypt counts from summary CSV...')
    vertex2d_df = collect_vertex2d_crypt_counts(VERTEX2D_SUMMARY)

    print(f'\nnode2d runs: {len(node2d_df)}, vertex2d runs: {len(vertex2d_df)}')
    print('\nPlotting...')
    plot_scatter(node2d_df, vertex2d_df, OUT_DIR)
