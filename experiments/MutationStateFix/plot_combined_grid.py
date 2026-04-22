#!/usr/bin/env python3
"""Combined 2×4 grid of all MutationStateFix figures.

Cols 0-1 — Vertex 2D      : row 0 = cell counts (300 Pa | 1300 Pa)
                             row 1 = quiescent fraction (300 Pa | 1300 Pa)
Cols 2-3 — Cell-centre 2D : same structure

Within each row/model-type pair the two stiffness panels share a y-axis.

All 8 panels share the same physical size.  One PNG and one SVG are saved.

Usage:
    python3 plot_combined_grid.py
    python3 plot_combined_grid.py --vertex-root /path/to/v --node-root /path/to/n
"""

import os
import glob
import argparse
import warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker
import vtk
from vtk.util.numpy_support import vtk_to_numpy

# ── paths ─────────────────────────────────────────────────────────────────────
ARCHIVE_BASE = '/home/orlando/Thesis/sim_results/mutation_state_fix_2026-04-10_14-15-05'
VERTEX_ROOT  = os.path.join(ARCHIVE_BASE, 'vertex2d')
NODE_ROOT    = os.path.join(ARCHIVE_BASE, 'node2d')

STIFFNESSES = [300, 1300]
DT          = 0.0005   # hours per simulation step

# Mutation state IDs
MUTATION_STATES = {0: 'SC', 4: 'TA', 5: 'PC', 6: 'EC'}
CELL_FULL_NAMES = {
    'SC': 'Stem Cell (SC)',
    'TA': 'Transit-Amplifying (TA)',
    'PC': 'Paneth Cell (PC)',
    'EC': 'Enterocyte (EC)',
}
CELL_COLOURS = {
    'SC': '#c8c5bd',
    'TA': '#ff0000',
    'PC': '#00ff00',
    'EC': '#0000ff',
}

# Each panel: 3.5 in wide × 3.5 in tall → total figure 14 × 7 inches
PANEL_W = 3.5
PANEL_H = 3.5

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), 'plots')


# ── style ─────────────────────────────────────────────────────────────────────
def setup_style():
    plt.rcParams.update({
        'font.family':       'serif',
        'mathtext.fontset':  'dejavuserif',
        'figure.dpi':        300,
        'savefig.dpi':       300,
        'font.size':         17,
        'axes.titlesize':    18,
        'axes.labelsize':    16,
        'legend.fontsize':   14,
        'xtick.labelsize':   14,
        'ytick.labelsize':   14,
        'lines.linewidth':   2.0,
        'axes.grid':         True,
        'grid.alpha':        0.3,
    })


# ── VTU helpers ───────────────────────────────────────────────────────────────
def _get_step(path):
    return int(os.path.basename(path).replace('results_', '').replace('.vtu', ''))


def find_results_dirs(root, stiffness):
    pattern = os.path.join(root, f'g{stiffness}_r*', 'CryptBudding',
                           '*', 'results_from_time_*')
    dirs = sorted(glob.glob(pattern))
    best = {}
    for d in dirs:
        run_key = os.path.dirname(d)
        phase   = int(os.path.basename(d).replace('results_from_time_', ''))
        if run_key not in best or phase > best[run_key][0]:
            best[run_key] = (phase, d)
    return [v[1] for v in sorted(best.values())]


# ── cell-count loading ────────────────────────────────────────────────────────
def load_count_timeseries(results_dir, node_model=False):
    vtu_files = sorted(glob.glob(os.path.join(results_dir, 'results_*.vtu')),
                       key=_get_step)
    if not vtu_files:
        return None, None

    reader = vtk.vtkXMLUnstructuredGridReader()
    times_list, counts_list = [], []

    for vtu in vtu_files:
        reader.SetFileName(vtu)
        reader.Update()
        mesh    = reader.GetOutput()
        data_obj = mesh.GetPointData() if node_model else mesh.GetCellData()
        arr = data_obj.GetArray('Mutation states')
        if arr is None:
            continue
        states = vtk_to_numpy(arr)
        row    = [np.sum(states == sid) for sid in sorted(MUTATION_STATES.keys())]
        times_list.append(_get_step(vtu) * DT)
        counts_list.append(row)

    if not times_list:
        return None, None
    return np.array(times_list), np.array(counts_list, dtype=float)


def load_count_group(root, stiffness, node_model=False):
    dirs = find_results_dirs(root, stiffness)
    if not dirs:
        print(f'  WARNING: no results dirs for g{stiffness} in {root}')
        return None, None, None, 0

    all_times, all_counts = [], []
    for d in dirs:
        t, c = load_count_timeseries(d, node_model=node_model)
        if t is not None:
            all_times.append(t)
            all_counts.append(c)
        else:
            print(f'  SKIP {d}: no VTU data')

    if not all_times:
        return None, None, None, 0

    ref_idx = max(range(len(all_times)), key=lambda i: len(all_times[i]))
    ref_t   = all_times[ref_idx]
    n_types = all_counts[0].shape[1]
    stacked = np.full((len(all_times), len(ref_t), n_types), np.nan)

    for i, (t, c) in enumerate(zip(all_times, all_counts)):
        for col in range(n_types):
            stacked[i, :, col] = (np.interp(ref_t, t, c[:, col])
                                  if len(t) > 1 else c[0, col])

    return ref_t, np.nanmean(stacked, axis=0), np.nanstd(stacked, axis=0), len(all_times)


# ── quiescent-fraction loading ────────────────────────────────────────────────
def load_quiescent_timeseries(results_dir, node_model=False):
    vtu_files = sorted(glob.glob(os.path.join(results_dir, 'results_*.vtu')),
                       key=_get_step)
    if not vtu_files:
        return None, None

    reader = vtk.vtkXMLUnstructuredGridReader()
    times_list, fracs_list = [], []

    for vtu in vtu_files:
        reader.SetFileName(vtu)
        reader.Update()
        mesh     = reader.GetOutput()
        data_obj = mesh.GetPointData() if node_model else mesh.GetCellData()
        mut_arr = data_obj.GetArray('Mutation states')
        ci_arr  = data_obj.GetArray('Contact inhibited')
        if mut_arr is None or ci_arr is None:
            continue
        states    = vtk_to_numpy(mut_arr).astype(int)
        quiescent = vtk_to_numpy(ci_arr) > 0.5
        row = []
        for sid in sorted(MUTATION_STATES.keys()):
            mask   = states == sid
            n_type = mask.sum()
            row.append(np.sum(quiescent[mask]) / n_type if n_type > 0 else np.nan)
        times_list.append(_get_step(vtu) * DT)
        fracs_list.append(row)

    if not times_list:
        return None, None
    return np.array(times_list), np.array(fracs_list, dtype=float)


def load_quiescent_group(root, stiffness, node_model=False):
    dirs = find_results_dirs(root, stiffness)
    if not dirs:
        print(f'  WARNING: no results dirs for g{stiffness} in {root}')
        return None, None, None, 0

    all_times, all_fracs = [], []
    for d in dirs:
        t, f = load_quiescent_timeseries(d, node_model=node_model)
        if t is not None:
            all_times.append(t)
            all_fracs.append(f)
        else:
            print(f'  SKIP {d}: no VTU data')

    if not all_times:
        return None, None, None, 0

    ref_idx = max(range(len(all_times)), key=lambda i: len(all_times[i]))
    ref_t   = all_times[ref_idx]
    n_types = all_fracs[0].shape[1]
    stacked = np.full((len(all_times), len(ref_t), n_types), np.nan)

    for i, (t, f) in enumerate(zip(all_times, all_fracs)):
        for col in range(n_types):
            col_data = f[:, col]
            valid    = ~np.isnan(col_data)
            if valid.sum() > 1:
                stacked[i, :, col] = np.interp(ref_t, t[valid], col_data[valid],
                                               left=np.nan, right=np.nan)
            elif valid.sum() == 1:
                stacked[i, :, col] = col_data[valid][0]

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return ref_t, np.nanmean(stacked, axis=0), np.nanstd(stacked, axis=0), len(all_times)


# ── per-axis drawing ──────────────────────────────────────────────────────────
def _draw_counts(ax, times, mean, std, title):
    for col_i, (sid, short) in enumerate(sorted(MUTATION_STATES.items())):
        colour = CELL_COLOURS[short]
        y_mean = mean[:, col_i]
        y_std  = std[:, col_i]
        ax.plot(times, y_mean, color=colour,
                label=CELL_FULL_NAMES[short], linewidth=1.8)
        ax.fill_between(times,
                        np.maximum(y_mean - y_std, 0),
                        y_mean + y_std,
                        color=colour, alpha=0.18)
    ax.set_title(title, pad=6)
    ax.set_xlabel('Time (h)')


def _draw_quiescent(ax, times, mean, std, title):
    for col_i, sid in enumerate(sorted(MUTATION_STATES.keys())):
        short  = MUTATION_STATES[sid]
        colour = CELL_COLOURS[short]
        y_mean = mean[:, col_i]
        y_std  = std[:, col_i]
        if np.all(np.isnan(y_mean)):
            continue
        ax.plot(times, y_mean, color=colour,
                label=CELL_FULL_NAMES[short], linewidth=1.8)
        ax.fill_between(times,
                        np.clip(y_mean - y_std, 0, 1),
                        np.clip(y_mean + y_std, 0, 1),
                        color=colour, alpha=0.18)
    ax.set_title(title, pad=6)
    ax.set_xlabel('Time (h)')
    ax.set_ylim(-0.05, 1.05)
    ax.yaxis.set_major_formatter(
        matplotlib.ticker.PercentFormatter(xmax=1, decimals=0))


def _no_data(ax, title):
    ax.text(0.5, 0.5, 'No data', transform=ax.transAxes,
            ha='center', va='center', color='grey')
    ax.set_title(title, pad=6)


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vertex-root', default=VERTEX_ROOT)
    parser.add_argument('--node-root',   default=NODE_ROOT)
    parser.add_argument('--output-dir',  default=OUTPUT_DIR)
    args = parser.parse_args()

    setup_style()
    os.makedirs(args.output_dir, exist_ok=True)

    # 2 rows × 5 GridSpec columns:
    #   [vertex300 | vertex1300 | <gap> | node300 | node1300]
    # Row 0 = cell counts, Row 1 = quiescent fraction
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(PANEL_W * 4, PANEL_H * 2))
    gs = GridSpec(
        2, 5,
        figure=fig,
        width_ratios=[1, 1, 0.3, 1, 1],
        wspace=0.05,   # tight within pairs; gap column provides inter-pair spacing
        hspace=0.10,
    )
    # gs columns: 0=vertex300, 1=vertex1300, 2=spacer, 3=node300, 4=node1300
    gs_cols = [0, 1, 3, 4]   # actual data columns (skipping spacer at index 2)
    axes = np.empty((2, 4), dtype=object)
    for r in range(2):
        for i, gc in enumerate(gs_cols):
            axes[r, i] = fig.add_subplot(gs[r, gc])

    # (row, col, root, node_model, stiffness, loader, drawer)
    panel_configs = [
        (0, 0, args.vertex_root, False, 300,  load_count_group,     _draw_counts),
        (0, 1, args.vertex_root, False, 1300, load_count_group,     _draw_counts),
        (1, 0, args.vertex_root, False, 300,  load_quiescent_group, _draw_quiescent),
        (1, 1, args.vertex_root, False, 1300, load_quiescent_group, _draw_quiescent),
        (0, 2, args.node_root,   True,  300,  load_count_group,     _draw_counts),
        (0, 3, args.node_root,   True,  1300, load_count_group,     _draw_counts),
        (1, 2, args.node_root,   True,  300,  load_quiescent_group, _draw_quiescent),
        (1, 3, args.node_root,   True,  1300, load_quiescent_group, _draw_quiescent),
    ]

    for row_i, col_i, root, node_model, stiff, loader, drawer in panel_configs:
        ax = axes[row_i, col_i]
        label = 'Vertex 2D' if not node_model else 'Cell-centre 2D'
        print(f'  Loading {label} g{stiff} Pa ({loader.__name__}) ...', end=' ', flush=True)
        t, mean, std, n = loader(root, stiff, node_model=node_model)
        title = f'{stiff} Pa'
        if t is None:
            print('SKIP')
            _no_data(ax, title)
        else:
            print(f'{n} replicates, {len(t)} timesteps')
            drawer(ax, t, mean, std, title)

    # Share y-axes within each stiffness pair (same row, same model)
    for row_i in range(2):
        axes[row_i, 0].get_shared_y_axes().join(axes[row_i, 0], axes[row_i, 1])
        axes[row_i, 2].get_shared_y_axes().join(axes[row_i, 2], axes[row_i, 3])
        for col_i in (1, 3):
            axes[row_i, col_i].tick_params(labelleft=False)
            axes[row_i, col_i].set_ylabel('')

    # Share x-axes within each column (counts top, quiescent bottom)
    for col_i in range(4):
        axes[0, col_i].get_shared_x_axes().join(axes[0, col_i], axes[1, col_i])
        axes[0, col_i].tick_params(labelbottom=False)
        axes[0, col_i].set_xlabel('')
        # Stiffness label only on top row; remove from bottom row
        axes[1, col_i].set_title('')

    # Row labels (left of col 0): replace individual y-axis titles
    for row_i, row_label in enumerate(['Cell counts', 'Proportion Quiescent']):
        axes[row_i, 0].annotate(
            row_label,
            xy=(0, 0.5), xycoords='axes fraction',
            xytext=(-0.30, 0.5), textcoords='axes fraction',
            ha='right', va='center',
            fontsize=plt.rcParams['axes.labelsize'],
            fontweight='bold',
            rotation=90,
        )

    # Column group labels (above cols 0-1 and cols 2-3)
    for col_i, group_label in [(0, 'Vertex 2D'), (2, 'Cell-centre 2D')]:
        ax_left  = axes[0, col_i]
        ax_right = axes[0, col_i + 1]
        # midpoint in figure coordinates between the two axes
        fig.canvas.draw()  # needed to resolve axes positions
        x_left  = ax_left.get_position().x0
        x_right = ax_right.get_position().x1
        x_mid   = (x_left + x_right) / 2
        y_top   = ax_left.get_position().y1
        fig.text(x_mid, y_top + 0.05, group_label,
                 ha='center', va='bottom',
                 fontsize=plt.rcParams['axes.titlesize'],
                 fontweight='bold')

    # Single shared legend from the first axes that has handles
    handles, labels_leg = None, None
    for ax in axes.flat:
        h, l = ax.get_legend_handles_labels()
        if h:
            handles, labels_leg = h, l
            break

    if handles:
        fig.legend(handles, labels_leg,
                   loc='lower center',
                   bbox_to_anchor=(0.46, 0),
                   ncol=4, frameon=True,
                   fontsize=plt.rcParams['legend.fontsize'],
                   borderpad=0.8, labelspacing=0.6, handlelength=2.0,
                   columnspacing=1.5)
        fig.subplots_adjust(bottom=0.18)

    stem = os.path.join(args.output_dir, 'combined_grid')
    for ext in ('svg', 'png'):
        out = f'{stem}.{ext}'
        fig.savefig(out, bbox_inches='tight')
        print(f'Saved: {out}')
    plt.close(fig)
    print('Done.')


if __name__ == '__main__':
    main()
