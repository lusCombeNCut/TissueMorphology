#!/usr/bin/env python3
"""Plot mean cell type counts over time for 300 Pa and 1300 Pa replicates.

Parses VTU files directly using the 'Mutation states' array.
  Vertex 2D: mutation states live in CellData
  Node 2D:   mutation states live in PointData

Mutation state IDs: SC=0, TA=4, PC=5, EC=6

Two figures: one for Vertex 2D, one for Node 2D.
Each figure has two subplots (300 Pa | 1300 Pa), each showing 4 coloured
lines (one per cell type) + shaded ±1 SD bands.

Vertex 2D data: /home/orlando/Thesis/sim_results/2d_sim_06042026/2d_sim_0904_Vertex_mutations
Node 2D data:   /home/orlando/Thesis/sim_results/2d_sim_06042026/16571895_node2d_2026-04-06_11-56-52
"""

import os
import glob
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy

# ── paths ─────────────────────────────────────────────────────────────────────
# Archive structure: mutation_state_fix_<TIMESTAMP>/{vertex2d,node2d}/g{Pa}_r{rep}/...
ARCHIVE_BASE = '/home/orlando/Thesis/sim_results/mutation_state_fix_2026-04-10_14-15-05'
VERTEX_ROOT  = os.path.join(ARCHIVE_BASE, 'vertex2d')
NODE_ROOT    = os.path.join(ARCHIVE_BASE, 'node2d')

STIFFNESSES = [300, 1300]   # Pa values to plot
DT          = 0.0005        # hours per timestep (step * DT = hours)

# Mutation state IDs → display label
MUTATION_STATES = {0: 'SC', 4: 'TA', 5: 'PC', 6: 'EC'}
CELL_FULL_NAMES = {
    'SC': 'Stem Cell (SC)',
    'TA': 'Transit-Amplifying (TA)',
    'PC': 'Paneth Cell (PC)',
    'EC': 'Enterocyte (EC)',
}
CELL_COLOURS = {
    'SC': '#87ffed',
    'TA': '#ff0000',
    'PC': '#00ff00',
    'EC': '#0000ff',
}

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), 'plots')


# ── style ─────────────────────────────────────────────────────────────────────
def setup_style():
    plt.rcParams.update({
        'font.family':                   'serif',
        'mathtext.fontset':              'dejavuserif',
        'figure.dpi':                    300,
        'savefig.dpi':                   300,
        'savefig.bbox':                  'tight',
        'font.size':                     12,
        'axes.titlesize':                13,
        'axes.labelsize':                12,
        'legend.fontsize':               10,
        'xtick.labelsize':               10,
        'ytick.labelsize':               10,
        'lines.linewidth':               1.8,
        'axes.grid':                     True,
        'grid.alpha':                    0.3,
    })


# ── VTU parsing ───────────────────────────────────────────────────────────────
def _get_step(path):
    return int(os.path.basename(path).replace('results_', '').replace('.vtu', ''))


def load_mutation_timeseries(results_dir, node_model=False):
    """Parse all VTU files in results_dir.

    Returns
    -------
    times : ndarray (T,) in hours
    counts : ndarray (T, 4) — columns ordered SC, TA, PC, EC
    """
    vtu_files = sorted(
        glob.glob(os.path.join(results_dir, 'results_*.vtu')),
        key=_get_step,
    )
    if not vtu_files:
        return None, None

    reader = vtk.vtkXMLUnstructuredGridReader()

    times_list  = []
    counts_list = []

    for vtu in vtu_files:
        reader.SetFileName(vtu)
        reader.Update()
        mesh = reader.GetOutput()

        if node_model:
            data_obj = mesh.GetPointData()
        else:
            data_obj = mesh.GetCellData()

        arr = data_obj.GetArray('Mutation states')
        if arr is None:
            continue

        states = vtk_to_numpy(arr)
        row = [np.sum(states == sid) for sid in sorted(MUTATION_STATES.keys())]
        times_list.append(_get_step(vtu) * DT)
        counts_list.append(row)

    if not times_list:
        return None, None

    return np.array(times_list), np.array(counts_list, dtype=float)


def find_results_dirs(root, stiffness):
    """Glob all results_from_time_* directories for a stiffness level."""
    pattern = os.path.join(root, f'g{stiffness}_r*', 'CryptBudding',
                           '*', 'results_from_time_*')
    dirs = sorted(glob.glob(pattern))
    # If multiple phases, keep only the highest-numbered per replicate run
    best = {}
    for d in dirs:
        run_key = os.path.dirname(d)   # up to the hash dir
        phase = int(os.path.basename(d).replace('results_from_time_', ''))
        if run_key not in best or phase > best[run_key][0]:
            best[run_key] = (phase, d)
    return [v[1] for v in sorted(best.values())]


def load_group(root, stiffness, node_model=False):
    """Load all replicates for one stiffness, return aggregated arrays.

    Returns
    -------
    times : ndarray (T,)
    mean  : ndarray (T, 4)
    std   : ndarray (T, 4)
    n     : int
    """
    dirs = find_results_dirs(root, stiffness)
    if not dirs:
        print(f'  WARNING: no results dirs found for g{stiffness} in {root}')
        return None, None, None, 0

    all_times, all_counts = [], []
    for d in dirs:
        t, c = load_mutation_timeseries(d, node_model=node_model)
        if t is not None:
            all_times.append(t)
            all_counts.append(c)
        else:
            print(f'  SKIP {d}: no VTU data')

    if not all_times:
        return None, None, None, 0

    # Interpolate onto the longest time grid
    ref_idx = max(range(len(all_times)), key=lambda i: len(all_times[i]))
    ref_t   = all_times[ref_idx]
    n_types = all_counts[0].shape[1]
    stacked = np.full((len(all_times), len(ref_t), n_types), np.nan)

    for i, (t, c) in enumerate(zip(all_times, all_counts)):
        for col in range(n_types):
            if len(t) > 1:
                stacked[i, :, col] = np.interp(ref_t, t, c[:, col])
            else:
                stacked[i, :, col] = c[0, col]

    return ref_t, np.nanmean(stacked, axis=0), np.nanstd(stacked, axis=0), len(all_times)


# ── plotting ──────────────────────────────────────────────────────────────────
def _fill_ax(ax, times, mean, std, title):
    """Draw 4 cell-type lines with SD fills on ax."""
    for col_i, (sid, label) in enumerate(sorted(MUTATION_STATES.items())):
        short = MUTATION_STATES[sid]
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
    ax.set_ylabel('Cell count')


def make_figure(root, model_label, node_model, output_path):
    """Build and save one figure (two subplots: 300 Pa | 1300 Pa)."""
    group_data = {}
    for g in STIFFNESSES:
        print(f'  Loading {model_label} g{g} Pa ...', end=' ', flush=True)
        times, mean, std, n = load_group(root, g, node_model=node_model)
        if times is None:
            print('SKIP')
            continue
        group_data[g] = (times, mean, std, n)
        print(f'{n} replicates, {len(times)} timesteps')

    if not group_data:
        print(f'  ERROR: no data for {model_label}')
        return

    n_panels = len(group_data)
    fig, axes = plt.subplots(1, n_panels, figsize=(4 * n_panels, 3),
                             sharey=False, sharex=True)
    if n_panels == 1:
        axes = [axes]

    for ax, g in zip(axes, sorted(group_data)):
        times, mean, std, n = group_data[g]
        _fill_ax(ax, times, mean, std, title=f'{g} Pa  (n={n})')

    # shared legend from first axes
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels,
               loc='lower center',
               bbox_to_anchor=(0.5, -0.12),
               ncol=4,
               frameon=True)
    # fig.suptitle(f'{model_label} — mean cell type counts (300 Pa vs 1300 Pa)',
    #              fontsize=14)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    plt.savefig(output_path, format='svg', bbox_inches='tight')
    print(f'  Saved: {output_path}')
    plt.savefig(output_path.replace('.svg', '.png'), dpi=300, bbox_inches='tight')
    print(f'  Saved: {output_path.replace(".svg", ".png")}')
    plt.close(fig)


# ── entry point ───────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', default=OUTPUT_DIR)
    parser.add_argument('--vertex-root', default=VERTEX_ROOT)
    parser.add_argument('--node-root',   default=NODE_ROOT)
    args = parser.parse_args()

    setup_style()
    os.makedirs(args.output_dir, exist_ok=True)

    print('\nVertex 2D:')
    make_figure(
        root=args.vertex_root,
        model_label='Vertex 2D',
        node_model=False,
        output_path=os.path.join(args.output_dir, 'vertex2d_cell_type_counts.svg'),
    )

    print('\nCell-centre 2D:')
    make_figure(
        root=args.node_root,
        model_label='Cell-centre 2D',
        node_model=True,
        output_path=os.path.join(args.output_dir, 'cell_centre_2d_cell_type_counts.svg'),
    )

    print('\nDone.')


if __name__ == '__main__':
    main()
