#!/usr/bin/env python3
"""Plot per-cell-type quiescent fraction over time.

For each cell type (SC, TA, PC, EC) computes the fraction of cells of that
type that are contact-inhibited at each timestep:

    Q_type(t) = (# cells of type with Contact inhibited == 1) /
                (# cells of that type)

Reads VTU files using:
  - 'Mutation states'    array (true type, fixed by CellTrueMutationStateWriter)
  - 'Contact inhibited'  array (1 if quiescent, 0 otherwise)

Requires the MutationStateFix dataset (with CellTrueMutationStateWriter).
Will NOT produce correct per-type results with old data — use
plot_contact_inhibition.py for overall quiescent fraction from .dat files.

Two figures: Vertex 2D and Node 2D.
Each figure: two subplots (300 Pa | 1300 Pa), four coloured lines per panel
             (one per cell type) with ±1 SD shaded bands across replicates.

Usage:
    python3 plot_quiescent_fraction.py
    python3 plot_quiescent_fraction.py --vertex-root /path/to/vertex --node-root /path/to/node
"""

import os
import glob
import argparse
import warnings
import numpy as np
import matplotlib
import matplotlib.ticker
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy

import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_utils import setup_style, figsize

# ── default paths (update after extracting the MutationStateFix archive) ─────
# Archive structure: mutation_state_fix_<TIMESTAMP>/{vertex2d,node2d}/g{Pa}_r{rep}/...
ARCHIVE_BASE = '/home/orlando/Thesis/sim_results/mutation_state_fix_2026-04-10_14-15-05'
VERTEX_ROOT  = os.path.join(ARCHIVE_BASE, 'vertex2d')
NODE_ROOT    = os.path.join(ARCHIVE_BASE, 'node2d')

STIFFNESSES = [300, 1300]
DT          = 0.0005   # hours per simulation step

# Mutation state IDs — must match CellTrueMutationStateWriter output
MUTATION_STATES = {0: 'SC', 4: 'TA', 5: 'PC', 6: 'EC'}
CELL_COLOURS = {
    'SC': '#87ffed',
    'TA': '#ff0000',
    'PC': '#00ff00',
    'EC': '#0000ff',
}
CELL_FULL_NAMES = {
    'SC': 'Stem Cell (SC)',
    'TA': 'Transit-Amplifying (TA)',
    'PC': 'Paneth Cell (PC)',
    'EC': 'Enterocyte (EC)',
}
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), 'plots')


# ── VTU parsing ───────────────────────────────────────────────────────────────
def _get_step(path):
    return int(os.path.basename(path).replace('results_', '').replace('.vtu', ''))


def load_quiescent_timeseries(results_dir, node_model=False):
    """Parse all VTU files in results_dir.

    Returns
    -------
    times  : ndarray (T,)  — simulation time in hours
    q_frac : ndarray (T, 4) — quiescent fraction per cell type,
             columns ordered SC, TA, PC, EC (matching sorted MUTATION_STATES)
             NaN where a cell type has zero cells at that timestep.
    """
    vtu_files = sorted(
        glob.glob(os.path.join(results_dir, 'results_*.vtu')),
        key=_get_step,
    )
    if not vtu_files:
        return None, None

    reader = vtk.vtkXMLUnstructuredGridReader()
    times_list  = []
    fracs_list  = []

    for vtu in vtu_files:
        reader.SetFileName(vtu)
        reader.Update()
        mesh = reader.GetOutput()

        data_obj = mesh.GetPointData() if node_model else mesh.GetCellData()

        mut_arr = data_obj.GetArray('Mutation states')
        ci_arr  = data_obj.GetArray('Contact inhibited')

        if mut_arr is None or ci_arr is None:
            continue

        states    = vtk_to_numpy(mut_arr).astype(int)
        quiescent = vtk_to_numpy(ci_arr) > 0.5   # bool

        row = []
        for sid in sorted(MUTATION_STATES.keys()):
            mask = states == sid
            n_type = mask.sum()
            row.append(np.sum(quiescent[mask]) / n_type if n_type > 0 else np.nan)

        times_list.append(_get_step(vtu) * DT)
        fracs_list.append(row)

    if not times_list:
        return None, None

    return np.array(times_list), np.array(fracs_list, dtype=float)


def find_results_dirs(root, stiffness):
    """Return sorted list of results_from_time_* dirs for a stiffness."""
    pattern = os.path.join(root, f'g{stiffness}_r*', 'CryptBudding',
                           '*', 'results_from_time_*')
    dirs = sorted(glob.glob(pattern))
    # Keep only the latest phase per replicate
    best = {}
    for d in dirs:
        run_key = os.path.dirname(d)
        phase   = int(os.path.basename(d).replace('results_from_time_', ''))
        if run_key not in best or phase > best[run_key][0]:
            best[run_key] = (phase, d)
    return [v[1] for v in sorted(best.values())]


def load_group(root, stiffness, node_model=False):
    """Aggregate quiescent fractions across all replicates.

    Returns
    -------
    times : ndarray (T,)
    mean  : ndarray (T, 4)
    std   : ndarray (T, 4)
    n     : int  — number of replicates loaded
    """
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

    # Interpolate onto the longest common time grid
    ref_idx  = max(range(len(all_times)), key=lambda i: len(all_times[i]))
    ref_t    = all_times[ref_idx]
    n_types  = all_fracs[0].shape[1]
    stacked  = np.full((len(all_times), len(ref_t), n_types), np.nan)

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


# ── plotting ──────────────────────────────────────────────────────────────────
def _fill_ax(ax, times, mean, std, title, n_reps):
    """Draw per-cell-type quiescent fraction lines + SD fill on ax."""
    for col_i, sid in enumerate(sorted(MUTATION_STATES.keys())):
        label  = MUTATION_STATES[sid]
        colour = CELL_COLOURS[label]
        y_mean = mean[:, col_i]
        y_std  = std[:, col_i]

        # Skip types that are NaN throughout (never present)
        if np.all(np.isnan(y_mean)):
            continue

        ax.plot(times, y_mean, color=colour,
                label=CELL_FULL_NAMES[label], linewidth=1.8)
        ax.fill_between(times,
                        np.clip(y_mean - y_std, 0, 1),
                        np.clip(y_mean + y_std, 0, 1),
                        color=colour, alpha=0.18)

    ax.set_title(f'{title}  (n={n_reps})', pad=6)
    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Quiescent fraction')
    ax.set_ylim(-0.05, 1.05)
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=1, decimals=0))


def make_figure(model_label, root, node_model=False):
    """Build and save a figure for one model (vertex2d or node2d)."""
    fig, axes = plt.subplots(1, 2, figsize=figsize(w_scale=1.4, h_scale=0.7),
                             sharey=True)

    any_data = False
    for ax, stiff in zip(axes, STIFFNESSES):
        t, mean, std, n = load_group(root, stiff, node_model=node_model)
        if t is None:
            ax.text(0.5, 0.5, 'No data', transform=ax.transAxes,
                    ha='center', va='center', color='grey')
            ax.set_title(f'{stiff} Pa', pad=6)
            continue
        _fill_ax(ax, t, mean, std, f'{stiff} Pa', n)
        any_data = True

    if not any_data:
        print(f'  No data found for {model_label}. Skipping figure.')
        plt.close(fig)
        return

    # Shared legend below both panels
    handles, labels = axes[0].get_legend_handles_labels()
    if not handles:
        handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center',
               bbox_to_anchor=(0.5, -0.18), ncol=4, frameon=True,
               fontsize=plt.rcParams['legend.fontsize'])

    # fig.suptitle(f'Per-cell-type quiescent fraction — {model_label}', y=1.02)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    stem = f'{model_label.lower().replace(" ", "_")}_quiescent_fraction'
    for ext in ('svg', 'png'):
        out = os.path.join(OUTPUT_DIR, f'{stem}.{ext}')
        fig.savefig(out, bbox_inches='tight')
        print(f'  Saved: {out}')
    plt.close(fig)


# ── entry point ───────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vertex-root', default=VERTEX_ROOT,
                        help='Root directory for vertex2d results')
    parser.add_argument('--node-root', default=NODE_ROOT,
                        help='Root directory for node2d results')
    args = parser.parse_args()

    setup_style()

    print('Vertex 2D:')
    make_figure('Vertex 2D', args.vertex_root, node_model=False)

    print('Node 2D:')
    make_figure('Node 2D', args.node_root,   node_model=True)


if __name__ == '__main__':
    main()
