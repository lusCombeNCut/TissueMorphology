#!/usr/bin/env python3
"""Shared analysis and plotting utilities for TissueMorphology parameter sweeps."""

import os
import sys
import re
import glob
import json
import csv
from collections import defaultdict

import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.cm as cm

HAS_MATPLOTLIB = True

# Default plot geometry (A4, 1-inch margins, 2 columns)
DEFAULT_FIG_WIDTH = 10   # inches
DEFAULT_FIG_HEIGHT = 8   # inches
DEFAULT_FONT_SCALE = 1.0

# Mutable globals set by add_common_args / apply_common_args
_FIG_W = DEFAULT_FIG_WIDTH
_FIG_H = DEFAULT_FIG_HEIGHT
_FONT_SCALE = DEFAULT_FONT_SCALE

_crypt_count_mod = None
for _subdir in ['ECMStiffnessSweep', 'CryptBudding', '.']:
    _candidate = os.path.join(os.path.dirname(os.path.abspath(__file__)), _subdir)
    if os.path.isfile(os.path.join(_candidate, 'simple_crypt_count.py')):
        sys.path.insert(0, _candidate)
        import simple_crypt_count as _crypt_count_mod
        break


def add_common_args(parser):
    """Add shared CLI arguments for figure sizing, font scaling, and units."""
    g = parser.add_argument_group('Plot formatting')
    g.add_argument('--fig-width', type=float, default=DEFAULT_FIG_WIDTH,
                   help='Figure width in inches (default: %(default)s)')
    g.add_argument('--fig-height', type=float, default=DEFAULT_FIG_HEIGHT,
                   help='Figure height in inches (default: %(default)s)')
    g.add_argument('--font-scale', type=float, default=DEFAULT_FONT_SCALE,
                   help='Multiply all font sizes by this factor (default: %(default)s)')

    u = parser.add_argument_group('Unit conversion overrides')
    u.add_argument('--cell-diameter', type=float, default=None,
                   help='Cell diameter in µm (overrides params.json)')
    u.add_argument('--ref-force', type=float, default=None,
                   help='Reference force in nN (overrides params.json)')
    u.add_argument('--ref-velocity', type=float, default=None,
                   help='Reference velocity in µm/h (overrides params.json)')
    return parser


def apply_common_args(args):
    """Apply shared CLI arguments to global plot settings."""
    global _FIG_W, _FIG_H, _FONT_SCALE
    _FIG_W = args.fig_width
    _FIG_H = args.fig_height
    _FONT_SCALE = args.font_scale


def figsize(w_scale=1.0, h_scale=1.0):
    """Return (width, height) tuple scaled from the global defaults."""
    return (_FIG_W * w_scale, _FIG_H * h_scale)


def setup_style():
    """Configure consistent plot aesthetics."""
    s = _FONT_SCALE
    plt.rcParams.update({
        'font.family': 'serif',
        'mathtext.fontset': 'dejavuserif',
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'figure.constrained_layout.use': True,
        'font.size': 22 * s,
        'axes.titlesize': 26 * s,
        'axes.labelsize': 24 * s,
        'legend.fontsize': 20 * s,
        'xtick.labelsize': 20 * s,
        'ytick.labelsize': 20 * s,
        'lines.linewidth': 1.8,
        'lines.markersize': 5,
        'axes.grid': True,
        'grid.alpha': 0.3,
    })


def _add_legend(ax, title=None, ncol=2):
    """Place legend below the axes with an optional title."""
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25),
              ncol=ncol, title=title, frameon=True)


def _param_label(value, unit):
    """Format a parameter value with its unit for legend entries."""
    if unit:
        return f'{value} {unit}'
    return f'{value}'


def _save_fig(output_path):
    """Save the current figure as SVG."""
    svg_path = os.path.splitext(output_path)[0] + '.svg'
    plt.savefig(svg_path, format='svg')


def get_param_colours(param_values):
    """Return a colour map for a list of parameter values."""
    n = len(param_values)
    if n <= 10:
        cmap = plt.cm.tab10
        return {v: cmap(i / max(n - 1, 1)) for i, v in enumerate(param_values)}
    cmap = plt.cm.viridis
    return {v: cmap(i / max(n - 1, 1)) for i, v in enumerate(param_values)}


def find_summary_csvs(base_dir):
    """Recursively find all crypt_summary.csv files under base_dir.

    When multiple results_from_time_* subdirectories exist, only the
    highest-numbered one is included.

    Returns list of (csv_path, params_json_path_or_None).
    """
    results = []
    for root, dirs, files in os.walk(base_dir, followlinks=True):
        if 'crypt_summary.csv' not in files:
            continue
        csv_path = os.path.join(root, 'crypt_summary.csv')

        basename = os.path.basename(root)
        m = re.match(r'results_from_time_(\d+)$', basename)
        if m:
            phase_time = int(m.group(1))
            parent = os.path.dirname(root)
            sibling_phases = [
                int(sm.group(1))
                for entry in os.listdir(parent)
                for sm in [re.match(r'results_from_time_(\d+)$', entry)]
                if sm
            ]
            if sibling_phases and phase_time < max(sibling_phases):
                continue

        params_path = None
        for parent in [root, os.path.dirname(root), os.path.dirname(os.path.dirname(root))]:
            candidate = os.path.join(parent, 'params.json')
            if os.path.isfile(candidate):
                params_path = candidate
                break
        results.append((csv_path, params_path))
    return results


def load_params_json(json_path):
    """Load a params.json file and return its contents as a dict."""
    with open(json_path, 'r') as f:
        return json.load(f)


def extract_param_from_json(params, param_path):
    """Extract a nested parameter value using a dot-separated path."""
    keys = param_path.split('.')
    val = params
    for k in keys:
        if isinstance(val, dict) and k in val:
            val = val[k]
        else:
            return None
    return val


def extract_param_from_dir(csv_path, param_name):
    """Infer the swept parameter value from the directory path."""
    path_str = csv_path.replace('\\', '/')

    m = re.search(r'/s([\d.]+)_r(\d+)/', path_str)
    if m and param_name in ('stiffness', 'ecmStiffness'):
        return float(m.group(1)), int(m.group(2))

    pattern = rf'{param_name}_([\d.eE+-]+)/run_(\d+)'
    m = re.search(pattern, path_str)
    if m:
        return float(m.group(1)), int(m.group(2))

    m = re.search(r'stiffness_([\d.]+)/run_(\d+)', path_str)
    if m:
        return float(m.group(1)), int(m.group(2))

    return None, None


def identify_model_type(csv_path, params=None):
    """Determine model type from path or params.json."""
    if params and 'simulation' in params:
        mt = params['simulation'].get('modelType')
        if mt:
            return mt
    path_lower = csv_path.lower()
    for mt in ['node2d', 'vertex2d', 'node3d', 'vertex3d']:
        if mt in path_lower:
            return mt
    return 'unknown'


def load_crypt_summary(csv_path):
    """Load crypt_summary.csv into a dict of numpy arrays."""
    data = defaultdict(list)
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            for key, val in row.items():
                key = key.strip()
                data[key].append(float(val))
    return {k: np.array(v) for k, v in data.items()}


def load_sweep_data(base_dir, param_name, param_json_path=None, model_filter=None):
    """Load all crypt_summary.csv files from a sweep and organise by
    (model_type, param_value, replicate).

    Returns:
        dict[model_type] -> dict[param_value] -> list of (run_number, data_dict)
    """
    csvs = find_summary_csvs(base_dir)
    if not csvs:
        print(f"  No crypt_summary.csv files found under {base_dir}")
        return {}

    sweep = defaultdict(lambda: defaultdict(list))

    for csv_path, params_path in csvs:
        params = load_params_json(params_path) if params_path else None
        model = identify_model_type(csv_path, params)
        if model_filter and model != model_filter:
            continue

        param_val = None
        run_num = None

        if params and param_json_path:
            param_val = extract_param_from_json(params, param_json_path)

        if param_val is None and params:
            param_val = extract_param_from_json(params,
                            'forces.ECMConfinementForce.ecmStiffness')

        if param_val is None:
            param_val, run_num = extract_param_from_dir(csv_path, param_name)

        if run_num is None and params:
            run_num = params.get('simulation', {}).get('runNumber', 0)
        if run_num is None:
            run_num = 0

        if param_val is None:
            print(f"  SKIP {csv_path}: could not determine {param_name}")
            continue

        data = load_crypt_summary(csv_path)
        data['_data_dir'] = os.path.dirname(csv_path)
        sweep[model][param_val].append((run_num, data))

    return dict(sweep)


def aggregate_timeseries(runs, column, time_col='time'):
    """Aggregate a time-series column across replicates onto a common time grid.

    Returns: times, mean, std, n
    """
    if not runs:
        return np.array([]), np.array([]), np.array([]), np.array([])

    all_times = []
    all_series = []
    for _, data in runs:
        if time_col in data and column in data:
            t = data[time_col]
            v = data[column]
            if len(t) > 0:
                all_times.append(t)
                all_series.append(v)

    if not all_times:
        return np.array([]), np.array([]), np.array([]), np.array([])

    ref_idx = max(range(len(all_times)), key=lambda i: len(all_times[i]))
    times = all_times[ref_idx]

    values = np.full((len(all_series), len(times)), np.nan)
    for i, (t, v) in enumerate(zip(all_times, all_series)):
        if len(t) > 1:
            values[i, :] = np.interp(times, t, v)
        elif len(t) == 1:
            values[i, :] = v[0]

    mean = np.nanmean(values, axis=0)
    std = np.nanstd(values, axis=0)
    n = np.sum(~np.isnan(values), axis=0)

    return times, mean, std, n


def list_vtu_files(data_dir):
    """Find and sort VTU files by timestep number."""
    vtu_files = glob.glob(os.path.join(data_dir, 'results_*.vtu'))
    def get_step(p):
        base = os.path.basename(p)
        return int(base.replace('results_', '').replace('.vtu', ''))
    return sorted(vtu_files, key=get_step)


def parse_vtu_cell_data(vtu_path, keys=None):
    """Parse cell data arrays from a VTU file using VTK.
    Returns dict mapping data array name -> numpy array.
    """
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_path)
    reader.Update()
    mesh = reader.GetOutput()
    result = {}
    cd = mesh.GetCellData()
    for i in range(cd.GetNumberOfArrays()):
        name = cd.GetArrayName(i)
        if keys is None or name in keys:
            result[name] = vtk_to_numpy(cd.GetArray(i))
    pd = mesh.GetPointData()
    for i in range(pd.GetNumberOfArrays()):
        name = pd.GetArrayName(i)
        if keys is None or name in keys:
            result[name] = vtk_to_numpy(pd.GetArray(i))
    return result


def get_vtu_timestep_data(data_dir, column_key, dt=None):
    """Extract a per-cell data column from each VTU timestep.

    Returns: times, values_per_timestep (list of arrays)
    """
    vtu_files = list_vtu_files(data_dir)
    if not vtu_files:
        return [], []

    times = []
    values = []
    for vtu in vtu_files:
        step = int(os.path.basename(vtu).replace('results_', '').replace('.vtu', ''))
        data = parse_vtu_cell_data(vtu, keys=[column_key])
        if column_key in data:
            t = step * dt if dt else float(step)
            times.append(t)
            values.append(data[column_key])

    return times, values


def list_ghost_vtp_files(data_dir):
    """Find and sort ghost ECM VTP files by timestep."""
    vtp_files = glob.glob(os.path.join(data_dir, 'ghost_ecm_*.vtp'))
    def get_step(p):
        base = os.path.basename(p)
        return int(base.replace('ghost_ecm_', '').replace('.vtp', ''))
    return sorted(vtp_files, key=get_step)


def parse_ghost_vtp(vtp_path):
    """Parse ghost node data from a VTP file using VTK."""
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtp_path)
    reader.Update()
    mesh = reader.GetOutput()
    n = mesh.GetNumberOfPoints()
    result = {'n_nodes': n}
    if n > 0:
        result['positions'] = vtk_to_numpy(mesh.GetPoints().GetData())
    pd = mesh.GetPointData()
    for i in range(pd.GetNumberOfArrays()):
        name = pd.GetArrayName(i)
        result[name] = vtk_to_numpy(pd.GetArray(i))
    return result


def get_ghost_node_timeseries(data_dir, dt=None):
    """Extract ghost node summary statistics over time from VTP files."""
    vtp_files = list_ghost_vtp_files(data_dir)
    if not vtp_files:
        return {}

    first_data = parse_ghost_vtp(vtp_files[0])
    init_positions = first_data.get('positions')

    times = []
    n_nodes = []
    mean_density = []
    mean_anisotropy = []
    mean_strain = []
    mean_displacement = []

    for vtp in vtp_files:
        step = int(os.path.basename(vtp).replace('ghost_ecm_', '').replace('.vtp', ''))
        t = step * dt if dt else float(step)
        data = parse_ghost_vtp(vtp)

        times.append(t)
        n_nodes.append(data['n_nodes'])

        if 'density' in data and len(data['density']) > 0:
            mean_density.append(np.mean(data['density']))
        else:
            mean_density.append(np.nan)

        if 'anisotropy' in data and len(data['anisotropy']) > 0:
            mean_anisotropy.append(np.mean(data['anisotropy']))
        else:
            mean_anisotropy.append(np.nan)

        if 'rest_length_strain' in data and len(data['rest_length_strain']) > 0:
            mean_strain.append(np.mean(data['rest_length_strain']))
        else:
            mean_strain.append(np.nan)

        if 'positions' in data and init_positions is not None:
            pos = data['positions']
            n_compare = min(len(pos), len(init_positions))
            if n_compare > 0:
                disp = np.linalg.norm(pos[:n_compare] - init_positions[:n_compare], axis=1)
                mean_displacement.append(np.mean(disp))
            else:
                mean_displacement.append(np.nan)
        else:
            mean_displacement.append(np.nan)

    result = {
        'times': np.array(times),
        'n_nodes': np.array(n_nodes),
        'mean_density': np.array(mean_density),
        'mean_anisotropy': np.array(mean_anisotropy),
    }
    if not all(np.isnan(mean_strain)):
        result['mean_strain'] = np.array(mean_strain)
    if not all(np.isnan(mean_displacement)):
        result['mean_displacement'] = np.array(mean_displacement)
    return result


def compute_contact_inhibition_timeseries(data_dir, dt=None):
    """Compute fraction of contact-inhibited cells at each VTU timestep."""
    vtu_files = list_vtu_files(data_dir)
    times = []
    fractions = []
    for vtu in vtu_files:
        step = int(os.path.basename(vtu).replace('results_', '').replace('.vtu', ''))
        data = parse_vtu_cell_data(vtu, keys=['Contact inhibited'])
        if 'Contact inhibited' in data:
            arr = data['Contact inhibited']
            t = step * dt if dt else float(step)
            times.append(t)
            fractions.append(np.mean(arr > 0.5) if len(arr) > 0 else 0.0)
    return np.array(times), np.array(fractions)


CELL_TYPE_NAMES = {0: 'Stem', 1: 'Transit', 2: 'Paneth', 3: 'Enterocyte'}


def compute_cell_type_timeseries(data_dir, dt=None):
    """Compute cell type fractions at each VTU timestep."""
    vtu_files = list_vtu_files(data_dir)
    times = []
    type_fracs = defaultdict(list)

    for vtu in vtu_files:
        step = int(os.path.basename(vtu).replace('results_', '').replace('.vtu', ''))
        data = parse_vtu_cell_data(vtu, keys=['cell_type_id'])
        if 'cell_type_id' not in data:
            continue
        arr = data['cell_type_id']
        t = step * dt if dt else float(step)
        times.append(t)
        n = len(arr)
        for type_id, name in CELL_TYPE_NAMES.items():
            frac = np.sum(np.abs(arr - type_id) < 0.5) / n if n > 0 else 0.0
            type_fracs[name].append(frac)

    result = {'times': np.array(times)}
    for name in CELL_TYPE_NAMES.values():
        result[name] = np.array(type_fracs[name])
    return result


def compute_lumen_force_timeseries(data_dir, dt=None):
    """Compute mean lumen force magnitude at each VTU timestep."""
    vtu_files = list_vtu_files(data_dir)
    times = []
    mean_forces = []

    for vtu in vtu_files:
        step = int(os.path.basename(vtu).replace('results_', '').replace('.vtu', ''))
        data = parse_vtu_cell_data(vtu,
                   keys=['lumen_force_x', 'lumen_force_y', 'lumen_force_z'])
        fx = data.get('lumen_force_x')
        fy = data.get('lumen_force_y')
        if fx is not None and fy is not None:
            fz = data.get('lumen_force_z', np.zeros_like(fx))
            mag = np.sqrt(fx**2 + fy**2 + fz**2)
            t = step * dt if dt else float(step)
            times.append(t)
            mean_forces.append(np.mean(mag))

    return np.array(times), np.array(mean_forces)


def compute_crypt_count_timeseries(data_dir, model_type, max_timesteps=30):
    """Compute crypt count and circularity at sampled timesteps for one run."""
    if _crypt_count_mod is None or model_type not in ('node2d', 'vertex2d'):
        return {}

    csv_path = os.path.join(data_dir, 'crypt_summary.csv')
    if not os.path.isfile(csv_path):
        return {}
    csv_data = load_crypt_summary(csv_path)
    csv_times = csv_data.get('time', np.array([]))
    if len(csv_times) == 0:
        return {}

    if model_type == 'node2d':
        all_files = sorted(glob.glob(os.path.join(data_dir, 'outline_*.vtp')))
        def _get_step(f):
            return int(os.path.basename(f).replace('outline_', '').replace('.vtp', ''))
    else:
        all_files = sorted(glob.glob(os.path.join(data_dir, 'results_*.vtu')))
        def _get_step(f):
            return int(os.path.basename(f).replace('results_', '').replace('.vtu', ''))

    all_files = [f for f in all_files if _get_step(f) >= 0]
    all_files.sort(key=_get_step)
    if not all_files:
        return {}

    n = len(all_files)
    if n > max_timesteps:
        indices = sorted(set(
            [0, n - 1] +
            list(np.linspace(0, n - 1, max_timesteps, dtype=int))
        ))
    else:
        indices = list(range(n))

    result_times = []
    result_crypts = []
    result_circ = []

    for idx in indices:
        f = all_files[idx]
        if idx >= len(csv_times):
            continue
        t = csv_times[idx]

        if model_type == 'node2d':
            boundary, _ = _crypt_count_mod.load_outline_from_vtp(f)
        else:
            boundary, _ = _crypt_count_mod.load_boundary_from_vtu(f)

        result = _crypt_count_mod.count_crypts_simple_method(
            boundary, boundary_is_ordered=True)
        result_times.append(t)
        result_crypts.append(result.num_crypts)
        result_circ.append(result.circularity)

    if not result_times:
        return {}

    return {
        'time': np.array(result_times),
        'num_crypts': np.array(result_crypts, dtype=float),
        'circularity': np.array(result_circ, dtype=float),
    }


def enrich_sweep_with_crypt_counts(sweep_data, model_type, max_timesteps=30):
    """Compute crypt count timeseries for every run in a sweep."""
    if _crypt_count_mod is None or model_type not in ('node2d', 'vertex2d'):
        return

    total_runs = sum(len(v) for v in sweep_data.values())
    done = 0

    for pval in sorted(sweep_data.keys()):
        for run_num, data in sweep_data[pval]:
            done += 1
            data_dir = data.get('_data_dir')
            if not data_dir:
                continue

            print(f"\r    Crypt counting: {done}/{total_runs} "
                  f"(param={pval}, run={run_num})", end='', flush=True)

            cc = compute_crypt_count_timeseries(data_dir, model_type, max_timesteps)
            if cc:
                data['num_crypts'] = cc['num_crypts']
                data['circularity'] = cc['circularity']
                data['_cc_time'] = cc['time']

    if total_runs > 0:
        print()


def save_crypt_count_csv(sweep_data, param_name, output_path, model_type=''):
    """Save per-run crypt count timeseries to CSV."""
    rows = []
    for pval in sorted(sweep_data.keys()):
        for run_num, data in sweep_data[pval]:
            cc_time = data.get('_cc_time')
            cc_n = data.get('num_crypts')
            cc_c = data.get('circularity')
            if cc_time is None or cc_n is None:
                continue
            for i in range(len(cc_time)):
                rows.append([model_type, pval, run_num, cc_time[i],
                             cc_n[i], cc_c[i]])

    if not rows:
        return

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['model', param_name, 'replicate', 'time',
                         'num_crypts', 'circularity'])
        writer.writerows(rows)
    print(f"  Saved crypt count CSV: {output_path}")


def plot_timeseries_by_param(sweep_data, column, param_name, param_unit,
                             ylabel, title, output_path, model_type=None,
                             time_col='time'):
    """Plot a time-series column with one line per parameter value (mean +/- SD)."""
    param_vals = sorted(sweep_data.keys())
    colours = get_param_colours(param_vals)

    fig, ax = plt.subplots(figsize=figsize())
    for pval in param_vals:
        runs = sweep_data[pval]
        times, mean, std, n = aggregate_timeseries(runs, column, time_col=time_col)
        if len(times) == 0:
            continue
        colour = colours[pval]
        label = _param_label(pval, param_unit)
        ax.plot(times, mean, color=colour, label=label)

    ax.set_xlabel('Time (hours)')
    ax.set_ylabel(ylabel)
    _add_legend(ax, title=param_name, ncol=2)
    _save_fig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_final_value_boxplot(sweep_data, column, param_name, param_unit,
                             ylabel, title, output_path, model_type=None):
    """Box plot of the final value of a column across replicates, per parameter."""
    param_vals = sorted(sweep_data.keys())

    final_values = []
    labels = []
    for pval in param_vals:
        vals = []
        for _, data in sweep_data[pval]:
            if column in data and len(data[column]) > 0:
                vals.append(data[column][-1])
        final_values.append(vals)
        labels.append(f'{pval}')

    box_colour = 'steelblue'
    fig, ax = plt.subplots(figsize=figsize())
    bp = ax.boxplot(final_values, positions=range(len(param_vals)),
                    widths=0.5, patch_artist=True,
                    medianprops=dict(color='black', linewidth=2))
    for patch in bp['boxes']:
        patch.set_facecolor(box_colour)
        patch.set_alpha(0.4)

    ax.set_xticks(range(len(param_vals)))
    ax.set_xticklabels(labels, rotation=45)
    ax.set_xlabel(f'{param_name} ({param_unit})' if param_unit else param_name)
    ax.set_ylabel(ylabel)
    _save_fig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_mean_vs_param(sweep_data, column, param_name, param_unit,
                       ylabel, title, output_path, model_type=None,
                       use_final=True, logx=False):
    """Mean +/- SD of a column value vs parameter, as errorbar plot."""
    param_vals = sorted(sweep_data.keys())
    means = []
    stds = []
    for pval in param_vals:
        vals = []
        for _, data in sweep_data[pval]:
            if column in data and len(data[column]) > 0:
                vals.append(data[column][-1] if use_final else np.mean(data[column]))
        means.append(np.mean(vals) if vals else np.nan)
        stds.append(np.std(vals) if vals else np.nan)

    fig, ax = plt.subplots(figsize=figsize())
    ax.errorbar(param_vals, means, yerr=stds, fmt='o-', capsize=5,
                markersize=8, linewidth=2, color='steelblue')
    ax.set_xlabel(f'{param_name} ({param_unit})' if param_unit else param_name)
    ax.set_ylabel(ylabel)
    if logx and all(v > 0 for v in param_vals):
        ax.set_xscale('log')
    _save_fig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_multi_model_comparison(all_sweep_data, column, param_name, param_unit,
                                ylabel, title, output_path, use_final=True,
                                logx=False):
    """Overlay multiple model types on one plot (mean +/- SD vs parameter)."""
    model_styles = {
        'node2d':   ('steelblue', 'o', '-'),
        'vertex2d': ('coral',     's', '-'),
        'node3d':   ('seagreen',  '^', '--'),
        'vertex3d': ('orchid',    'D', '--'),
    }

    fig, ax = plt.subplots(figsize=figsize())
    for model, sweep_data in sorted(all_sweep_data.items()):
        colour, marker, ls = model_styles.get(model, ('gray', 'x', ':'))
        param_vals = sorted(sweep_data.keys())
        means = []
        stds = []
        for pval in param_vals:
            vals = []
            for _, data in sweep_data[pval]:
                if column in data and len(data[column]) > 0:
                    vals.append(data[column][-1] if use_final
                                else np.mean(data[column]))
            means.append(np.mean(vals) if vals else np.nan)
            stds.append(np.std(vals) if vals else np.nan)

        ax.errorbar(param_vals, means, yerr=stds, fmt=f'{marker}{ls}',
                    color=colour, capsize=5, markersize=8, linewidth=2,
                    label=model)

    ax.set_xlabel(f'{param_name} ({param_unit})' if param_unit else param_name)
    ax.set_ylabel(ylabel)
    if logx:
        ax.set_xscale('log')
    _add_legend(ax, title='Model', ncol=2)
    _save_fig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def print_summary_table(sweep_data, param_name, model_type=''):
    """Print text summary of sweep results."""
    param_vals = sorted(sweep_data.keys())
    header = f"{model_type}: Summary by {param_name}" if model_type else f"Summary by {param_name}"
    print(f"\n{'=' * 70}")
    print(f"  {header}")
    print(f"{'=' * 70}")
    print(f"  {param_name:>12}  {'N':>4}  {'cells':>8}  {'mean_r':>8}  "
          f"{'var_r':>8}  {'r_range':>8}")
    print(f"  {'-'*12}  {'-'*4}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*8}")

    for pval in param_vals:
        runs = sweep_data[pval]
        n = len(runs)
        final_cells = [d['num_cells'][-1] for _, d in runs
                       if 'num_cells' in d and len(d['num_cells']) > 0]
        final_r = [d['mean_r'][-1] for _, d in runs
                   if 'mean_r' in d and len(d['mean_r']) > 0]
        final_var = [d['var_r'][-1] for _, d in runs
                     if 'var_r' in d and len(d['var_r']) > 0]
        final_range = [d['r_range'][-1] for _, d in runs
                       if 'r_range' in d and len(d['r_range']) > 0]

        fc = f"{np.mean(final_cells):.1f}" if final_cells else "N/A"
        fr = f"{np.mean(final_r):.3f}" if final_r else "N/A"
        fv = f"{np.mean(final_var):.4f}" if final_var else "N/A"
        frng = f"{np.mean(final_range):.3f}" if final_range else "N/A"
        print(f"  {pval:12.3f}  {n:4d}  {fc:>8}  {fr:>8}  {fv:>8}  {frng:>8}")

    print(f"{'=' * 70}\n")


def save_summary_csv(sweep_data, param_name, output_path, model_type=''):
    """Save summary statistics to CSV."""
    param_vals = sorted(sweep_data.keys())
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['model', param_name, 'replicate', 'final_cells',
                         'final_mean_r', 'final_var_r', 'final_r_range',
                         'final_max_vol', 'final_min_vol',
                         'final_num_crypts', 'final_circularity'])
        for pval in param_vals:
            for run_num, data in sweep_data[pval]:
                row = [model_type, pval, run_num]
                for col in ['num_cells', 'mean_r', 'var_r', 'r_range',
                            'max_vol', 'min_vol',
                            'num_crypts', 'circularity']:
                    if col in data and len(data[col]) > 0:
                        row.append(data[col][-1])
                    else:
                        row.append('')
                writer.writerow(row)
    print(f"  Saved CSV: {output_path}")


def generate_standard_plots(sweep_data, param_name, param_unit, model_type,
                            plots_dir, logx=False, crypt_count=True):
    """Generate the standard set of time-series and summary plots."""
    setup_style()
    prefix = model_type

    ts_plots = [
        ('num_cells', 'Number of Cells', 'Cell Proliferation Over Time'),
        ('mean_r', 'Mean Radial Distance (CD)', 'Organoid Growth Over Time'),
        ('var_r', 'Radial Variance (CD\u00b2)', 'Radial Variance Over Time (Budding Indicator)'),
        ('r_range', 'Radial Range (CD)', 'Budding Amplitude Over Time'),
        ('max_vol', 'Max Cell Volume', 'Maximum Cell Volume Over Time'),
        ('min_vol', 'Min Cell Volume', 'Minimum Cell Volume Over Time'),
    ]

    for col, ylabel, title in ts_plots:
        path = os.path.join(plots_dir, f'{prefix}_{col}_timeseries.svg')
        plot_timeseries_by_param(sweep_data, col, param_name, param_unit,
                                 ylabel, title, path, model_type)

    summary_plots = [
        ('num_cells', 'Final Cell Count', 'Final Cell Count vs ' + param_name),
        ('mean_r', 'Final Mean Radius (CD)', 'Final Organoid Size vs ' + param_name),
        ('var_r', 'Final Radial Variance', 'Final Radial Variance vs ' + param_name),
        ('r_range', 'Final Radial Range (CD)', 'Final Budding Amplitude vs ' + param_name),
    ]

    for col, ylabel, title in summary_plots:
        path = os.path.join(plots_dir, f'{prefix}_{col}_boxplot.svg')
        plot_final_value_boxplot(sweep_data, col, param_name, param_unit,
                                ylabel, title, path, model_type)
        path = os.path.join(plots_dir, f'{prefix}_{col}_vs_{param_name}.svg')
        plot_mean_vs_param(sweep_data, col, param_name, param_unit,
                           ylabel, title, path, model_type, logx=logx)

    if crypt_count and _crypt_count_mod and model_type in ('node2d', 'vertex2d'):
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
                path = os.path.join(plots_dir, f'{prefix}_{col}_timeseries.svg')
                plot_timeseries_by_param(sweep_data, col, param_name, param_unit,
                                         ylabel, title, path, model_type,
                                         time_col='_cc_time')

            for col, ylabel, title in [
                ('circularity', 'Final Circularity', 'Final Circularity vs ' + param_name),
            ]:
                path = os.path.join(plots_dir, f'{prefix}_{col}_boxplot.svg')
                plot_final_value_boxplot(sweep_data, col, param_name, param_unit,
                                         ylabel, title, path, model_type)
                path = os.path.join(plots_dir, f'{prefix}_{col}_vs_{param_name}.svg')
                plot_mean_vs_param(sweep_data, col, param_name, param_unit,
                                    ylabel, title, path, model_type, logx=logx)

            for col, ylabel, title in [
                ('num_crypts', 'Number of Crypts', 'Crypt Count Over Time'),
            ]:
                path = os.path.join(plots_dir, f'{prefix}_{col}_timeseries.svg')
                plot_timeseries_by_param(sweep_data, col, param_name, param_unit,
                                         ylabel, title, path, model_type,
                                         time_col='_cc_time')

            for col, ylabel, title in [
                ('num_crypts', 'Final Crypt Count', 'Final Crypt Count vs ' + param_name),
            ]:
                path = os.path.join(plots_dir, f'{prefix}_{col}_boxplot.svg')
                plot_final_value_boxplot(sweep_data, col, param_name, param_unit,
                                         ylabel, title, path, model_type)
                path = os.path.join(plots_dir, f'{prefix}_{col}_vs_{param_name}.svg')
                plot_mean_vs_param(sweep_data, col, param_name, param_unit,
                                    ylabel, title, path, model_type, logx=logx)

            cc_csv = os.path.join(plots_dir, f'{prefix}_crypt_count_timeseries.csv')
            save_crypt_count_csv(sweep_data, param_name, cc_csv, model_type)

    print_summary_table(sweep_data, param_name, model_type)
    csv_path = os.path.join(plots_dir, f'{prefix}_summary.csv')
    save_summary_csv(sweep_data, param_name, csv_path, model_type)
