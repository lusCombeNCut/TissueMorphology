#!/usr/bin/env python3
"""
timestep_summary.py

Summarise the number of timesteps completed for each parameter configuration
in the CryptBudding stiffness sweep.  This detects simulations that stopped
early (e.g. due to numerical instability) vs those that ran to completion.

Reads timestep information from:
  1. results.pvd  — counts <DataSet> entries (most reliable)
  2. VTU files     — counts results_*.vtu files as fallback

Also optionally fixes broken PVD files (missing XML closing tags) in-place.

Usage:
  python timestep_summary.py --base-dir /path/to/merged --model vertex3d
  python timestep_summary.py --base-dir /path/to/merged --model vertex3d --fix-pvd
"""

import os
import sys
import glob
import re
import argparse
import csv
from collections import defaultdict

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


def count_timesteps_pvd(pvd_path):
    """Count the number of timestep entries in a PVD file.
    Returns (n_timesteps, max_timestep_value, is_complete)."""
    if not os.path.exists(pvd_path):
        return 0, 0.0, False

    n = 0
    max_t = 0.0
    has_closing = False
    with open(pvd_path, 'r') as f:
        for line in f:
            m = re.search(r'timestep="([^"]+)"', line)
            if m:
                n += 1
                t = float(m.group(1))
                if t > max_t:
                    max_t = t
            if '</VTKFile>' in line:
                has_closing = True
    return n, max_t, has_closing


def count_vtu_files(data_dir):
    """Count results_*.vtu files (excluding *_tri_* triangle meshes)."""
    vtu = glob.glob(os.path.join(data_dir, 'results_*.vtu'))
    vtu = [f for f in vtu if '_tri_' not in os.path.basename(f)]
    return len(vtu)


def get_max_vtu_timestep(data_dir):
    """Extract the maximum timestep number from VTU filenames."""
    vtu = glob.glob(os.path.join(data_dir, 'results_*.vtu'))
    vtu = [f for f in vtu if '_tri_' not in os.path.basename(f)]
    max_t = 0
    for f in vtu:
        m = re.search(r'results_(\d+)\.vtu$', os.path.basename(f))
        if m:
            t = int(m.group(1))
            if t > max_t:
                max_t = t
    return max_t


def fix_pvd_file(pvd_path):
    """Add missing XML closing tags to a PVD file if needed."""
    if not os.path.exists(pvd_path):
        return False

    with open(pvd_path, 'r') as f:
        content = f.read()

    if '</VTKFile>' in content:
        return False  # already OK

    # Strip trailing whitespace/newlines
    content = content.rstrip()

    # Add closing tags
    if '</Collection>' not in content:
        content += '\n    </Collection>\n'
    if '</VTKFile>' not in content:
        content += '</VTKFile>\n'

    with open(pvd_path, 'w') as f:
        f.write(content)
    return True


def analyse_timesteps(base_dir, model_type, fix_pvd=False):
    """Scan all runs for a model type and collect timestep statistics.

    Returns list of dicts with per-run info.
    """
    model_dir = os.path.join(base_dir, 'CryptBudding', model_type)
    if not os.path.isdir(model_dir):
        print(f"  Directory not found: {model_dir}")
        return []

    records = []
    stiffness_dirs = sorted(glob.glob(os.path.join(model_dir, 'stiffness_*')))

    for s_dir in stiffness_dirs:
        try:
            stiffness = float(os.path.basename(s_dir).replace('stiffness_', ''))
        except ValueError:
            continue

        run_dirs = sorted(glob.glob(os.path.join(s_dir, 'run_*')))
        for r_dir in run_dirs:
            run_num = int(os.path.basename(r_dir).replace('run_', ''))

            # Find all results_from_time_* subdirectories
            subdirs = sorted(glob.glob(os.path.join(r_dir, 'results_from_time_*')))
            if not subdirs:
                subdirs = [r_dir]

            total_timesteps = 0
            max_time = 0.0
            pvd_complete = True
            n_phases = len(subdirs)
            pvd_fixed = False

            for sd in subdirs:
                pvd_path = os.path.join(sd, 'results.pvd')
                n_ts, max_t, complete = count_timesteps_pvd(pvd_path)

                if n_ts == 0:
                    # Fall back to counting VTU files
                    n_ts = count_vtu_files(sd)
                    max_t = float(get_max_vtu_timestep(sd))

                total_timesteps += n_ts
                if max_t > max_time:
                    max_time = max_t
                if not complete:
                    pvd_complete = False

                if fix_pvd and not complete and os.path.exists(pvd_path):
                    if fix_pvd_file(pvd_path):
                        pvd_fixed = True

            records.append({
                'stiffness': stiffness,
                'run': run_num,
                'total_timesteps': total_timesteps,
                'max_time': max_time,
                'n_phases': n_phases,
                'pvd_complete': pvd_complete,
                'pvd_fixed': pvd_fixed,
            })

    return records


def print_summary(records, model_type, expected_max_time=None):
    """Print formatted summary tables."""
    if not records:
        print("  No records found.")
        return

    # Group by stiffness
    by_stiffness = defaultdict(list)
    for r in records:
        by_stiffness[r['stiffness']].append(r)

    # Detect expected max time from the data
    if expected_max_time is None:
        all_max = [r['max_time'] for r in records]
        expected_max_time = max(all_max)

    print(f"\n{'=' * 90}")
    print(f"  Timestep Summary: {model_type}  (expected end time: {expected_max_time})")
    print(f"{'=' * 90}")

    # Per-run detail
    print(f"\n  {'Stiffness':>10}  {'Run':>4}  {'Timesteps':>10}  {'Max Time':>10}  {'Phases':>6}  {'PVD OK':>6}  {'Status':>12}")
    print(f"  {'-'*10}  {'-'*4}  {'-'*10}  {'-'*10}  {'-'*6}  {'-'*6}  {'-'*12}")

    early_stopped = []
    for s in sorted(by_stiffness.keys()):
        runs = sorted(by_stiffness[s], key=lambda x: x['run'])
        for r in runs:
            status = 'COMPLETE' if r['max_time'] >= expected_max_time * 0.95 else 'EARLY STOP'
            pvd_ok = 'Yes' if r['pvd_complete'] else 'NO'
            print(f"  {s:10.1f}  {r['run']:4d}  {r['total_timesteps']:10d}  "
                  f"{r['max_time']:10.1f}  {r['n_phases']:6d}  {pvd_ok:>6}  {status:>12}")
            if status == 'EARLY STOP':
                early_stopped.append(r)

    # Summary statistics per stiffness
    print(f"\n  {'=' * 80}")
    print(f"  Summary Statistics per Stiffness")
    print(f"  {'=' * 80}")
    print(f"  {'Stiffness':>10}  {'N':>4}  {'Mean TS':>10}  {'SD TS':>10}  "
          f"{'Min TS':>8}  {'Max TS':>8}  {'Complete':>8}  {'Early':>6}")
    print(f"  {'-'*10}  {'-'*4}  {'-'*10}  {'-'*10}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*6}")

    for s in sorted(by_stiffness.keys()):
        runs = by_stiffness[s]
        ts_list = [r['total_timesteps'] for r in runs]
        n_complete = sum(1 for r in runs if r['max_time'] >= expected_max_time * 0.95)
        n_early = len(runs) - n_complete

        if HAS_NUMPY:
            mean_ts = np.mean(ts_list)
            sd_ts = np.std(ts_list)
        else:
            mean_ts = sum(ts_list) / len(ts_list)
            var_ts = sum((x - mean_ts) ** 2 for x in ts_list) / len(ts_list)
            sd_ts = var_ts ** 0.5

        print(f"  {s:10.1f}  {len(runs):4d}  {mean_ts:10.1f}  {sd_ts:10.1f}  "
              f"{min(ts_list):8d}  {max(ts_list):8d}  {n_complete:8d}  {n_early:6d}")

    # Overall
    all_ts = [r['total_timesteps'] for r in records]
    total_complete = sum(1 for r in records if r['max_time'] >= expected_max_time * 0.95)
    total_early = len(records) - total_complete
    broken_pvd = sum(1 for r in records if not r['pvd_complete'])

    print(f"\n  Overall: {len(records)} runs, {total_complete} complete, "
          f"{total_early} early-stopped, {broken_pvd} broken PVD files")

    if early_stopped:
        print(f"\n  Early-stopped runs:")
        for r in early_stopped:
            print(f"    stiffness={r['stiffness']:.1f} run={r['run']} "
                  f"stopped at time={r['max_time']:.1f} ({r['total_timesteps']} timesteps)")

    print(f"{'=' * 90}\n")
    return expected_max_time


def save_timestep_csv(records, output_path):
    """Save per-run timestep data to CSV."""
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['stiffness', 'run', 'total_timesteps', 'max_time',
                         'n_phases', 'pvd_complete', 'status'])
        max_time_all = max(r['max_time'] for r in records) if records else 0
        for r in sorted(records, key=lambda x: (x['stiffness'], x['run'])):
            status = 'complete' if r['max_time'] >= max_time_all * 0.95 else 'early_stop'
            writer.writerow([r['stiffness'], r['run'], r['total_timesteps'],
                             r['max_time'], r['n_phases'], r['pvd_complete'], status])
    print(f"  Saved timestep CSV: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Summarise timestep completion per parameter configuration')
    parser.add_argument('--base-dir', type=str, required=True,
                        help='Base output directory (contains CryptBudding/<model>/...)')
    parser.add_argument('--model', nargs='+', default=None,
                        choices=['node2d', 'vertex2d', 'node3d', 'vertex3d'],
                        help='Model type(s) to analyse')
    parser.add_argument('-o', '--output-dir', type=str,
                        default='timestep_analysis_output',
                        help='Output directory for CSV')
    parser.add_argument('--fix-pvd', action='store_true',
                        help='Fix broken PVD files by adding missing XML closing tags')

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    all_models = ['node2d', 'vertex2d', 'node3d', 'vertex3d']
    if args.model:
        models = args.model
    else:
        models = []
        for m in all_models:
            d = os.path.join(args.base_dir, 'CryptBudding', m)
            if os.path.isdir(d):
                models.append(m)

    if not models:
        print(f"No model directories found under {args.base_dir}/CryptBudding/")
        sys.exit(1)

    print(f"Models to analyse: {', '.join(models)}")
    print(f"Base directory:    {args.base_dir}")
    if args.fix_pvd:
        print(f"  ** Will fix broken PVD files in-place **")

    for model_type in models:
        print(f"\n--- Analysing {model_type} ---")
        records = analyse_timesteps(args.base_dir, model_type, fix_pvd=args.fix_pvd)

        if records:
            print_summary(records, model_type)
            csv_path = os.path.join(args.output_dir, f'{model_type}_timestep_summary.csv')
            save_timestep_csv(records, csv_path)

            if args.fix_pvd:
                n_fixed = sum(1 for r in records if r.get('pvd_fixed'))
                if n_fixed:
                    print(f"  Fixed {n_fixed} broken PVD files")

    print("\nDone.")


if __name__ == '__main__':
    main()
