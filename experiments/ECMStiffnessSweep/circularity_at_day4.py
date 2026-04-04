#!/usr/bin/env python3
"""
circularity_at_day4.py — Extract organoid circularity at day 4 of proliferation.

Scans a simulation output directory, finds the VTU/VTP file closest to
t = 96 h (day 4 of proliferation), computes circularity = 4πA/P², and
writes one CSV row per run.

Usage:
  python3 circularity_at_day4.py --data-dir 2d_sim_04042026/16564370_vertex2d_...
  python3 circularity_at_day4.py --data-dir 2d_sim_04042026  --model node2d
  python3 circularity_at_day4.py --data-dir 2d_sim_04042026  -o results.csv
"""

import os
import sys
import argparse
import csv
import glob

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_utils import (
    find_summary_csvs,
    load_params_json,
    extract_param_from_json,
    identify_model_type,
)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from simple_crypt_count import (
    load_outline_at_time,
    load_final_vertex_boundary,
    count_crypts_simple_method,
)

# Day 4 of proliferation in simulation-time hours
TARGET_TIME_H = 96.0

STIFFNESS_JSON_PATH = 'forces.GhostNodeECM.ViscoelasticECM.ecmShearModulusPa'


def circularity_for_run(data_dir, model_type, dt):
    """
    Load the boundary at TARGET_TIME_H and return circularity.

    Args:
        data_dir:   Path to a results_from_time_* directory.
        model_type: 'node2d' or 'vertex2d'.
        dt:         Simulation timestep in hours.

    Returns:
        circularity (float) or None on failure.
    """
    if model_type == 'node2d':
        boundary, _ = load_outline_at_time(data_dir, TARGET_TIME_H, dt)
    elif model_type == 'vertex2d':
        boundary, _ = load_final_vertex_boundary(data_dir,
                                                  target_time=TARGET_TIME_H,
                                                  dt=dt)
    else:
        return None

    result = count_crypts_simple_method(boundary, boundary_is_ordered=True)
    return result.circularity


def main():
    parser = argparse.ArgumentParser(
        description='Compute organoid circularity at day 4 of proliferation.')
    parser.add_argument('--data-dir', required=True,
                        help='Root of simulation output (searched recursively)')
    parser.add_argument('--model', default='all',
                        choices=['node2d', 'vertex2d', 'all'],
                        help='Model type to process (default: all)')
    parser.add_argument('-o', '--output', default=None,
                        help='Output CSV path (default: <data-dir>/circularity_day4.csv)')
    args = parser.parse_args()

    out_path = args.output or os.path.join(args.data_dir, 'circularity_day4.csv')
    model_filter = None if args.model == 'all' else args.model

    csvs = find_summary_csvs(args.data_dir)
    if not csvs:
        print(f'No crypt_summary.csv files found under {args.data_dir}')
        sys.exit(1)

    rows = []

    for csv_path, params_path in sorted(csvs):
        params = load_params_json(params_path) if params_path else None
        model = identify_model_type(csv_path, params)

        if model_filter and model != model_filter:
            continue
        if model not in ('node2d', 'vertex2d'):
            continue

        stiffness = extract_param_from_json(params, STIFFNESS_JSON_PATH) if params else None
        run_number = params.get('simulation', {}).get('runNumber', 0) if params else 0
        dt = params.get('simulation', {}).get('dt', 0.001) if params else 0.001

        data_dir = os.path.dirname(csv_path)   # results_from_time_* dir

        try:
            circ = circularity_for_run(data_dir, model, dt)
            status = 'ok'
        except Exception as e:
            circ = None
            status = str(e)
            print(f'  SKIP {data_dir}: {e}')

        rows.append({
            'model': model,
            'stiffness_Pa': stiffness,
            'run': run_number,
            'circularity': f'{circ:.6f}' if circ is not None else '',
            'status': status,
        })
        if circ is not None:
            print(f'  {model}  G={stiffness} Pa  run={run_number}  '
                  f'circularity={circ:.4f}')

    if not rows:
        print('No matching runs found.')
        sys.exit(1)

    # Sort by model, stiffness, run for clean output
    rows.sort(key=lambda r: (r['model'], r['stiffness_Pa'] or 0, r['run']))

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f,
                                fieldnames=['model', 'stiffness_Pa', 'run',
                                            'circularity', 'status'])
        writer.writeheader()
        writer.writerows(rows)

    print(f'\nWrote {len(rows)} rows to {out_path}')


if __name__ == '__main__':
    main()
