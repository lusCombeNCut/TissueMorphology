#!/usr/bin/env python3
"""
debug_circularity.py — Generate crypt-count debug plots for 2D simulation runs.

For each matched run, loads the organoid boundary at a given simulation time,
runs the full crypt-counting pipeline, and saves a 2×3 debug figure showing:
  1. Input boundary
  2. Fourier-smoothed boundary + circularity
  3. Curvature along boundary
  4. Curvature-scaled normals (inside/outside detection)
  5. All detected crypt sections (before filtering)
  6. Final filtered crypts

Usage:
  # All vertex2d runs at day 4:
  python3 debug_circularity.py --data-dir 2d_sim_04042026/16564370_vertex2d_...

  # One specific stiffness:
  python3 debug_circularity.py --data-dir 2d_sim_04042026 --stiffness 1300 --model vertex2d

  # Different time point:
  python3 debug_circularity.py --data-dir 2d_sim_04042026 --time 48
"""

import os
import sys
import glob
import argparse

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
    plot_debug_analysis,
)

STIFFNESS_JSON_PATH = 'forces.GhostNodeECM.ViscoelasticECM.ecmShearModulusPa'


def main():
    parser = argparse.ArgumentParser(
        description='Generate crypt-count debug plots at a chosen simulation time.')
    parser.add_argument('--data-dir', required=True,
                        help='Root of simulation output (searched recursively)')
    parser.add_argument('--model', default='all',
                        choices=['node2d', 'vertex2d', 'all'],
                        help='Model type to process (default: all)')
    parser.add_argument('--stiffness', type=float, default=None,
                        help='Only process this stiffness value (Pa)')
    parser.add_argument('--run', type=int, default=None,
                        help='Only process this run number')
    parser.add_argument('--time', type=float, default=96.0,
                        help='Simulation time in hours (default: 96 = day 4)')
    parser.add_argument('-o', '--output-dir', default=None,
                        help='Output directory for debug plots (default: <data-dir>/debug_plots)')
    args = parser.parse_args()

    out_dir = args.output_dir or os.path.join(args.data_dir, 'debug_plots')
    os.makedirs(out_dir, exist_ok=True)
    model_filter = None if args.model == 'all' else args.model

    csvs = find_summary_csvs(args.data_dir)
    if not csvs:
        print(f'No crypt_summary.csv files found under {args.data_dir}')
        sys.exit(1)

    n_saved = 0

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

        if args.stiffness is not None and stiffness != args.stiffness:
            continue
        if args.run is not None and run_number != args.run:
            continue

        data_dir = os.path.dirname(csv_path)   # results_from_time_* dir
        tag = f'{model}_G{int(stiffness) if stiffness else "?"}Pa_r{run_number}'
        print(f'Processing {tag} ...')

        try:
            if model == 'node2d':
                boundary, cell_types = load_outline_at_time(data_dir, args.time, dt)
            else:
                boundary, cell_types = load_final_vertex_boundary(
                    data_dir, target_time=args.time, dt=dt)

            result = count_crypts_simple_method(boundary, boundary_is_ordered=True)

            out_path = os.path.join(out_dir, f'debug_{tag}_t{int(args.time)}h.png')
            plot_debug_analysis(
                result,
                cell_types=cell_types,
                output_path=out_path,
                title_prefix=f'{tag}  t={args.time}h',
            )
            print(f'  circularity={result.circularity:.4f}  crypts={result.num_crypts}'
                  f'  -> {os.path.basename(out_path)}')
            n_saved += 1

        except Exception as e:
            print(f'  ERROR: {e}')

    print(f'\nSaved {n_saved} debug plots to {out_dir}')


if __name__ == '__main__':
    main()
