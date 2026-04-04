#!/usr/bin/env python3
"""
Re-run crypt counting analysis with fixed code on the ECM stiffness sweep data.

Fixes applied:
  1. Direct Fourier reconstruction (replaces broken chain-code method)
  2. Simulation-calibrated parameters (replaces paper parameters)
  3. Directed-edge boundary extraction for vertex2d (fixes self-intersections)
  4. vertex2d uses t=35h to avoid late-stage folding artefacts

Handles the flat s{stiffness}_r{run} directory structure.

Usage:
  python3 reanalyse_crypt_counts.py [--debug-plots]
"""

import os
import sys
import glob
import csv
import json
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from simple_crypt_count import (
    count_crypts_simple_method, SIMULATION_PARAMS,
    load_final_outline, load_final_vertex_boundary,
    load_outline_at_time, plot_crypt_outline
)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ===========================================================================
# Configuration
# ===========================================================================
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    'results-download', '03-31_08_56_51-merged')

MODEL_DIRS = {
    'node2d':   os.path.join(BASE, '16535750_node2d_2026-03-31_08-56-51'),
    'vertex2d': os.path.join(BASE, '16535751_vertex2d_2026-03-31_08-56-51'),
}

ANALYSIS_DIR = os.path.join(BASE, 'analysis_ecm_stiffness_2026-03-31_08-56-51')
OUTPUT_DIR = os.path.join(ANALYSIS_DIR, 'crypt_analysis_output')

STIFFNESSES = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 35.0, 50.0, 70.0, 100.0]
N_RUNS = 10

# vertex2d target time (hours) — avoids late-stage folding
VERTEX2D_TARGET_TIME = 35.0
VERTEX2D_DT = 0.001


def find_results_dir(model_dir, stiffness, run):
    """Find the results_from_time_0 directory for a given simulation."""
    sim_dir = os.path.join(model_dir, f's{stiffness}_r{run}')
    results_dir = os.path.join(sim_dir, 'results_from_time_0')
    if os.path.isdir(results_dir):
        return results_dir
    if os.path.isdir(sim_dir):
        return sim_dir
    return None


def analyse_simulation(model_type, results_dir, params):
    """Run crypt counting on a single simulation."""
    if model_type == 'node2d':
        boundary, cell_types = load_final_outline(results_dir)
        result = count_crypts_simple_method(boundary, params, boundary_is_ordered=True)
    elif model_type == 'vertex2d':
        boundary, cell_types = load_final_vertex_boundary(
            results_dir, target_time=VERTEX2D_TARGET_TIME, dt=VERTEX2D_DT)
        result = count_crypts_simple_method(boundary, params, boundary_is_ordered=True)
    else:
        raise ValueError(f"Unsupported model type: {model_type}")
    return result


def run_analysis(models=None, debug_plots=False):
    """Run crypt counting analysis for all 2D models."""
    if models is None:
        models = ['node2d', 'vertex2d']

    params = SIMULATION_PARAMS.copy()
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    if debug_plots:
        debug_dir = os.path.join(OUTPUT_DIR, 'debug_plots')
        os.makedirs(debug_dir, exist_ok=True)

    all_results = []

    for model_type in models:
        model_dir = MODEL_DIRS.get(model_type)
        if not model_dir or not os.path.isdir(model_dir):
            print(f"  SKIP {model_type}: directory not found")
            continue

        print(f"\n--- Analysing {model_type} ---")
        if model_type == 'vertex2d':
            print(f"  (using t={VERTEX2D_TARGET_TIME}h to avoid late-stage folding)")

        for stiffness in STIFFNESSES:
            for run in range(N_RUNS):
                results_dir = find_results_dir(model_dir, stiffness, run)
                if results_dir is None:
                    print(f"  SKIP {model_type}/s{stiffness}_r{run}: not found")
                    continue

                try:
                    result = analyse_simulation(model_type, results_dir, params)
                    all_results.append({
                        'model': model_type,
                        'stiffness': stiffness,
                        'replicate': run,
                        'num_crypts': result.num_crypts,
                        'circularity': result.circularity,
                        'sphericity': '',
                    })
                    print(f"  {model_type}/s{stiffness}_r{run}: "
                          f"{result.num_crypts} crypts, circ={result.circularity:.4f}")

                    if debug_plots:
                        debug_path = os.path.join(
                            debug_dir,
                            f'{model_type}_stiffness_{stiffness}_run_{run}_crypt_outline.png')
                        title = f'{model_type} - ECM Stiffness = {stiffness} - run_{run}'
                        if model_type == 'vertex2d':
                            title += f' (t={VERTEX2D_TARGET_TIME}h)'
                        plot_crypt_outline(result, output_path=debug_path, title=title)

                except Exception as e:
                    print(f"  ERROR {model_type}/s{stiffness}_r{run}: {e}")
                    all_results.append({
                        'model': model_type,
                        'stiffness': stiffness,
                        'replicate': run,
                        'num_crypts': 0,
                        'circularity': float('nan'),
                        'sphericity': '',
                    })

    # Save CSV
    csv_path = os.path.join(OUTPUT_DIR, 'all_crypt_counts.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['model', 'stiffness', 'replicate',
                                                'num_crypts', 'circularity', 'sphericity'])
        writer.writeheader()
        writer.writerows(all_results)
    print(f"\nSaved: {csv_path}")

    return all_results


def main():
    parser = argparse.ArgumentParser(description='Re-run crypt counting with fixed code')
    parser.add_argument('--debug-plots', action='store_true',
                        help='Generate debug plots for each simulation')
    parser.add_argument('--models', nargs='+', default=['node2d', 'vertex2d'],
                        help='Model types to analyse')
    args = parser.parse_args()

    print("=" * 60)
    print("Re-analysing crypt counts with fixed code")
    print(f"Parameters: {SIMULATION_PARAMS}")
    print("=" * 60)

    run_analysis(models=args.models, debug_plots=args.debug_plots)

    print("\nDone!")


if __name__ == '__main__':
    main()
