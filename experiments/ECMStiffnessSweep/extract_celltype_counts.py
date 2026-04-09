#!/usr/bin/env python3
"""
Extract cell type counts from all VTU files across all runs for both models.

Reads cell_type_id (0=Stem, 1=TA, 2=Differentiated) from:
  - vertex2d: CellData (51 VTU files per run, full timeseries)
  - node2d:   PointData (results_0.vtu only, t=0)

Timestep→hours conversion: time_hours = pvd_timestep * 0.001  (dt = 0.001 h)

Output: celltype_counts_all_runs.csv
Columns: model, stiffness_Pa, replicate, time_hours, stem, ta, differentiated, total
"""

import os
import re
import csv
import sys
import xml.etree.ElementTree as ET

try:
    import vtk
except ImportError:
    print('ERROR: vtk not found. Install with: pip3 install vtk')
    sys.exit(1)

# ── Config ─────────────────────────────────────────────────────────────────────
NODE2D_DIR   = '/home/orlando/Thesis/sim_results/2d_sim_04042026/16564369_node2d_2026-04-04_12-50-05'
VERTEX2D_DIR = '/home/orlando/Thesis/sim_results/2d_sim_04042026/16564370_vertex2d_2026-04-04_12-50-05'

OUT_DIR = '/home/orlando/Thesis/sim_results/2d_sim_04042026/plots_vertex2d'
CSV_OUT = os.path.join(OUT_DIR, 'celltype_counts_all_runs.csv')

# PVD timestep units → hours
DT_HOURS = 0.001

# cell_type_id encoding
STEM  = 0
TA    = 1
DIFF  = 2

RUN_RE = re.compile(r'^g(\d+)_r(\d+)$')

# ── VTU reading ────────────────────────────────────────────────────────────────

def read_pvd(pvd_path):
    """Return list of (timestep_int, vtu_filename) from a PVD file."""
    tree = ET.parse(pvd_path)
    root = tree.getroot()
    entries = []
    for ds in root.iter('DataSet'):
        ts   = int(float(ds.attrib['timestep']))
        file = ds.attrib['file']
        entries.append((ts, file))
    return entries


def count_cell_types_from_vtu(vtu_path):
    """
    Read a VTU file and return (stem, ta, diff) counts.
    Checks CellData first (vertex2d), then PointData (node2d).
    Returns (0, 0, 0) if cell_type_id array not found.
    """
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_path)
    reader.Update()
    ug = reader.GetOutput()

    arr = ug.GetCellData().GetArray('cell_type_id')
    if arr is None:
        arr = ug.GetPointData().GetArray('cell_type_id')
    if arr is None:
        return (0, 0, 0)

    stem, ta, diff = 0, 0, 0
    for i in range(arr.GetNumberOfTuples()):
        v = int(arr.GetValue(i))
        if   v == STEM: stem += 1
        elif v == TA:   ta   += 1
        elif v == DIFF: diff += 1

    return (stem, ta, diff)


# ── Run discovery ──────────────────────────────────────────────────────────────

def iter_runs(base_dir, model_name):
    """
    Yield (stiffness_Pa, replicate_idx, results_dir) for every valid run.
    """
    for run_name in sorted(os.listdir(base_dir)):
        m = RUN_RE.match(run_name)
        if not m:
            continue
        stiffness, rep = int(m.group(1)), int(m.group(2))
        cb_dir = os.path.join(base_dir, run_name, 'CryptBudding')
        if not os.path.isdir(cb_dir):
            continue
        hash_dirs = [d for d in os.listdir(cb_dir)
                     if os.path.isdir(os.path.join(cb_dir, d))]
        if not hash_dirs:
            continue
        results_dir = os.path.join(cb_dir, hash_dirs[0], 'results_from_time_0')
        if not os.path.isdir(results_dir):
            continue
        yield stiffness, rep, results_dir


# ── Main extraction ────────────────────────────────────────────────────────────

def extract_all(writer):
    for model_name, base_dir in [('vertex2d', VERTEX2D_DIR),
                                  ('node2d',   NODE2D_DIR)]:
        runs = list(iter_runs(base_dir, model_name))
        print(f'\n{model_name}: {len(runs)} runs found')

        for run_idx, (stiffness, rep, results_dir) in enumerate(runs):
            pvd_path = os.path.join(results_dir, 'results.pvd')
            if not os.path.exists(pvd_path):
                print(f'  SKIP (no PVD): {results_dir}')
                continue

            entries = read_pvd(pvd_path)
            n_files = len(entries)

            print(f'  [{run_idx+1:3d}/{len(runs)}] g{stiffness}_r{rep}  '
                  f'({n_files} VTU files)', end='', flush=True)

            for ts, vtu_name in entries:
                vtu_path = os.path.join(results_dir, vtu_name)
                if not os.path.exists(vtu_path):
                    continue
                stem, ta, diff = count_cell_types_from_vtu(vtu_path)
                writer.writerow({
                    'model':         model_name,
                    'stiffness_Pa':  stiffness,
                    'replicate':     rep,
                    'time_hours':    round(ts * DT_HOURS, 6),
                    'stem':          stem,
                    'ta':            ta,
                    'differentiated': diff,
                    'total':         stem + ta + diff,
                })

            print('  OK')


if __name__ == '__main__':
    os.makedirs(OUT_DIR, exist_ok=True)

    fieldnames = ['model', 'stiffness_Pa', 'replicate', 'time_hours',
                  'stem', 'ta', 'differentiated', 'total']

    with open(CSV_OUT, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        extract_all(writer)

    print(f'\nSaved: {CSV_OUT}')

    # quick sanity check
    import subprocess
    result = subprocess.run(['wc', '-l', CSV_OUT], capture_output=True, text=True)
    print(f'Rows (including header): {result.stdout.strip()}')
