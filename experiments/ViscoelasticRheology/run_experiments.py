#!/usr/bin/env python3
"""
Unified runner for viscoelastic shear experiments.

Subcommands
-----------
  sweep   Run a G sweep (100-1900 Pa) for 2D, 3D, or both.
  compare Run 2D and 3D at specific G values (default: 300, 800, 1300).

Output naming convention (matches plot_results.py expectations):
  Timeseries:    VC_DIR/SimpleShear_{dim}_G{g}/shear_data.csv
  Sweep summary: VC_DIR/shear_sweep_{dim}.csv

Examples
--------
  python run_experiments.py sweep --dim 2d
  python run_experiments.py sweep --dim 3d
  python run_experiments.py sweep --dim both
  python run_experiments.py compare --g 300 800 1300
  python run_experiments.py compare --g 300 800 1300 --dim 2d
"""

from __future__ import annotations

import argparse
import os
import numpy as np

from common import (
    VC_DIR, EXE_2D, EXE_3D,
    build_test, run_shear_test, sweep_csv_path,
    compute_metrics, print_metrics,
)


def _exe(dim: str) -> str:
    return EXE_2D if dim == "2D" else EXE_3D


def _build_target(dim: str) -> str:
    return ("TestViscoelasticShearExperiment" if dim == "2D"
            else "TestViscoelasticShearExperiment3d")


def _run_and_collect(dim: str, g_values: list[int]) -> list[dict]:
    """Run shear tests for each G value, return list of result dicts."""
    exe = _exe(dim)
    results = []
    for G in g_values:
        out_dir = f"SimpleShear_{dim}_G{G}"
        print(f"  {dim}  G = {G:5d} Pa  ... ", end="", flush=True)
        res = run_shear_test(exe, G, out_dir)
        if res:
            ratio = res["G_relax"] / G
            print(f"G_relax = {res['G_relax']:7.1f} Pa  (ratio = {ratio:.3f})")
            results.append(res)
        else:
            print("FAILED")
    return results


def _save_sweep_csv(results: list[dict], path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("G_target_Pa,G_relax_Pa,G_peak_Pa\n")
        for r in results:
            f.write(f"{r['G_target']},{r['G_relax']:.4f},{r['G_peak']:.4f}\n")
    print(f"CSV saved: {path}")


# ── Sweep ─────────────────────────────────────────────────────────────────

def cmd_sweep(args):
    dims = ["2D", "3D"] if args.dim == "both" else [args.dim.upper()]
    g_values = list(range(args.g_min, args.g_max + 1, args.g_step))

    for dim in dims:
        build_test(_build_target(dim))

    for dim in dims:
        print(f"\n--- {dim} sweep ({len(g_values)} points) ---")
        results = _run_and_collect(dim, g_values)
        if not results:
            print(f"No {dim} results collected!")
            continue

        _save_sweep_csv(results, sweep_csv_path(dim))

        data = np.array([[r["G_target"], r["G_relax"]] for r in results])
        m = compute_metrics(data[:, 0], data[:, 1])
        print_metrics(m, len(results), f"{g_values[0]}-{g_values[-1]} Pa")


# ── Compare ───────────────────────────────────────────────────────────────

def cmd_compare(args):
    dims = ["2D", "3D"] if args.dim == "both" else [args.dim.upper()]
    g_values = args.g

    for dim in dims:
        build_test(_build_target(dim))

    all_results = {}
    for dim in dims:
        print(f"\n--- {dim} at G = {g_values} Pa ---")
        all_results[dim] = {
            r["G_target"]: r for r in _run_and_collect(dim, g_values)
        }

    # Summary table
    print(f"\n{'='*65}")
    header = f"{'G_target':>10}"
    for dim in dims:
        header += f"  {'G_relax '+dim:>12}  {'Ratio '+dim:>10}"
    print(header)
    print("-" * 65)
    for G in g_values:
        row = f"{G:>10}"
        for dim in dims:
            r = all_results[dim].get(G)
            if r:
                row += f"  {r['G_relax']:>12.1f}  {100*r['G_relax']/G:>9.1f}%"
            else:
                row += f"  {'---':>12}  {'---':>10}"
        print(row)
    print("=" * 65)


# ── CLI ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Run viscoelastic shear experiments")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("sweep", help="Sweep G across a range")
    p.add_argument("--dim", choices=["2d", "3d", "both"], default="2d")
    p.add_argument("--g-min", type=int, default=100)
    p.add_argument("--g-max", type=int, default=1900)
    p.add_argument("--g-step", type=int, default=200)
    p.set_defaults(func=cmd_sweep)

    p = sub.add_parser("compare", help="Run 2D/3D at specific G values")
    p.add_argument("--g", type=int, nargs="+", default=[300, 800, 1300])
    p.add_argument("--dim", choices=["2d", "3d", "both"], default="both")
    p.set_defaults(func=cmd_compare)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
