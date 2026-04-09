#!/usr/bin/env python3
"""
compare_shear_2d_3d.py
======================
Runs the 2D and 3D simple-shear relaxation tests at three target shear moduli
(300, 800, 1300 Pa), reads the resulting CSV files, and produces a side-by-side
SVG figure suitable for inclusion in the thesis.

Layout
------
  Row 1 (top):    Shear stress vs time — 2D (left) | 3D (right)
  Row 2 (bottom): G_recovered / G_target bar chart  — 2D (left) | 3D (right)

Usage
-----
  cd /home/orlando/Thesis
  python3 projects/TissueMorphology/experiments/ViscoelasticRheology/compare_shear_2d_3d.py
"""

import os
import subprocess
import sys
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── Paths ──────────────────────────────────────────────────────────────────
BUILD_DIR   = "/home/orlando/Thesis"
TEST_OUTPUT = os.environ.get("CHASTE_TEST_OUTPUT",
                             "/home/orlando/Thesis/testoutput")
VISCO_DIR   = os.path.join(TEST_OUTPUT, "ViscoelasticCompression")
OUT_SVG     = os.path.join(VISCO_DIR, "shear_2d_3d_comparison.svg")

# ── Experiment parameters ─────────────────────────────────────────────────
G_VALUES  = [300, 800, 1300]          # Pa
COLORS    = ["#2176AE", "#57B8FF", "#B66D0D"]   # blue, light-blue, brown
LABELS    = ["300 Pa (soft)", "800 Pa", "1300 Pa (stiff)"]

ARM_RATIO = 0.2
TAU       = 1.0

# ── Build helper ───────────────────────────────────────────────────────────

def build_test(test_name: str) -> bool:
    print(f"\n{'='*60}")
    print(f"Building {test_name} ...")
    result = subprocess.run(
        ["make", "-j4", test_name],
        cwd=BUILD_DIR,
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print("BUILD FAILED:\n", result.stderr[-2000:])
        return False
    print(f"  Build OK")
    return True


def run_test(test_name: str, g_pa: float, out_dir: str, skip_vtp: bool = True):
    """Run a test with env-var parameters; return parsed SHEAR_RESULT dict."""
    env = os.environ.copy()
    env["SHEAR_G_PA"]       = str(g_pa)
    env["SHEAR_ARM_RATIO"]  = str(ARM_RATIO)
    env["SHEAR_TAU"]        = str(TAU)
    env["SHEAR_OUTPUT_DIR"] = out_dir
    env["SHEAR_SKIP_VTP"]   = "1" if skip_vtp else "0"

    exe = os.path.join(BUILD_DIR, test_name)
    result = subprocess.run([exe], env=env, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  FAILED (G={g_pa} Pa): return code {result.returncode}")
        print(result.stderr[-1000:])
        return None

    # Parse machine-readable summary line
    for line in result.stdout.splitlines():
        if line.startswith("SHEAR_RESULT,"):
            parts = line.split(",")
            return {
                "G_target":   float(parts[1]),
                "G_relax":    float(parts[2]),
                "G_peak":     float(parts[3]),
                "E0":         float(parts[4]),
                "E1":         float(parts[5]),
                "tau":        float(parts[6]),
                "eta":        float(parts[7]),
                "monotonic":  parts[8].strip() == "1",
            }
    print(f"  WARNING: no SHEAR_RESULT line found for G={g_pa} Pa")
    return None


def csv_path(sub_dir: str) -> str:
    return os.path.join(VISCO_DIR, sub_dir, "shear_data.csv")


def read_csv(fpath: str) -> dict:
    """Return dict of column_name -> numpy array."""
    with open(fpath, newline="") as f:
        reader = csv.DictReader(f)
        cols = {k: [] for k in reader.fieldnames}
        for row in reader:
            for k in reader.fieldnames:
                cols[k].append(float(row[k]))
    return {k: np.array(v) for k, v in cols.items()}


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    os.makedirs(VISCO_DIR, exist_ok=True)

    # ── 1. Build both tests ────────────────────────────────────────────────
    tests = {
        "2D": "projects/TissueMorphology/test/ECMForces/TestViscoelasticShearExperiment",
        "3D": "projects/TissueMorphology/test/ECMForces/TestViscoelasticShearExperiment3d",
    }
    # Use executable name convention: Test<Name> → build target is same
    exe_2d = os.path.join(BUILD_DIR, "projects/TissueMorphology/test/TestViscoelasticShearExperiment")
    exe_3d = os.path.join(BUILD_DIR, "projects/TissueMorphology/test/TestViscoelasticShearExperiment3d")

    if not build_test("TestViscoelasticShearExperiment"):
        sys.exit(1)
    if not build_test("TestViscoelasticShearExperiment3d"):
        sys.exit(1)

    # ── 2. Run experiments ─────────────────────────────────────────────────
    results_2d = {}
    results_3d = {}

    for G in G_VALUES:
        sub_2d = f"SimpleShear_2D_G{G}"
        sub_3d = f"SimpleShear_3D_G{G}"

        print(f"\n--- 2D  G = {G} Pa  (dir: {sub_2d}) ---")
        r2d = run_test(exe_2d, G, sub_2d, skip_vtp=True)
        if r2d:
            results_2d[G] = r2d
            print(f"    G_relax = {r2d['G_relax']:.1f} Pa  "
                  f"({100*r2d['G_relax']/G:.1f}% of target)")

        print(f"\n--- 3D  G = {G} Pa  (dir: {sub_3d}) ---")
        r3d = run_test(exe_3d, G, sub_3d, skip_vtp=True)
        if r3d:
            results_3d[G] = r3d
            print(f"    G_relax = {r3d['G_relax']:.1f} Pa  "
                  f"({100*r3d['G_relax']/G:.1f}% of target)")

    # ── 3. Plot ────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("Simple Shear Stress-Relaxation: 2D vs 3D Ghost Node Network",
                 fontsize=13, fontweight="bold", y=0.98)

    ax_ts_2d, ax_ts_3d = axes[0, 0], axes[0, 1]   # time-series row
    ax_bar_2d, ax_bar_3d = axes[1, 0], axes[1, 1]  # bar-chart row

    # ── Time-series panels ─────────────────────────────────────────────────
    for dim_label, ax, sub_fmt, res_dict in [
        ("2D (8-connected square)",  ax_ts_2d, "SimpleShear_2D_G{}", results_2d),
        ("3D (18-connected cubic)",  ax_ts_3d, "SimpleShear_3D_G{}", results_3d),
    ]:
        ax.set_title(dim_label, fontsize=11)
        ax.set_xlabel("Time (h)")
        ax.set_ylabel("Shear stress (Pa)")

        for G, col, lab in zip(G_VALUES, COLORS, LABELS):
            fpath = csv_path(sub_fmt.format(G))
            if not os.path.exists(fpath):
                print(f"  WARNING: CSV not found: {fpath}")
                continue
            df = read_csv(fpath)
            ax.plot(df["time"], df["shear_stress_Pa"],
                    color=col, linewidth=1.4, label=lab)

        # Mark t_load with a dashed vertical line
        ax.axvline(x=1.0, color="gray", linestyle="--", linewidth=0.8,
                   label="End of loading")
        ax.legend(fontsize=8, loc="upper right")
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(
            lambda v, _: f"{v:.0f}" if v >= 1 else f"{v:.1f}"))
        ax.grid(True, alpha=0.3, linewidth=0.5)

    # ── Recovery bar charts ────────────────────────────────────────────────
    bar_w = 0.35
    x_pos = np.arange(len(G_VALUES))

    for dim_label, ax, res_dict in [
        ("2D",  ax_bar_2d, results_2d),
        ("3D",  ax_bar_3d, results_3d),
    ]:
        ax.set_title(f"G recovery — {dim_label}", fontsize=11)
        ax.set_xlabel("Target G (Pa)")
        ax.set_ylabel("G_relax / G_target")
        ax.set_xticks(x_pos)
        ax.set_xticklabels([str(G) for G in G_VALUES])
        ax.axhline(y=1.0, color="black", linewidth=0.8, linestyle="--")
        ax.set_ylim(0, 1.25)
        ax.grid(True, axis="y", alpha=0.3, linewidth=0.5)

        ratios = []
        bar_cols = []
        for G, col in zip(G_VALUES, COLORS):
            r = res_dict.get(G)
            ratio = (r["G_relax"] / G) if r else 0.0
            ratios.append(ratio)
            bar_cols.append(col)

        bars = ax.bar(x_pos, ratios, width=bar_w * 2, color=bar_cols,
                      edgecolor="k", linewidth=0.6, alpha=0.85)
        for bar, ratio in zip(bars, ratios):
            if ratio > 0:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        ratio + 0.02, f"{100*ratio:.1f}%",
                        ha="center", va="bottom", fontsize=9)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(OUT_SVG, format="svg", dpi=150, bbox_inches="tight")
    print(f"\nSaved: {OUT_SVG}")

    # ── Summary table ──────────────────────────────────────────────────────
    print("\n" + "="*65)
    print(f"{'G_target':>10}  {'G_relax 2D':>12}  {'Ratio 2D':>10}  "
          f"{'G_relax 3D':>12}  {'Ratio 3D':>10}")
    print("-"*65)
    for G in G_VALUES:
        r2 = results_2d.get(G)
        r3 = results_3d.get(G)
        gr2 = f"{r2['G_relax']:.1f}" if r2 else "---"
        gr3 = f"{r3['G_relax']:.1f}" if r3 else "---"
        ra2 = f"{100*r2['G_relax']/G:.1f}%" if r2 else "---"
        ra3 = f"{100*r3['G_relax']/G:.1f}%" if r3 else "---"
        print(f"{G:>10}  {gr2:>12}  {ra2:>10}  {gr3:>12}  {ra3:>10}")
    print("="*65)


if __name__ == "__main__":
    main()
