#!/usr/bin/env python3
"""
sweep_3d_and_compare.py
=======================
1. Runs TestViscoelasticShearExperiment3d for G = 100–1900 Pa (200 Pa steps).
2. Reads existing 2D sweep results from shear_sweep_results.csv.
3. Generates TWO output figures:

   (a) shear_sweep_2d_3d_scatter.svg
       Two side-by-side panels (2D left, 3D right), each a measured-vs-target
       scatter plot matching the style of shear_sweep_validation.png.
       Serif font, no overall title, same axis limits.

   (b) shear_sweep_timeseries.svg
       Two side-by-side panels (G=300 Pa left, G=1300 Pa right).
       Each shows shear stress (Pa) vs time for both 2D and 3D.
       Uses the CSV data written by the compare_shear_2d_3d.py run.

Usage
-----
    cd /home/orlando/Thesis
    python3 projects/TissueMorphology/experiments/ViscoelasticRheology/sweep_3d_and_compare.py
"""

import os
import sys
import csv
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Paths ──────────────────────────────────────────────────────────────────
BUILD_DIR   = os.path.expanduser("~/Thesis")
TEST_OUTPUT = os.environ.get("CHASTE_TEST_OUTPUT",
                             os.path.expanduser("~/Thesis/testoutput"))
VC_DIR      = os.path.join(TEST_OUTPUT, "ViscoelasticCompression")
EXE_2D      = os.path.join(BUILD_DIR,
    "projects/TissueMorphology/test/TestViscoelasticShearExperiment")
EXE_3D      = os.path.join(BUILD_DIR,
    "projects/TissueMorphology/test/TestViscoelasticShearExperiment3d")
CSV_2D_SWEEP = os.path.join(VC_DIR, "shear_sweep_results.csv")

OUT_SCATTER   = os.path.join(VC_DIR, "shear_sweep_2d_3d_scatter.svg")
OUT_TIMESERIES = os.path.join(VC_DIR, "shear_timeseries_2d_3d.svg")

G_VALUES   = list(range(100, 2000, 200))   # 100, 300, 500, …, 1900
ARM_RATIO  = 0.2
TAU        = 1.0
TS_G_VALS  = [300, 1300]                    # G values for timeseries panel

# ── Shared plot style (matching shear_sweep_validation.png) ────────────────
RCPARAMS = {
    "font.family":      "serif",
    "mathtext.fontset": "dejavuserif",
    "axes.labelsize":   11,
    "xtick.labelsize":  10,
    "ytick.labelsize":  10,
    "legend.fontsize":  10,
}

# ── Helpers ────────────────────────────────────────────────────────────────

def build(target: str) -> None:
    print(f"Building {target} ...", flush=True)
    r = subprocess.run(["make", "-j4", target], cwd=BUILD_DIR,
                       capture_output=True, text=True)
    if r.returncode != 0:
        print("BUILD FAILED:\n", r.stderr[-2000:])
        sys.exit(1)
    print("  OK")


def run_test(exe: str, g_pa: int, sub_dir: str):
    env = os.environ.copy()
    env["SHEAR_G_PA"]       = str(g_pa)
    env["SHEAR_ARM_RATIO"]  = str(ARM_RATIO)
    env["SHEAR_TAU"]        = str(TAU)
    env["SHEAR_OUTPUT_DIR"] = sub_dir
    env["SHEAR_SKIP_VTP"]   = "1"
    r = subprocess.run([exe], env=env, capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        print(f"  FAILED (G={g_pa}): {r.stderr[-500:]}")
        return None
    for line in r.stdout.splitlines():
        if line.startswith("SHEAR_RESULT,"):
            p = line.split(",")
            return {"G_target": float(p[1]), "G_relax": float(p[2]),
                    "G_peak": float(p[3])}
    return None


def read_csv_cols(path: str, *col_names: str) -> dict:
    """Read named columns from a CSV file, return dict of np.arrays."""
    cols = {c: [] for c in col_names}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            for c in col_names:
                cols[c].append(float(row[c]))
    return {c: np.array(v) for c, v in cols.items()}


def ts_csv(dim: str, g: int) -> str:
    """Path to the timeseries CSV written by compare_shear_2d_3d.py."""
    tag = "2D" if dim == "2D" else "3D"
    return os.path.join(VC_DIR, f"SimpleShear_{tag}_G{g}", "shear_data.csv")


# ── 1. Build ───────────────────────────────────────────────────────────────
build("TestViscoelasticShearExperiment")
build("TestViscoelasticShearExperiment3d")

# ── 2. Run 3D sweep ────────────────────────────────────────────────────────
print("\nRunning 3D shear sweep ...")
results_3d = []
for G in G_VALUES:
    sub = f"SimpleShear_3D_G{G}Pa"
    print(f"  G = {G:5d} Pa  ... ", end="", flush=True)
    res = run_test(EXE_3D, G, sub)
    if res:
        results_3d.append(res)
        ratio = res["G_relax"] / G
        print(f"G_relax = {res['G_relax']:.1f} Pa  ({100*ratio:.1f}%)")
    else:
        print("FAILED")

# ── 3. Load 2D sweep results ───────────────────────────────────────────────
print(f"\nLoading 2D sweep from {CSV_2D_SWEEP} ...")
if not os.path.exists(CSV_2D_SWEEP):
    print("2D sweep CSV not found – running 2D sweep first ...")
    results_2d_raw = []
    for G in G_VALUES:
        sub = f"SimpleShear_G{G}Pa"
        print(f"  G = {G:5d} Pa  ... ", end="", flush=True)
        res = run_test(EXE_2D, G, sub)
        if res:
            results_2d_raw.append(res)
            print(f"G_relax = {res['G_relax']:.1f} Pa")
    data_2d = np.array([[r["G_target"], r["G_relax"], r["G_peak"]]
                        for r in results_2d_raw])
else:
    data_2d = np.genfromtxt(CSV_2D_SWEEP, delimiter=",", skip_header=1)

G_target_2d = data_2d[:, 0]
G_relax_2d  = data_2d[:, 1]
G_peak_2d   = data_2d[:, 2]

data_3d = np.array([[r["G_target"], r["G_relax"], r["G_peak"]]
                    for r in results_3d])
G_target_3d = data_3d[:, 0]
G_relax_3d  = data_3d[:, 1]
G_peak_3d   = data_3d[:, 2]

# ── 4. Make sure timeseries CSVs exist for 300 and 1300 Pa ────────────────
for G in TS_G_VALS:
    for dim, exe in [("2D", EXE_2D), ("3D", EXE_3D)]:
        path = ts_csv(dim, G)
        if not os.path.exists(path):
            sub = f"SimpleShear_{'2D' if dim == '2D' else '3D'}_G{G}"
            print(f"Generating timeseries CSV for {dim} G={G} Pa ...")
            run_test(exe, G, sub)

# ── 5. Figure A: side-by-side scatter plots ────────────────────────────────
print("\nGenerating scatter comparison figure ...")
os.makedirs(VC_DIR, exist_ok=True)

plt.rcParams.update(RCPARAMS)

fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.3), sharey=True)

lim = (0, 2100)
ref = [0, 2000]

panels = [
    ("2D lattice", G_target_2d, G_relax_2d, G_peak_2d),
    ("3D lattice", G_target_3d, G_relax_3d, G_peak_3d),
]

for ax, (title, gt, gr, gp) in zip(axes, panels):
    errors = gr - gt
    rmse   = float(np.sqrt(np.mean(errors**2)))
    ratio  = float(np.mean(gr / gt))

    ax.plot(ref, ref, "k--", lw=1, label="$y = x$")
    ax.scatter(gt, gr, c="tab:blue",   s=55, zorder=3,
               label=r"$G_\mathrm{relax}$")
    ax.scatter(gt, gp, c="tab:orange", s=40, marker="^", zorder=3,
               alpha=0.85, label=r"$G_\mathrm{peak}$")

    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Target $G$ (Pa)")
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left")
    # ax.text(0.97, 0.05,
    #         f"RMSE = {rmse:.0f} Pa\nmean ratio = {ratio:.3f}",
    #         transform=ax.transAxes, ha="right", va="bottom",
    #         fontsize=9, family="serif",
    #         bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.8", lw=0.6))

axes[0].set_ylabel("Measured $G$ (Pa)")

plt.tight_layout()
plt.savefig(OUT_SCATTER, format="svg", bbox_inches="tight")
print(f"  Saved: {OUT_SCATTER}")

# ── 6. Figure B: stress timeseries at 300 and 1300 Pa ─────────────────────
print("Generating timeseries figure ...")

# Line style: colour encodes G value, linestyle encodes 2D/3D
COLORS_G  = {300: "#2176AE", 1300: "#B66D0D"}   # blue / brown
STYLES_DIM = {"2D": "-", "3D": "--"}

fig2, axes2 = plt.subplots(1, 2, figsize=(8, 3.5), sharey=False)

panel_titles = {300: "$G = 300$ Pa (soft)",
                1300: "$G = 1300$ Pa (stiff)"}

for ax, G in zip(axes2, TS_G_VALS):
    ax.set_title(panel_titles[G], fontsize=11)
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Shear stress (Pa)")
    ax.grid(True, alpha=0.3)
    ax.axvline(x=1.0, color="gray", lw=0.7, ls=":", label="End of loading")

    for dim in ["2D", "3D"]:
        path = ts_csv(dim, G)
        if not os.path.exists(path):
            print(f"  WARNING: CSV not found: {path}")
            continue
        d = read_csv_cols(path, "time", "shear_stress_Pa")
        lbl = f"{dim} ({G} Pa)"
        ax.plot(d["time"], d["shear_stress_Pa"],
                color=COLORS_G[G],
                ls=STYLES_DIM[dim],
                lw=1.6,
                label=f"{dim}")

    ax.set_xlim(0, 11)
    # ax.set_ylim(top=0)
    ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig(OUT_TIMESERIES, format="svg", bbox_inches="tight")
print(f"  Saved: {OUT_TIMESERIES}")

# ── 7. Console summary ─────────────────────────────────────────────────────
print("\n" + "=" * 65)
print(f"{'G (Pa)':>8}  {'G_relax 2D':>12}  {'Ratio 2D':>9}  "
      f"{'G_relax 3D':>12}  {'Ratio 3D':>9}")
print("-" * 65)
for r3 in results_3d:
    G  = int(r3["G_target"])
    # find matching 2D row
    idx = np.where(G_target_2d == G)[0]
    if len(idx):
        gr2 = G_relax_2d[idx[0]]
        ra2 = f"{100*gr2/G:.1f}%"
    else:
        gr2, ra2 = float("nan"), "---"
    gr3 = r3["G_relax"]
    ra3 = f"{100*gr3/G:.1f}%"
    print(f"{G:>8}  {gr2:>12.1f}  {ra2:>9}  {gr3:>12.1f}  {ra3:>9}")
print("=" * 65)
