#!/usr/bin/env python3
"""
Plot shear modulus sweep validation from a CSV produced by sweep_shear_modulus.py.

Usage:
    python3 plot_shear_sweep.py [--csv PATH] [--out PATH]

Defaults:
    --csv  $CHASTE_TEST_OUTPUT/ViscoelasticCompression/shear_sweep_results.csv
    --out  $CHASTE_TEST_OUTPUT/ViscoelasticCompression/shear_sweep_validation.png
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Defaults ──────────────────────────────────────────────────────────
OUT_BASE = os.environ.get("CHASTE_TEST_OUTPUT",
                          os.path.expanduser("~/Thesis/testoutput"))
_vc_dir  = os.path.join(OUT_BASE, "ViscoelasticCompression")

csv_path = os.path.join(_vc_dir, "shear_sweep_results.csv")
fig_path = os.path.join(_vc_dir, "shear_sweep_validation.svg")

# ── Parse CLI args ────────────────────────────────────────────────────
args = sys.argv[1:]
i = 0
while i < len(args):
    if args[i] == "--csv" and i + 1 < len(args):
        csv_path = args[i + 1]; i += 2
    elif args[i] == "--out" and i + 1 < len(args):
        fig_path = args[i + 1]; i += 2
    else:
        i += 1

# ── Load data ─────────────────────────────────────────────────────────
if not os.path.exists(csv_path):
    print(f"ERROR: CSV not found at {csv_path}")
    sys.exit(1)

data     = np.genfromtxt(csv_path, delimiter=",", skip_header=1)
G_target = data[:, 0]
G_relax  = data[:, 1]
G_peak   = data[:, 2]

errors     = G_relax - G_target
rmse       = float(np.sqrt(np.mean(errors**2)))
mae        = float(np.mean(np.abs(errors)))
mean_ratio = float(np.mean(G_relax / G_target))

print(f"Loaded {len(G_target)} points from {csv_path}")
print(f"  RMSE = {rmse:.1f} Pa  |  MAE = {mae:.1f} Pa  |  mean ratio = {mean_ratio:.4f}")

# ── Plot ──────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":      "serif",
    "mathtext.fontset": "dejavuserif",
    "axes.labelsize":   11,
    "xtick.labelsize":  10,
    "ytick.labelsize":  10,
    "legend.fontsize":  10,
})

fig, ax = plt.subplots(figsize=(5, 5))

ax.plot([0, 2000], [0, 2000], "k--", lw=1, label="$y = x$")
ax.scatter(G_target, G_relax, c="tab:blue",   s=60, zorder=3,
           label=r"$G_\mathrm{relax}$")
ax.scatter(G_target, G_peak,  c="tab:orange", s=40, marker="^", zorder=3,
           alpha=0.8, label=r"$G_\mathrm{peak}$")
ax.set_xlabel("Target $G$ (Pa)")
ax.set_ylabel("Measured $G$ (Pa)")
ax.legend()
ax.set_xlim(0, 2100)
ax.set_ylim(0, 2100)
ax.set_aspect("equal")
ax.grid(True, alpha=0.3)

plt.tight_layout()
os.makedirs(os.path.dirname(fig_path), exist_ok=True)
plt.savefig(fig_path)
print(f"Plot saved: {fig_path}")
