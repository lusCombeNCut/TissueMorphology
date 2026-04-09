#!/usr/bin/env python3
"""
Sweep shear modulus validation across G = 100-1900 Pa.

Runs the TestViscoelasticShearExperiment for each target G, collects the
measured G_relax (equilibrium shear modulus), and produces:
  1. A scatter plot: measured G vs target G with y=x reference line
  2. MSE and RMSE between target and measured
  3. CSV summary at testoutput/ViscoelasticCompression/shear_sweep_results.csv

Usage:
    python3 sweep_shear_modulus.py [--build-dir /home/orlando/Thesis]
"""

import os
import sys
import subprocess
import re
import numpy as np

# ── Configuration ──────────────────────────────────────────────────
G_VALUES = list(range(100, 2000, 200))  # 100, 300, 500, ..., 1900
ARM_RATIO = 0.2
TAU = 1.0

BUILD_DIR = os.environ.get("CHASTE_BUILD_DIR",
                           os.path.expanduser("~/Thesis"))
TEST_EXE = os.path.join(BUILD_DIR,
    "projects/TissueMorphology/test/TestViscoelasticShearExperiment")
OUT_BASE = os.environ.get("CHASTE_TEST_OUTPUT",
                          os.path.expanduser("~/Thesis/testoutput"))

# ── Parse CLI args ─────────────────────────────────────────────────
for i, arg in enumerate(sys.argv[1:], 1):
    if arg == "--build-dir" and i < len(sys.argv) - 1:
        BUILD_DIR = sys.argv[i + 1]
        TEST_EXE = os.path.join(BUILD_DIR,
            "projects/TissueMorphology/test/TestViscoelasticShearExperiment")

# ── Build first ────────────────────────────────────────────────────
print(f"Building TestViscoelasticShearExperiment...")
ret = subprocess.run(["make", "-j4", "TestViscoelasticShearExperiment"],
                     cwd=BUILD_DIR, capture_output=True, text=True)
if ret.returncode != 0:
    print(f"Build failed:\n{ret.stderr}")
    sys.exit(1)
print("Build OK.\n")

# ── Run sweep ──────────────────────────────────────────────────────
results = []  # (G_target, G_relax_measured, G_peak_measured)

for G in G_VALUES:
    out_dir = f"SimpleShear_G{G}Pa"
    env = os.environ.copy()
    env["SHEAR_G_PA"]       = str(G)
    env["SHEAR_ARM_RATIO"]  = str(ARM_RATIO)
    env["SHEAR_TAU"]        = str(TAU)
    env["SHEAR_OUTPUT_DIR"] = out_dir
    env["SHEAR_SKIP_VTP"]   = "1"

    print(f"Running G = {G:5d} Pa ... ", end="", flush=True)

    proc = subprocess.run([TEST_EXE], env=env, capture_output=True, text=True,
                          timeout=600)

    if proc.returncode != 0:
        print(f"FAILED (exit {proc.returncode})")
        # Try to extract what we can
        print(proc.stderr[-500:] if proc.stderr else "no stderr")
        continue

    # Parse SHEAR_RESULT line
    for line in proc.stdout.split("\n"):
        if line.startswith("SHEAR_RESULT,"):
            parts = line.strip().split(",")
            G_target  = float(parts[1])
            G_relax   = float(parts[2])
            G_peak    = float(parts[3])
            results.append((G_target, G_relax, G_peak))
            ratio = G_relax / G_target if G_target > 0 else 0
            print(f"G_relax = {G_relax:7.1f} Pa  (ratio = {ratio:.3f})")
            break
    else:
        print("FAILED (no SHEAR_RESULT line)")
        print(proc.stdout[-500:])

if len(results) == 0:
    print("No results collected!")
    sys.exit(1)

# ── Compute error metrics ──────────────────────────────────────────
data = np.array(results)
G_target  = data[:, 0]
G_relax   = data[:, 1]
G_peak    = data[:, 2]

errors = G_relax - G_target
mse  = np.mean(errors**2)
rmse = np.sqrt(mse)
mae  = np.mean(np.abs(errors))
mean_ratio = np.mean(G_relax / G_target)

print(f"\n{'='*60}")
print(f"SWEEP SUMMARY ({len(results)} points, G = {G_VALUES[0]}-{G_VALUES[-1]} Pa)")
print(f"{'='*60}")
print(f"  MSE  = {mse:.1f} Pa²")
print(f"  RMSE = {rmse:.1f} Pa")
print(f"  MAE  = {mae:.1f} Pa")
print(f"  Mean G_relax/G_target = {mean_ratio:.4f}")
print(f"  Mean error = {np.mean(errors):.1f} Pa ({np.mean(errors/G_target)*100:.1f}%)")
print()

# ── Save CSV ───────────────────────────────────────────────────────
csv_path = os.path.join(OUT_BASE, "ViscoelasticCompression",
                        "shear_sweep_results.csv")
os.makedirs(os.path.dirname(csv_path), exist_ok=True)
with open(csv_path, "w") as f:
    f.write("G_target_Pa,G_relax_Pa,G_peak_Pa,error_Pa,ratio\n")
    for gt, gr, gp in results:
        f.write(f"{gt},{gr:.4f},{gp:.4f},{gr-gt:.4f},{gr/gt:.6f}\n")
print(f"CSV saved: {csv_path}")

# ── Plot ───────────────────────────────────────────────────────────
plot_script = os.path.join(os.path.dirname(__file__), "plot_shear_sweep.py")
subprocess.run([sys.executable, plot_script, "--csv", csv_path], check=True)
