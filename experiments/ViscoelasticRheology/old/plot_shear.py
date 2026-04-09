"""
Plot shear stress-strain and stress-relaxation curves from the simple shear
test of the viscoelastic ghost node ECM.

Directly compares recovered G to the input shear modulus G=800 Pa.

Usage:
    MPLBACKEND=Agg python3 plot_shear.py [path/to/shear_data.csv]
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# ── Locate CSV ──────────────────────────────────────────────────────────────
if len(sys.argv) > 1:
    csv_path = sys.argv[1]
else:
    chaste_out = os.environ.get("CHASTE_TEST_OUTPUT",
                                os.path.expanduser("~/Thesis/testoutput"))
    csv_path = os.path.join(chaste_out, "ViscoelasticCompression",
                            "SimpleShear_2D", "shear_data.csv")

if not os.path.exists(csv_path):
    print(f"ERROR: CSV not found at {csv_path}")
    sys.exit(1)

data = np.genfromtxt(csv_path, delimiter=",", names=True)

class df:
    time             = data["time"]
    shear_disp       = data["shear_displacement"]
    shear_strain     = data["shear_strain"]
    shear_stress     = data["shear_stress"]
    shear_stress_Pa  = data["shear_stress_Pa"]

# ── Input parameters ─────────────────────────────────────────────────────
G_input = 800.0   # Pa (the target shear modulus)
t_load_end = 1.0
gamma_max = float(df.shear_strain[df.time <= t_load_end][-1])

# ── Split phases ─────────────────────────────────────────────────────────
mask_load  = df.time <= t_load_end
mask_relax = df.time >  t_load_end

# ── Fit exponential to relaxation ────────────────────────────────────────
relax_time       = df.time[mask_relax] - t_load_end
relax_stress_Pa  = df.shear_stress_Pa[mask_relax]

def sls_relax(t, tau_inf, delta_tau, tau_time):
    return tau_inf + delta_tau * np.exp(-t / tau_time)

try:
    p0 = [relax_stress_Pa[-1], relax_stress_Pa[0] - relax_stress_Pa[-1], 1.0]
    popt, _ = curve_fit(sls_relax, relax_time, relax_stress_Pa, p0=p0)
    tau_inf_fit, delta_tau_fit, tau_fit = popt
    fit_ok = True
except Exception as e:
    fit_ok = False
    print(f"Fit failed: {e}")

# ── Compute recovered G ──────────────────────────────────────────────────
peak_stress_Pa  = float(df.shear_stress_Pa[mask_load].max())
final_stress_Pa = float(relax_stress_Pa[-1])

G_peak    = peak_stress_Pa / gamma_max
G_relaxed = final_stress_Pa / gamma_max
if fit_ok:
    G_relaxed_fit = tau_inf_fit / gamma_max

# ── Figure ───────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
fig.suptitle("Simple Shear Test — Viscoelastic Ghost Node ECM\n"
             f"Input G = {G_input:.0f} Pa  |  E₀ = 464 sim  |  E₁ = 23.2 sim  |  τ = 1 h",
             fontsize=11)

# Panel A: Shear stress–strain (loading)
ax = axes[0]
load_strain     = df.shear_strain[mask_load]
load_stress_Pa  = df.shear_stress_Pa[mask_load]
ax.plot(load_strain * 100, load_stress_Pa, "C0", lw=2, label="Simulation")

# Overlay theoretical lines for G_input and G_measured
gamma_line = np.linspace(0, float(load_strain.max()), 100)
ax.plot(gamma_line * 100, G_input * gamma_line, "C3--", lw=1.5,
        label=f"G = {G_input:.0f} Pa (input)")
ax.plot(gamma_line * 100, G_relaxed * gamma_line, "C2--", lw=1.5,
        label=f"G = {G_relaxed:.0f} Pa (measured, relaxed)")

ax.set_xlabel("Shear strain  γ  (%)")
ax.set_ylabel("Shear stress  τ  (Pa)")
ax.set_title("(A) Shear Stress–Strain (loading)")
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)

# Panel B: Shear stress vs time (full protocol)
ax = axes[1]
ax.plot(df.time, df.shear_stress_Pa, "C0", lw=1.5, label="Simulation")
ax.axvline(t_load_end, color="k", ls="--", lw=0.8, label="End loading")
ax.axhline(G_input * gamma_max, color="C3", ls=":", lw=1,
           label=f"G_input × γ = {G_input * gamma_max:.1f} Pa")
ax.set_xlabel("Time  (h)")
ax.set_ylabel("Shear stress  τ  (Pa)")
ax.set_title("(B) Shear Stress vs Time")
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)

# Panel C: Relaxation with fit
ax = axes[2]
ax.plot(relax_time, relax_stress_Pa, "C0", lw=1.5, label="Simulation")
if fit_ok:
    t_fine = np.linspace(0, float(relax_time[-1]), 500)
    ax.plot(t_fine, sls_relax(t_fine, *popt), "C1--", lw=1.5, label="SLS fit")
    ax.axhline(tau_inf_fit, color="grey", ls=":", lw=1,
               label=f"τ_∞ = {tau_inf_fit:.1f} Pa")
    fit_label = (f"Fit: τ(t) = {tau_inf_fit:.1f} + {delta_tau_fit:.1f}·exp(-t/{tau_fit:.3f})\n"
                 f"τ_fit = {tau_fit:.3f} h  (input τ = 1.000 h)\n"
                 f"G_relax(fit) = {G_relaxed_fit:.0f} Pa")
    ax.text(0.97, 0.97, fit_label, transform=ax.transAxes,
            fontsize=7, ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
ax.set_xlabel("Time after hold  (h)")
ax.set_ylabel("Shear stress  τ  (Pa)")
ax.set_title("(C) Shear Relaxation with SLS Fit")
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)

plt.tight_layout()
out_png = os.path.join(os.path.dirname(csv_path), "shear_relaxation.png")
plt.savefig(out_png, dpi=150, bbox_inches="tight")
print(f"Saved: {out_png}")

# ── Summary ──────────────────────────────────────────────────────────────
print("\n=== Shear Modulus Recovery Summary ===")
print(f"Input G:             {G_input:.0f} Pa")
print(f"Max shear strain:    {gamma_max:.4f} ({gamma_max*100:.1f}%)")
print(f"Peak shear stress:   {peak_stress_Pa:.1f} Pa")
print(f"Final shear stress:  {final_stress_Pa:.1f} Pa")
print(f"G_peak  (σ_peak/γ):  {G_peak:.0f} Pa   (ratio to input: {G_peak/G_input:.1f}x)")
print(f"G_relax (σ_final/γ): {G_relaxed:.0f} Pa   (ratio to input: {G_relaxed/G_input:.1f}x)")
if fit_ok:
    print(f"G_relax (from fit):  {G_relaxed_fit:.0f} Pa   (ratio to input: {G_relaxed_fit/G_input:.1f}x)")
    print(f"Fitted τ:            {tau_fit:.4f} h  (input: 1.0000 h)")
print(f"\nIf G_relax ≈ G_input → conversion is correct (G maps to relaxed modulus)")
print(f"If G_relax ≈ 21×G_input → the /SLS_ARM_RATIO in the formula needs review")
