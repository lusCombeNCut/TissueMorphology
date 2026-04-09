"""
Plot stress-strain and stress-relaxation curves from the uniaxial compression
test of the viscoelastic ghost node ECM.

Usage:
    python plot_compression.py [path/to/compression_data.csv]

Default CSV location: $CHASTE_TEST_OUTPUT/ViscoelasticCompression/compression_data.csv
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
                            "StressRelaxation_2D", "compression_data.csv")

if not os.path.exists(csv_path):
    print(f"ERROR: CSV not found at {csv_path}")
    sys.exit(1)

data = np.genfromtxt(csv_path, delimiter=",", names=True)
# Columns: time, displacement, strain, stress, total_force_y

class df:  # simple namespace
    time         = data["time"]
    displacement = data["displacement"]
    strain       = data["strain"]
    stress       = data["stress"]

# ── Physical conversion ─────────────────────────────────────────────────────
# kF = 0.001 N/m per sim unit  (10 nN traction, 10 µm/h, 10 µm cell diameter)
kF = 0.001     # N/m per sim unit
spacing = 1.0  # sim units = 10 µm
width   = 19.0 # sim units

# Convert stress to Pa: stress [sim] * kF / (spacing [m] * 1 [depth, 2D])
# spacing_m = spacing * 10e-6 = 10e-6 m
spacing_m = spacing * 10e-6
stress_Pa = df.stress * kF / spacing_m

# ── Split phases ─────────────────────────────────────────────────────────────
t_load_end = 1.0
mask_load  = df.time <= t_load_end
mask_relax = df.time >  t_load_end

# ── Fit exponential to relaxation phase ─────────────────────────────────────
relax_time     = df.time[mask_relax] - t_load_end
relax_stress_Pa = stress_Pa[mask_relax]

def sls_relax(t, sigma_inf, delta_sigma, tau):
    return sigma_inf + delta_sigma * np.exp(-t / tau)

try:
    p0 = [relax_stress_Pa[-1], relax_stress_Pa[0] - relax_stress_Pa[-1], 1.0]
    popt, _ = curve_fit(sls_relax, relax_time, relax_stress_Pa, p0=p0)
    sigma_inf_fit, delta_sigma_fit, tau_fit = popt
    fit_ok = True
    fit_label = (f"Fit: σ(t) = {sigma_inf_fit:.0f} + {delta_sigma_fit:.0f}·exp(-t/{tau_fit:.3f})\n"
                 f"τ_fit = {tau_fit:.3f} h  (input τ = 1.000 h)")
except Exception as e:
    fit_ok = False
    fit_label = f"Fit failed: {e}"

# ── Figure ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
fig.suptitle("Viscoelastic Ghost Node ECM — Uniaxial Confined Compression\n"
             "G = 800 Pa  |  E₀ = 464 sim  |  E₁ = 23.2 sim  |  τ = 1 h",
             fontsize=11)

# Panel A: Stress–strain curve (loading only)
ax = axes[0]
load_strain    = df.strain[mask_load]
load_stress_Pa = stress_Pa[mask_load]
ax.plot(load_strain * 100, load_stress_Pa, "C0", lw=2)
ax.set_xlabel("Engineering strain  ε  (%)")
ax.set_ylabel("Compressive stress  σ  (Pa)")
ax.set_title("(A) Stress–Strain (loading)")
ax.grid(True, alpha=0.3)

# Panel B: Stress vs time (full protocol)
ax = axes[1]
ax.plot(df.time, stress_Pa, "C0", lw=1.5, label="Simulation")
ax.axvline(t_load_end, color="k", ls="--", lw=0.8, label="End loading")
ax.set_xlabel("Time  (h)")
ax.set_ylabel("Compressive stress  σ  (Pa)")
ax.set_title("(B) Stress vs Time (full protocol)")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Panel C: Relaxation phase with exponential fit
ax = axes[2]
ax.plot(relax_time, relax_stress_Pa, "C0", lw=1.5, label="Simulation")
if fit_ok:
    t_fine = np.linspace(0, relax_time[-1], 500)
    ax.plot(t_fine, sls_relax(t_fine, *popt), "C1--", lw=1.5, label="SLS fit")
    ax.axhline(sigma_inf_fit, color="grey", ls=":", lw=1, label=f"σ_∞ = {sigma_inf_fit:.0f} Pa")
ax.set_xlabel("Time after hold  (h)")
ax.set_ylabel("Compressive stress  σ  (Pa)")
ax.set_title("(C) Stress Relaxation with SLS Fit")
ax.legend(fontsize=7)
ax.text(0.97, 0.97, fit_label, transform=ax.transAxes,
        fontsize=7, ha="right", va="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
ax.grid(True, alpha=0.3)

plt.tight_layout()
out_png = os.path.join(os.path.dirname(csv_path), "stress_relaxation.png")
plt.savefig(out_png, dpi=150, bbox_inches="tight")
print(f"Saved: {out_png}")
plt.show()

# ── Summary ──────────────────────────────────────────────────────────────────
print("\n=== Summary ===")
print(f"Peak stress:    {float(load_stress_Pa.max()):.1f} Pa")
print(f"Final stress:   {float(relax_stress_Pa[-1]):.1f} Pa")
print(f"Stress ratio:   {float(relax_stress_Pa[-1])/float(load_stress_Pa.max()):.4f}  "
      f"(step-strain theory: {464/(464+23.2):.4f})")
if fit_ok:
    print(f"Fitted τ:       {tau_fit:.4f} h  (input: 1.0000 h)")
