"""
Shared utilities for viscoelastic rheology experiments.

Provides: path constants, build/run helpers, CSV I/O, SLS fitting,
plot styling, and error metrics.
"""

from __future__ import annotations

import os
import sys
import subprocess
import numpy as np

# ── Paths ─────────────────────────────────────────────────────────────────
BUILD_DIR = os.environ.get("CHASTE_BUILD_DIR",
                           os.path.expanduser("~/Thesis"))
TEST_OUTPUT = os.environ.get("CHASTE_TEST_OUTPUT",
                             os.path.expanduser("~/Thesis/testoutput"))
VC_DIR = os.path.join(TEST_OUTPUT, "ViscoelasticCompression")

EXE_2D = os.path.join(BUILD_DIR,
    "projects/TissueMorphology/test/TestViscoelasticShearExperiment")
EXE_3D = os.path.join(BUILD_DIR,
    "projects/TissueMorphology/test/TestViscoelasticShearExperiment3d")

# ── Default experiment parameters ─────────────────────────────────────────
ARM_RATIO = 0.2
TAU = 1.0

# ── Plot style ────────────────────────────────────────────────────────────
RCPARAMS = {
    "font.family":      "serif",
    "mathtext.fontset": "dejavuserif",
    "axes.labelsize":   11,
    "xtick.labelsize":  10,
    "ytick.labelsize":  10,
    "legend.fontsize":  10,
}

# ── Canonical output paths ────────────────────────────────────────────────

def shear_csv_path(dim: str, g_pa: int) -> str:
    """Path to the timeseries CSV: VC_DIR/SimpleShear_{dim}_G{g}/shear_data.csv"""
    return os.path.join(VC_DIR, f"SimpleShear_{dim}_G{g_pa}", "shear_data.csv")


def sweep_csv_path(dim: str) -> str:
    """Path to the sweep summary CSV: VC_DIR/shear_sweep_{dim}.csv"""
    return os.path.join(VC_DIR, f"shear_sweep_{dim}.csv")


# ── Build & run ───────────────────────────────────────────────────────────

def build_test(target: str) -> None:
    """Build a Chaste test target. Exits on failure."""
    print(f"Building {target} ...", flush=True)
    r = subprocess.run(["make", "-j4", target], cwd=BUILD_DIR,
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(f"BUILD FAILED:\n{r.stderr[-2000:]}")
        sys.exit(1)
    print("  OK")


def run_shear_test(exe: str, g_pa: float, out_dir: str,
                   skip_vtp: bool = True, timeout: int = 600) -> dict | None:
    """Run a shear experiment and return the parsed SHEAR_RESULT dict, or None."""
    env = os.environ.copy()
    env["SHEAR_G_PA"]       = str(g_pa)
    env["SHEAR_ARM_RATIO"]  = str(ARM_RATIO)
    env["SHEAR_TAU"]        = str(TAU)
    env["SHEAR_OUTPUT_DIR"] = out_dir
    env["SHEAR_SKIP_VTP"]   = "1" if skip_vtp else "0"

    r = subprocess.run([exe], env=env, capture_output=True, text=True,
                       timeout=timeout)
    if r.returncode != 0:
        print(f"  FAILED (G={g_pa} Pa): exit {r.returncode}")
        if r.stderr:
            print(r.stderr[-500:])
        return None

    for line in r.stdout.splitlines():
        if line.startswith("SHEAR_RESULT,"):
            p = line.strip().split(",")
            return {
                "G_target": float(p[1]),
                "G_relax":  float(p[2]),
                "G_peak":   float(p[3]),
            }

    print(f"  WARNING: no SHEAR_RESULT line for G={g_pa} Pa")
    return None


# ── CSV I/O ───────────────────────────────────────────────────────────────

def load_csv(path: str) -> dict[str, np.ndarray]:
    """Read a CSV into {column_name: numpy array}."""
    data = np.genfromtxt(path, delimiter=",", names=True)
    return {name: data[name] for name in data.dtype.names}


# ── SLS analytical solutions ──────────────────────────────────────────────

def sls_ramp_stress(t, G, arm_ratio, tau, gamma_dot):
    """Analytical SLS stress under constant shear strain rate.

    For an SLS model (spring G in parallel with Maxwell arm G1-eta),
    under constant strain rate gamma_dot:
        sigma(t) = G * gamma_dot * t  +  G1 * gamma_dot * tau * (1 - exp(-t/tau))

    where G1 = arm_ratio * G.

    Reference: Tirella et al., J. Biomed. Mater. Res. A, 2014, Appendix Eq. A4a
               (inverse Laplace of SLS transfer function under ramp input).
    """
    G1 = arm_ratio * G
    return G * gamma_dot * t + G1 * gamma_dot * tau * (1.0 - np.exp(-t / tau))


def sls_hold_stress(t, sigma_at_hold, G, gamma_max, tau):
    """Analytical SLS stress during hold (constant strain) after ramp loading.

    sigma(t) = G * gamma_max + (sigma_hold - G * gamma_max) * exp(-t/tau)

    where t is time since hold began and sigma_hold is the stress at the
    end of loading.

    Reference: Serra-Aguila et al., Viscoelastic models revisited, 2019,
               Eq. 6a (relaxation modulus of SLS).
    """
    sigma_eq = G * gamma_max
    return sigma_eq + (sigma_at_hold - sigma_eq) * np.exp(-t / tau)


def sls_full_protocol(time, G, arm_ratio, tau, gamma_dot, t_load):
    """Analytical SLS stress for full ramp-then-hold protocol.

    Returns stress array matching the time array. Uses sls_ramp_stress
    during loading and sls_hold_stress during relaxation.
    """
    gamma_max = gamma_dot * t_load
    stress = np.empty_like(time)

    mask_load = time <= t_load
    mask_hold = time > t_load

    stress[mask_load] = sls_ramp_stress(time[mask_load], G, arm_ratio, tau, gamma_dot)

    sigma_at_hold = sls_ramp_stress(t_load, G, arm_ratio, tau, gamma_dot)
    stress[mask_hold] = sls_hold_stress(
        time[mask_hold] - t_load, sigma_at_hold, G, gamma_max, tau)

    return stress


# ── SLS fitting ───────────────────────────────────────────────────────────

def sls_relax(t, sigma_inf, delta_sigma, tau):
    """Standard linear solid relaxation: sigma_inf + delta_sigma * exp(-t/tau)."""
    return sigma_inf + delta_sigma * np.exp(-t / tau)


def fit_sls(time: np.ndarray, stress: np.ndarray):
    """Fit SLS relaxation model. Returns (popt, fit_ok) where popt = (sigma_inf, delta_sigma, tau)."""
    from scipy.optimize import curve_fit
    try:
        p0 = [stress[-1], stress[0] - stress[-1], 1.0]
        popt, _ = curve_fit(sls_relax, time, stress, p0=p0)
        return popt, True
    except Exception as e:
        print(f"SLS fit failed: {e}")
        return None, False


# ── Error metrics ─────────────────────────────────────────────────────────

def compute_metrics(target: np.ndarray, measured: np.ndarray) -> dict:
    """Compute RMSE, MAE, mean ratio, and mean error between target and measured."""
    errors = measured - target
    return {
        "mse":        float(np.mean(errors**2)),
        "rmse":       float(np.sqrt(np.mean(errors**2))),
        "mae":        float(np.mean(np.abs(errors))),
        "mean_ratio": float(np.mean(measured / target)),
        "mean_error": float(np.mean(errors)),
        "mean_error_pct": float(np.mean(errors / target) * 100),
    }


def print_metrics(m: dict, n_points: int, g_range: str = "") -> None:
    """Print a formatted summary of sweep error metrics."""
    print(f"\n{'='*60}")
    print(f"SWEEP SUMMARY ({n_points} points{', G = ' + g_range if g_range else ''})")
    print(f"{'='*60}")
    print(f"  MSE  = {m['mse']:.1f} Pa\u00b2")
    print(f"  RMSE = {m['rmse']:.1f} Pa")
    print(f"  MAE  = {m['mae']:.1f} Pa")
    print(f"  Mean G_relax/G_target = {m['mean_ratio']:.4f}")
    print(f"  Mean error = {m['mean_error']:.1f} Pa ({m['mean_error_pct']:.1f}%)")
    print()
