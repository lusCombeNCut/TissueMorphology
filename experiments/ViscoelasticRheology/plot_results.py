#!/usr/bin/env python3
"""
Unified plotter for viscoelastic rheology experiment results.

Works entirely from CSV data — no simulations are re-run.
Expects the naming convention produced by run_experiments.py:
  Timeseries:    VC_DIR/SimpleShear_{dim}_G{g}/shear_data.csv
  Sweep summary: VC_DIR/shear_sweep_{dim}.csv

Subcommands
-----------
  sweep-scatter       Measured-vs-target scatter for a single sweep CSV.
  sweep-compare       Side-by-side 2D vs 3D sweep scatter plots.
  timeseries-compare  2D vs 3D shear stress vs time at specific G values.
  comparison-grid     2x2 grid: time-series (top) + recovery bars (bottom).
  shear-detail        3-panel shear analysis with SLS fit (single test).
  compression-detail  3-panel compression analysis with SLS fit (single test).

Examples
--------
  python plot_results.py sweep-scatter
  python plot_results.py sweep-compare
  python plot_results.py timeseries-compare --g 300 1300
  python plot_results.py comparison-grid --g 300 800 1300
  python plot_results.py shear-detail --csv path/to/shear_data.csv
  python plot_results.py compression-detail --csv path/to/compression_data.csv
"""

import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from common import (
    VC_DIR, RCPARAMS, ARM_RATIO, TAU,
    load_csv, shear_csv_path, sweep_csv_path,
    sls_relax, fit_sls, compute_metrics,
    sls_ramp_stress, sls_hold_stress, sls_full_protocol,
)


def _apply_style():
    plt.rcParams.update(RCPARAMS)


def _save(fig, path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"Saved: {path}")
    plt.close(fig)


def _load_sweep(path: str):
    """Load a sweep CSV, return (G_target, G_relax, G_peak) arrays."""
    d = load_csv(path)
    return d["G_target_Pa"], d["G_relax_Pa"], d["G_peak_Pa"]


def _require(path: str):
    if not os.path.exists(path):
        print(f"ERROR: CSV not found: {path}")
        sys.exit(1)


# ── sweep-scatter ─────────────────────────────────────────────────────────

def cmd_sweep_scatter(args):
    _apply_style()
    csv_path = args.csv or sweep_csv_path(args.dim.upper())
    _require(csv_path)

    G_target, G_relax, G_peak = _load_sweep(csv_path)

    m = compute_metrics(G_target, G_relax)
    print(f"Loaded {len(G_target)} points from {csv_path}")
    print(f"  RMSE = {m['rmse']:.1f} Pa  |  MAE = {m['mae']:.1f} Pa  "
          f"|  mean ratio = {m['mean_ratio']:.4f}")

    fig, ax = plt.subplots(figsize=(3.5, 2.5))
    ax.plot([0, 2000], [0, 2000], "k--", lw=1, label="$y = x$")
    ax.scatter(G_target, G_relax, c="tab:blue", s=60, zorder=3,
               label=r"$G_\mathrm{relax}$")
    ax.scatter(G_target, G_peak, c="tab:orange", s=40, marker="^", zorder=3,
               alpha=0.8, label=r"$G_\mathrm{peak}$")
    ax.set_xlabel("Target $G$ (Pa)")
    ax.set_ylabel("Measured $G$ (Pa)")
    ax.legend()
    ax.set_xlim(0, 2100)
    ax.set_ylim(0, 2100)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    _save(fig, args.out or os.path.join(VC_DIR, "shear_sweep_validation.svg"))


# ── sweep-compare ─────────────────────────────────────────────────────────

def cmd_sweep_compare(args):
    _apply_style()

    csv_2d = args.csv_2d or sweep_csv_path("2D")
    csv_3d = args.csv_3d or sweep_csv_path("3D")
    _require(csv_2d)
    _require(csv_3d)

    panels = [
        ("2D lattice", *_load_sweep(csv_2d)),
        ("3D lattice", *_load_sweep(csv_3d)),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.3), sharey=True)

    for ax, (title, gt, gr, gp) in zip(axes, panels):
        ax.plot([0, 4000], [0, 4000], "k--", lw=1, label="$y = x$")
        ax.scatter(gt, gr, c="tab:blue", s=55, zorder=3,
                   label=r"$G_\mathrm{relax}$")
        ax.scatter(gt, gp, c="tab:orange", s=40, marker="^", zorder=3,
                   alpha=0.85, label=r"$G_\mathrm{peak}$")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Target $G$ (Pa)")
        ax.set_xlim(0, 2100)
        ax.set_ylim(0, 4000)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper left")

    axes[0].set_ylabel("Measured $G$ (Pa)")
    plt.tight_layout()

    _save(fig, args.out or os.path.join(VC_DIR, "shear_sweep_2d_3d_scatter.svg"))


# ── timeseries-compare ────────────────────────────────────────────────────

def cmd_timeseries_compare(args):
    _apply_style()

    g_values = args.g
    arm_ratio = args.arm_ratio
    tau = args.tau

    # First G value: solid teal; second G value: dashed orange
    STYLE_G = [
        {"color": "#1B998B", "ls": "-"},
        {"color": "#9B2226", "ls": "--"},
    ]
    STYLE_ANALYTICAL = [
        {"color": "#1B998B"},
        {"color": "#9B2226"},
    ]

    t_load = 1.0
    gamma_dot = 0.05 / t_load

    fig, axes = plt.subplots(1, 2, figsize=(8, 3.5), sharey=True)

    for ax, dim in zip(axes, ["2D", "3D"]):
        ax.set_title(f"{dim} lattice", fontsize=11)
        ax.set_xlabel("Time (h)")
        ax.grid(True, alpha=0.3)
        ax.axvline(x=1.0, color="gray", lw=0.7, ls=":", label="End of loading")

        for i, G in enumerate(g_values):
            style = STYLE_G[i] if i < len(STYLE_G) else {"color": f"C{i}", "ls": "-"}
            path = shear_csv_path(dim, G)
            if not os.path.exists(path):
                print(f"  WARNING: {path} not found, skipping")
                continue
            d = load_csv(path)
            ax.plot(d["time"], d["shear_stress_Pa"],
                    lw=1.6, label=f"{G} Pa", **style)

            # Analytical overlay
            t_max = min(float(d["time"][-1]), 3.0)
            t_an = np.linspace(0, t_max, 1000)
            sigma_an = sls_full_protocol(t_an, G, arm_ratio, tau, gamma_dot, t_load)
            a_style = STYLE_ANALYTICAL[i] if i < len(STYLE_ANALYTICAL) else {"color": f"C{i}"}
            ax.plot(t_an, sigma_an, ls=":", lw=1.2, alpha=0.9,
                    label=f"{G} Pa (analytical)", **a_style)

        ax.set_xlim(0, 3)
        ax.legend(fontsize=8)

    axes[0].set_ylabel("Shear stress (Pa)")
    plt.tight_layout()
    _save(fig, args.out or os.path.join(VC_DIR, "shear_timeseries_2d_3d.svg"))


# ── comparison-grid ───────────────────────────────────────────────────────

def cmd_comparison_grid(args):
    _apply_style()

    g_values = args.g
    COLORS = ["#2176AE", "#57B8FF", "#B66D0D"]
    while len(COLORS) < len(g_values):
        COLORS.append(f"C{len(COLORS)}")
    LABELS = [f"{G} Pa" for G in g_values]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("Simple Shear Stress-Relaxation: 2D vs 3D Ghost Node Network",
                 fontsize=13, fontweight="bold", y=0.98)

    ax_ts_2d, ax_ts_3d = axes[0, 0], axes[0, 1]
    ax_bar_2d, ax_bar_3d = axes[1, 0], axes[1, 1]

    # Timeseries panels
    for dim_label, ax, dim in [
        ("2D (8-connected square)", ax_ts_2d, "2D"),
        ("3D (18-connected cubic)", ax_ts_3d, "3D"),
    ]:
        ax.set_title(dim_label, fontsize=11)
        ax.set_xlabel("Time (h)")
        ax.set_ylabel("Shear stress (Pa)")

        for G, col, lab in zip(g_values, COLORS, LABELS):
            path = shear_csv_path(dim, G)
            if not os.path.exists(path):
                print(f"  WARNING: {path} not found, skipping")
                continue
            d = load_csv(path)
            ax.plot(d["time"], d["shear_stress_Pa"],
                    color=col, linewidth=1.4, label=lab)

        ax.axvline(x=1.0, color="gray", linestyle="--", linewidth=0.8,
                   label="End of loading")
        ax.legend(fontsize=8, loc="upper right")
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(
            lambda v, _: f"{v:.0f}" if v >= 1 else f"{v:.1f}"))
        ax.grid(True, alpha=0.3, linewidth=0.5)

    # Recovery bar charts
    x_pos = np.arange(len(g_values))

    for dim_label, ax, dim in [("2D", ax_bar_2d, "2D"), ("3D", ax_bar_3d, "3D")]:
        ax.set_title(f"G recovery \u2014 {dim_label}", fontsize=11)
        ax.set_xlabel("Target G (Pa)")
        ax.set_ylabel("G_relax / G_target")
        ax.set_xticks(x_pos)
        ax.set_xticklabels([str(G) for G in g_values])
        ax.axhline(y=1.0, color="black", linewidth=0.8, linestyle="--")
        ax.set_ylim(0, 1.25)
        ax.grid(True, axis="y", alpha=0.3, linewidth=0.5)

        ratios = []
        bar_cols = []
        for G, col in zip(g_values, COLORS):
            path = shear_csv_path(dim, G)
            if os.path.exists(path):
                d = load_csv(path)
                gamma_max = d["shear_strain"][d["time"] <= 1.0].max()
                G_relax = d["shear_stress_Pa"][-1] / gamma_max
                ratios.append(G_relax / G)
            else:
                ratios.append(0.0)
            bar_cols.append(col)

        bars = ax.bar(x_pos, ratios, width=0.7, color=bar_cols,
                      edgecolor="k", linewidth=0.6, alpha=0.85)
        for bar, ratio in zip(bars, ratios):
            if ratio > 0:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        ratio + 0.02, f"{100*ratio:.1f}%",
                        ha="center", va="bottom", fontsize=9)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    _save(fig, args.out or os.path.join(VC_DIR, "shear_2d_3d_comparison.svg"))


# ── shear-detail ──────────────────────────────────────────────────────────

def cmd_shear_detail(args):
    csv_path = args.csv or os.path.join(VC_DIR, "SimpleShear_2D", "shear_data.csv")
    _require(csv_path)

    d = load_csv(csv_path)
    G_input = args.g_input
    arm_ratio = args.arm_ratio
    tau = args.tau
    t_load_end = 1.0

    mask_load = d["time"] <= t_load_end
    mask_relax = d["time"] > t_load_end
    gamma_max = float(d["shear_strain"][mask_load][-1])
    gamma_dot = gamma_max / t_load_end

    relax_time = d["time"][mask_relax] - t_load_end
    relax_stress = d["shear_stress_Pa"][mask_relax]

    popt, fit_ok = fit_sls(relax_time, relax_stress)

    peak_stress = float(d["shear_stress_Pa"][mask_load].max())
    final_stress = float(relax_stress[-1])
    G_peak = peak_stress / gamma_max
    G_relaxed = final_stress / gamma_max
    G_relaxed_fit = popt[0] / gamma_max if fit_ok else None

    # Analytical SLS solution for full protocol
    t_analytical = np.linspace(0, float(d["time"][-1]), 2000)
    sigma_analytical = sls_full_protocol(
        t_analytical, G_input, arm_ratio, tau, gamma_dot, t_load_end)

    # Analytical during loading only (for stress-strain panel)
    t_load_fine = np.linspace(0, t_load_end, 500)
    sigma_load_analytical = sls_ramp_stress(
        t_load_fine, G_input, arm_ratio, tau, gamma_dot)
    gamma_load_analytical = gamma_dot * t_load_fine

    # Figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    fig.suptitle(f"Simple Shear Test \u2014 Viscoelastic Ghost Node ECM\n"
                 f"Input G = {G_input:.0f} Pa", fontsize=11)

    # Panel A: Stress-strain (loading)
    ax = axes[0]
    load_strain = d["shear_strain"][mask_load]
    load_stress = d["shear_stress_Pa"][mask_load]
    ax.plot(load_strain * 100, load_stress, "C0", lw=2, label="Simulation")
    ax.plot(gamma_load_analytical * 100, sigma_load_analytical,
            "C1--", lw=1.5, label="SLS analytical")
    gamma_line = np.linspace(0, float(load_strain.max()), 100)
    ax.plot(gamma_line * 100, G_input * gamma_line, "C3:", lw=1.2,
            label=f"Elastic ($G$ = {G_input:.0f} Pa)")
    ax.set_xlabel("Shear strain  $\\gamma$  (%)")
    ax.set_ylabel("Shear stress  $\\tau$  (Pa)")
    ax.set_title("(A) Shear Stress\u2013Strain (loading)")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # Panel B: Stress vs time (with analytical overlay)
    ax = axes[1]
    ax.plot(d["time"], d["shear_stress_Pa"], "C0", lw=1.5, label="Simulation")
    ax.plot(t_analytical, sigma_analytical, "C1--", lw=1.5,
            label="SLS analytical")
    ax.axvline(t_load_end, color="k", ls="--", lw=0.8, label="End loading")
    ax.set_xlabel("Time  (h)")
    ax.set_ylabel("Shear stress  $\\tau$  (Pa)")
    ax.set_title("(B) Shear Stress vs Time")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # Panel C: Relaxation with fit and analytical
    ax = axes[2]
    ax.plot(relax_time, relax_stress, "C0", lw=1.5, label="Simulation")

    # Analytical relaxation
    sigma_at_hold = sls_ramp_stress(t_load_end, G_input, arm_ratio, tau, gamma_dot)
    t_relax_fine = np.linspace(0, float(relax_time[-1]), 500)
    sigma_relax_analytical = sls_hold_stress(
        t_relax_fine, sigma_at_hold, G_input, gamma_max, tau)
    ax.plot(t_relax_fine, sigma_relax_analytical, "C1--", lw=1.5,
            label="SLS analytical")

    if fit_ok:
        t_fine = np.linspace(0, float(relax_time[-1]), 500)
        ax.plot(t_fine, sls_relax(t_fine, *popt), "C2:", lw=1.5, label="SLS fit")
        ax.axhline(popt[0], color="grey", ls=":", lw=1,
                   label=f"$\\tau_\\infty$ = {popt[0]:.1f} Pa")
        fit_label = (f"Fit: $\\tau(t)$ = {popt[0]:.1f} + {popt[1]:.1f}\u00b7exp(-t/{popt[2]:.3f})\n"
                     f"$\\tau_{{fit}}$ = {popt[2]:.3f} h  (input $\\tau$ = {tau:.3f} h)\n"
                     f"$G_{{relax}}$(fit) = {G_relaxed_fit:.0f} Pa")
        ax.text(0.97, 0.97, fit_label, transform=ax.transAxes,
                fontsize=7, ha="right", va="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    ax.set_xlabel("Time after hold  (h)")
    ax.set_ylabel("Shear stress  $\\tau$  (Pa)")
    ax.set_title("(C) Shear Relaxation with SLS Fit")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    _save(fig, args.out or os.path.join(os.path.dirname(csv_path), "shear_relaxation.png"))

    # Summary
    print(f"\n=== Shear Modulus Recovery Summary ===")
    print(f"Input G:             {G_input:.0f} Pa")
    print(f"Max shear strain:    {gamma_max:.4f} ({gamma_max*100:.1f}%)")
    print(f"Peak shear stress:   {peak_stress:.1f} Pa")
    print(f"Final shear stress:  {final_stress:.1f} Pa")
    print(f"G_peak  (\u03c3_peak/\u03b3):  {G_peak:.0f} Pa   (ratio: {G_peak/G_input:.1f}x)")
    print(f"G_relax (\u03c3_final/\u03b3): {G_relaxed:.0f} Pa   (ratio: {G_relaxed/G_input:.1f}x)")
    if fit_ok:
        print(f"G_relax (from fit):  {G_relaxed_fit:.0f} Pa   (ratio: {G_relaxed_fit/G_input:.1f}x)")
        print(f"Fitted \u03c4:            {popt[2]:.4f} h  (input: 1.0000 h)")


# ── compression-detail ────────────────────────────────────────────────────

def cmd_compression_detail(args):
    csv_path = args.csv or os.path.join(VC_DIR, "StressRelaxation_2D", "compression_data.csv")
    _require(csv_path)

    d = load_csv(csv_path)
    # Physical conversion: sim stress -> Pa
    kF = 0.001       # N/m per sim unit
    spacing_m = 1.0 * 10e-6  # 1 sim unit = 10 um
    stress_Pa = d["stress"] * kF / spacing_m

    t_load_end = 1.0
    mask_load = d["time"] <= t_load_end
    mask_relax = d["time"] > t_load_end

    relax_time = d["time"][mask_relax] - t_load_end
    relax_stress = stress_Pa[mask_relax]

    popt, fit_ok = fit_sls(relax_time, relax_stress)
    fit_label = (f"Fit: \u03c3(t) = {popt[0]:.0f} + {popt[1]:.0f}\u00b7exp(-t/{popt[2]:.3f})\n"
                 f"\u03c4_fit = {popt[2]:.3f} h  (input \u03c4 = 1.000 h)"
                 if fit_ok else "Fit failed")

    # Figure
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    fig.suptitle("Viscoelastic Ghost Node ECM \u2014 Uniaxial Confined Compression\n"
                 "G = 800 Pa  |  E\u2080 = 464 sim  |  E\u2081 = 23.2 sim  |  \u03c4 = 1 h",
                 fontsize=11)

    ax = axes[0]
    ax.plot(d["strain"][mask_load] * 100, stress_Pa[mask_load], "C0", lw=2)
    ax.set_xlabel("Engineering strain  \u03b5  (%)")
    ax.set_ylabel("Compressive stress  \u03c3  (Pa)")
    ax.set_title("(A) Stress\u2013Strain (loading)")
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.plot(d["time"], stress_Pa, "C0", lw=1.5, label="Simulation")
    ax.axvline(t_load_end, color="k", ls="--", lw=0.8, label="End loading")
    ax.set_xlabel("Time  (h)")
    ax.set_ylabel("Compressive stress  \u03c3  (Pa)")
    ax.set_title("(B) Stress vs Time (full protocol)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[2]
    ax.plot(relax_time, relax_stress, "C0", lw=1.5, label="Simulation")
    if fit_ok:
        t_fine = np.linspace(0, relax_time[-1], 500)
        ax.plot(t_fine, sls_relax(t_fine, *popt), "C1--", lw=1.5, label="SLS fit")
        ax.axhline(popt[0], color="grey", ls=":", lw=1,
                   label=f"\u03c3_\u221e = {popt[0]:.0f} Pa")
    ax.set_xlabel("Time after hold  (h)")
    ax.set_ylabel("Compressive stress  \u03c3  (Pa)")
    ax.set_title("(C) Stress Relaxation with SLS Fit")
    ax.legend(fontsize=7)
    ax.text(0.97, 0.97, fit_label, transform=ax.transAxes,
            fontsize=7, ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    _save(fig, args.out or os.path.join(os.path.dirname(csv_path), "stress_relaxation.png"))

    peak = float(stress_Pa[mask_load].max())
    final = float(relax_stress[-1])
    print(f"\n=== Summary ===")
    print(f"Peak stress:    {peak:.1f} Pa")
    print(f"Final stress:   {final:.1f} Pa")
    print(f"Stress ratio:   {final/peak:.4f}  "
          f"(step-strain theory: {464/(464+23.2):.4f})")
    if fit_ok:
        print(f"Fitted \u03c4:       {popt[2]:.4f} h  (input: 1.0000 h)")


# ── CLI ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Plot viscoelastic rheology results (no re-running)")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("sweep-scatter", help="Measured-vs-target scatter")
    p.add_argument("--csv", help="Path to sweep CSV")
    p.add_argument("--dim", default="2D", help="Dimension label for default path")
    p.add_argument("--out")
    p.set_defaults(func=cmd_sweep_scatter)

    p = sub.add_parser("sweep-compare", help="Side-by-side 2D vs 3D scatter")
    p.add_argument("--csv-2d")
    p.add_argument("--csv-3d")
    p.add_argument("--out")
    p.set_defaults(func=cmd_sweep_compare)

    p = sub.add_parser("timeseries-compare", help="2D vs 3D stress vs time + analytical")
    p.add_argument("--g", type=int, nargs="+", default=[300, 1300])
    p.add_argument("--arm-ratio", type=float, default=ARM_RATIO)
    p.add_argument("--tau", type=float, default=TAU)
    p.add_argument("--out")
    p.set_defaults(func=cmd_timeseries_compare)

    p = sub.add_parser("comparison-grid", help="2x2 grid: timeseries + bars")
    p.add_argument("--g", type=int, nargs="+", default=[300, 800, 1300])
    p.add_argument("--out")
    p.set_defaults(func=cmd_comparison_grid)

    p = sub.add_parser("shear-detail", help="3-panel shear with SLS fit + analytical")
    p.add_argument("--csv")
    p.add_argument("--g-input", type=float, default=800.0)
    p.add_argument("--arm-ratio", type=float, default=ARM_RATIO)
    p.add_argument("--tau", type=float, default=TAU)
    p.add_argument("--out")
    p.set_defaults(func=cmd_shear_detail)

    p = sub.add_parser("compression-detail", help="3-panel compression with SLS fit")
    p.add_argument("--csv")
    p.add_argument("--out")
    p.set_defaults(func=cmd_compression_detail)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
