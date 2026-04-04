#!/usr/bin/env python3
"""
Debug script for crypt counting analysis.

Loads sample simulation data, runs the SimpleCryptCount pipeline step-by-step,
and produces diagnostic images showing where/why crypt detection fails.

Usage:
  python debug_crypt_counting.py
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path

# Add this directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from simple_crypt_count import (
    load_outline_from_vtp, load_final_outline, load_final_vertex_boundary,
    compute_chain_code, fourier_reconstruct, direct_fourier_reconstruct,
    compute_curvature, compute_normals,
    find_crypt_sections, filter_crypts, compute_circularity, compute_section_area,
    compute_arc_length, count_crypts_simple_method, DEFAULT_PARAMS,
    SIMULATION_PARAMS, PAPER_PARAMS, _extract_outer_boundary
)

# ===========================================================================
# Configuration
# ===========================================================================
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    'results-download', '03-31_08_56_51-merged')
NODE2D_DIR = os.path.join(BASE, '16535750_node2d_2026-03-31_08-56-51')
VERTEX2D_DIR = os.path.join(BASE, '16535751_vertex2d_2026-03-31_08-56-51')

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'debug_output')
os.makedirs(OUTPUT_DIR, exist_ok=True)


def load_node2d_outline(stiffness, run):
    """Load the final outline for a node2d simulation."""
    sim_dir = os.path.join(NODE2D_DIR, f's{stiffness}_r{run}', 'results_from_time_0')
    return load_final_outline(sim_dir)


def load_vertex2d_boundary(stiffness, run, target_time=35.0):
    """Load the boundary for a vertex2d simulation at t=35h (avoids late-stage folding)."""
    sim_dir = os.path.join(VERTEX2D_DIR, f's{stiffness}_r{run}', 'results_from_time_0')
    return load_final_vertex_boundary(sim_dir, target_time=target_time, dt=0.001)


# ===========================================================================
# Diagnostic 1: Inspect raw boundary data and chain code quality
# ===========================================================================
def diagnose_chain_code(boundary, title=""):
    """Show how chain code discretization affects the boundary."""
    chain = compute_chain_code(boundary)

    # Compute actual angles between consecutive points
    diff = np.vstack([np.diff(boundary, axis=0), boundary[0] - boundary[-1]])
    actual_angles = np.arctan2(diff[:, 1], diff[:, 0])
    discretized_angles = chain * np.pi / 4

    # Compute step lengths
    step_lengths = np.linalg.norm(diff, axis=1)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. Raw boundary with point indices
    ax = axes[0, 0]
    ax.plot(boundary[:, 0], boundary[:, 1], 'b-o', markersize=3)
    ax.plot([boundary[-1, 0], boundary[0, 0]], [boundary[-1, 1], boundary[0, 1]], 'b-o', markersize=3)
    for i in range(0, len(boundary), max(1, len(boundary) // 20)):
        ax.annotate(str(i), boundary[i], fontsize=6)
    ax.set_title(f'Raw boundary ({len(boundary)} points)')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    # 2. Step lengths
    ax = axes[0, 1]
    ax.bar(range(len(step_lengths)), step_lengths, width=1.0)
    ax.set_title(f'Step lengths (mean={step_lengths.mean():.3f}, std={step_lengths.std():.3f})')
    ax.set_xlabel('Point index')
    ax.set_ylabel('Distance to next point')
    ax.axhline(y=step_lengths.mean(), color='r', linestyle='--', label='mean')
    ax.legend()

    # 3. Actual vs discretized angles
    ax = axes[1, 0]
    ax.plot(actual_angles, 'b-', alpha=0.7, label='Actual angles')
    ax.plot(discretized_angles, 'r.', markersize=3, alpha=0.7, label='Discretized (8-dir)')
    ax.set_title('Angle discretization error')
    ax.set_xlabel('Point index')
    ax.set_ylabel('Angle (radians)')
    ax.legend()

    # 4. Angle error
    ax = axes[1, 1]
    error = np.abs(actual_angles - discretized_angles)
    error = np.minimum(error, 2 * np.pi - error)  # Handle wraparound
    ax.hist(error, bins=30)
    ax.set_title(f'Discretization error (mean={np.degrees(error.mean()):.1f}°, max={np.degrees(error.max()):.1f}°)')
    ax.set_xlabel('Angle error (radians)')
    ax.set_ylabel('Count')

    fig.suptitle(f'{title} — Chain Code Diagnostics', fontsize=14, fontweight='bold')
    plt.tight_layout()
    path = os.path.join(OUTPUT_DIR, f'{title.replace(" ", "_")}_chain_code.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")


# ===========================================================================
# Diagnostic 2: Compare Fourier reconstruction quality
# ===========================================================================
def diagnose_fourier(boundary, title=""):
    """Compare Fourier reconstruction with original boundary."""
    # Reconstruct with current method (chain code based, normalized)
    smooth_cc = fourier_reconstruct(boundary, n_harmonics=25, n_points=1000, normalize=True)

    # Reconstruct WITHOUT normalization
    smooth_cc_nonorm = fourier_reconstruct(boundary, n_harmonics=25, n_points=1000, normalize=False)

    # Direct Fourier fit (bypass chain code)
    smooth_direct = direct_fourier_reconstruct(boundary, n_harmonics=25, n_points=1000)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    boundary_closed = np.vstack([boundary, boundary[0:1]])

    for ax, smooth, label in [(axes[0], smooth_cc, 'Chain code + normalized'),
                               (axes[1], smooth_cc_nonorm, 'Chain code, no normalize'),
                               (axes[2], smooth_direct, 'Direct Fourier (no chain code)')]:
        ax.plot(boundary_closed[:, 0], boundary_closed[:, 1], 'b--', lw=1, alpha=0.5, label='Original')
        if len(smooth) > 0:
            smooth_closed = np.vstack([smooth, smooth[0:1]])
            ax.plot(smooth_closed[:, 0], smooth_closed[:, 1], 'r-', lw=1.5, label='Reconstructed')
            circ = compute_circularity(smooth)
            ax.set_title(f'{label}\nCircularity={circ:.4f}, {len(smooth)} pts')
        else:
            ax.set_title(f'{label}\nFAILED')
        ax.set_aspect('equal')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    fig.suptitle(f'{title} — Fourier Reconstruction Comparison', fontsize=14, fontweight='bold')
    plt.tight_layout()
    path = os.path.join(OUTPUT_DIR, f'{title.replace(" ", "_")}_fourier.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")


# ===========================================================================
# Diagnostic 3: Full pipeline comparison
# ===========================================================================
def diagnose_full_pipeline(boundary, title="", params=None):
    """
    Run the full crypt detection pipeline with both chain-code and direct Fourier,
    showing each step side by side.
    """
    if params is None:
        params = DEFAULT_PARAMS.copy()

    # Method 1: Current (chain code)
    smooth_cc = fourier_reconstruct(boundary, n_harmonics=25, n_points=1000, normalize=True)
    curv_cc = compute_curvature(smooth_cc)
    norm_cc = compute_normals(smooth_cc)

    # Method 2: Direct Fourier
    smooth_df = direct_fourier_reconstruct(boundary, n_harmonics=25, n_points=1000)
    curv_df = compute_curvature(smooth_df)
    norm_df = compute_normals(smooth_df)

    fig, axes = plt.subplots(2, 4, figsize=(24, 12))

    for row, (smooth, curv, norms, method_name) in enumerate([
        (smooth_cc, curv_cc, norm_cc, 'Chain Code Fourier'),
        (smooth_df, curv_df, norm_df, 'Direct Fourier')
    ]):
        if len(smooth) == 0:
            for col in range(4):
                axes[row, col].set_title(f'{method_name}: FAILED')
            continue

        smooth_closed = np.vstack([smooth, smooth[0:1]])
        circ = compute_circularity(smooth)

        # Col 0: Smooth boundary
        ax = axes[row, 0]
        boundary_closed = np.vstack([boundary, boundary[0:1]])
        ax.plot(boundary_closed[:, 0], boundary_closed[:, 1], 'b--', lw=1, alpha=0.4, label='Original')
        ax.plot(smooth_closed[:, 0], smooth_closed[:, 1], 'k-', lw=1, label='Smooth')
        ax.set_title(f'{method_name}\nCircularity={circ:.4f}')
        ax.set_aspect('equal')
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

        # Col 1: Curvature
        ax = axes[row, 1]
        vmax = max(abs(curv.min()), abs(curv.max()), 0.01)
        sc = ax.scatter(smooth[:, 0], smooth[:, 1], c=curv, cmap='RdBu_r',
                       s=3, vmin=-vmax, vmax=vmax)
        plt.colorbar(sc, ax=ax, label='Curvature')
        ax.set_title(f'Curvature\nmin={curv.min():.4f}, max={curv.max():.4f}')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)

        # Col 2: Curvature-scaled normals + inpolygon
        ax = axes[row, 2]
        scaled_endpoints = smooth + curv[:, np.newaxis] * norms
        polygon_path = Path(smooth)
        inside = polygon_path.contains_points(scaled_endpoints)

        ax.plot(smooth_closed[:, 0], smooth_closed[:, 1], 'k-', lw=0.5)
        step = max(1, len(smooth) // 80)
        for i in range(0, len(smooth), step):
            pt = smooth[i]
            ep = scaled_endpoints[i]
            color = 'green' if inside[i] else 'red'
            ax.plot([pt[0], ep[0]], [pt[1], ep[1]], color=color, lw=0.7, alpha=0.6)

        n_inside = np.sum(inside)
        ax.set_title(f'Normal tips: {n_inside} inside, {len(inside)-n_inside} outside')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)

        # Col 3: Detected sections + filtering
        ax = axes[row, 3]
        sections, inside_mask = find_crypt_sections(smooth, curv, norms, params, return_mask=True)
        valid_crypts = filter_crypts(smooth, sections, params)

        ax.plot(smooth_closed[:, 0], smooth_closed[:, 1], 'k-', lw=1)
        ax.fill(smooth_closed[:, 0], smooth_closed[:, 1], color='lightgray', alpha=0.3)

        # Show all sections
        total_area = compute_section_area(smooth, list(range(len(smooth))))
        total_arc = compute_arc_length(smooth_closed)

        colors = plt.cm.tab10(np.linspace(0, 1, max(len(sections), 1)))
        section_info = []
        for i, sec in enumerate(sections):
            pts = smooth[sec]
            ax.plot(pts[:, 0], pts[:, 1], color=colors[i % 10], lw=3, alpha=0.7)
            sec_area = compute_section_area(smooth, sec)
            sec_arc = compute_arc_length(pts)
            na = sec_area / total_area if total_area > 0 else 0
            nl = sec_arc / total_arc if total_arc > 0 else 0
            section_info.append(f'S{i+1}: A={na:.4f} L={nl:.4f} n={len(sec)}')

        # Highlight valid crypts
        for crypt in valid_crypts:
            pts = smooth[crypt['indices']]
            ax.fill(pts[:, 0], pts[:, 1], color='lime', alpha=0.5, edgecolor='darkgreen', lw=2)

        info = '\n'.join(section_info[:15])
        ax.text(0.02, 0.98, info, transform=ax.transAxes, fontsize=6,
                va='top', ha='left', family='monospace',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))

        filter_text = f"Filters: A∈[{params['min_area']:.3f},{params['max_area']:.3f}) L≥{params['min_arc_length']:.3f}"
        ax.set_title(f'{len(sections)} sections → {len(valid_crypts)} crypts\n{filter_text}')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)

    fig.suptitle(f'{title} — Full Pipeline Comparison', fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    path = os.path.join(OUTPUT_DIR, f'{title.replace(" ", "_")}_pipeline.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")


# ===========================================================================
# Diagnostic 4: Vertex2d boundary extraction issue
# ===========================================================================
def diagnose_vertex2d_boundary(stiffness, run, title=""):
    """Diagnose the vertex2d boundary extraction at t=35h."""
    import glob

    sim_dir = os.path.join(VERTEX2D_DIR, f's{stiffness}_r{run}', 'results_from_time_0')

    # Load at t=35h to avoid late-stage folding
    try:
        boundary, cell_types = load_final_vertex_boundary(sim_dir, target_time=35.0, dt=0.001)
        print(f"  Loaded boundary: {len(boundary)} points")
    except Exception as e:
        print(f"  ERROR loading boundary: {e}")
        import traceback; traceback.print_exc()
        boundary = None

    # Also load cell positions at the same time for comparison
    from analyse_crypt_budding import load_positions_from_vtu
    target_ts = int(35.0 / 0.001)
    vtu_files = [f for f in glob.glob(os.path.join(sim_dir, 'results_*.vtu'))
                 if 'ecm_grid' not in f]
    best_vtu = min(vtu_files,
                   key=lambda p: abs(int(os.path.basename(p).replace('results_', '').replace('.vtu', '')) - target_ts))
    print(f"  Using VTU: {os.path.basename(best_vtu)} (target t=35h)")
    positions = load_positions_from_vtu(best_vtu, dim=2)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # 1. Raw cell positions (centroids)
    ax = axes[0]
    ax.scatter(positions[:, 0], positions[:, 1], s=20, c='blue', alpha=0.6)
    ax.set_title(f'Cell centroids ({len(positions)} cells)')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    # 2. Extracted boundary
    ax = axes[1]
    if boundary is not None:
        boundary_closed = np.vstack([boundary, boundary[0:1]])
        ax.plot(boundary_closed[:, 0], boundary_closed[:, 1], 'r-', lw=1.5)
        ax.scatter(boundary[:, 0], boundary[:, 1], s=5, c='red', zorder=5)
        # Check for self-intersections
        n_crossings = count_self_intersections(boundary)
        ax.set_title(f'Extracted boundary ({len(boundary)} pts)\n{n_crossings} self-intersections')
    else:
        ax.set_title('Boundary extraction FAILED')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    # 3. Convex hull + alpha shape from cell centroids
    ax = axes[2]
    from simple_crypt_count import extract_boundary_alpha_shape
    alpha_boundary = extract_boundary_alpha_shape(positions)
    if len(alpha_boundary) > 0:
        ab_closed = np.vstack([alpha_boundary, alpha_boundary[0:1]])
        ax.plot(ab_closed[:, 0], ab_closed[:, 1], 'g-', lw=1.5, label='Alpha shape')
    ax.scatter(positions[:, 0], positions[:, 1], s=10, c='blue', alpha=0.4)
    ax.set_title(f'Alpha shape from centroids ({len(alpha_boundary)} pts)')
    ax.set_aspect('equal')
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.suptitle(f'{title} — Vertex2d Boundary Diagnosis', fontsize=14, fontweight='bold')
    plt.tight_layout()
    path = os.path.join(OUTPUT_DIR, f'{title.replace(" ", "_")}_vertex2d_boundary.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")

    return boundary, positions


def count_self_intersections(boundary):
    """Count approximate self-intersections in a polygon."""
    n = len(boundary)
    crossings = 0
    for i in range(n):
        p1 = boundary[i]
        p2 = boundary[(i + 1) % n]
        for j in range(i + 2, n):
            if j == (i - 1) % n or (i == 0 and j == n - 1):
                continue  # Skip adjacent edges
            p3 = boundary[j]
            p4 = boundary[(j + 1) % n]
            if segments_intersect(p1, p2, p3, p4):
                crossings += 1
    return crossings


def segments_intersect(p1, p2, p3, p4):
    """Check if line segments (p1,p2) and (p3,p4) intersect."""
    d1 = direction(p3, p4, p1)
    d2 = direction(p3, p4, p2)
    d3 = direction(p1, p2, p3)
    d4 = direction(p1, p2, p4)

    if ((d1 > 0 and d2 < 0) or (d1 < 0 and d2 > 0)) and \
       ((d3 > 0 and d4 < 0) or (d3 < 0 and d4 > 0)):
        return True
    return False


def direction(pi, pj, pk):
    return (pk[0] - pi[0]) * (pj[1] - pi[1]) - (pj[0] - pi[0]) * (pk[1] - pi[1])


# ===========================================================================
# Diagnostic 5: Test with relaxed parameters
# ===========================================================================
def diagnose_parameter_sensitivity(boundary, title=""):
    """Test crypt detection with different parameter settings."""
    smooth = direct_fourier_reconstruct(boundary, n_harmonics=25, n_points=1000)
    curv = compute_curvature(smooth)
    norms = compute_normals(smooth)

    param_sets = [
        ("Paper (Montes-Olivas)", PAPER_PARAMS.copy()),
        ("Simulation calibrated", SIMULATION_PARAMS.copy()),
        ("Relaxed area", {**SIMULATION_PARAMS, 'min_area': 0.001, 'max_area': 0.5}),
        ("Tight area", {**SIMULATION_PARAMS, 'min_area': 0.005, 'min_arc_length': 0.02}),
    ]

    fig, axes = plt.subplots(1, len(param_sets), figsize=(6 * len(param_sets), 6))

    smooth_closed = np.vstack([smooth, smooth[0:1]])
    total_area = compute_section_area(smooth, list(range(len(smooth))))
    total_arc = compute_arc_length(smooth_closed)

    for idx, (label, params) in enumerate(param_sets):
        ax = axes[idx]
        sections, inside_mask = find_crypt_sections(smooth, curv, norms, params, return_mask=True)
        valid = filter_crypts(smooth, sections, params)

        ax.plot(smooth_closed[:, 0], smooth_closed[:, 1], 'k-', lw=1)
        ax.fill(smooth_closed[:, 0], smooth_closed[:, 1], color='lightgray', alpha=0.3)

        # Show all sections in light color
        for sec in sections:
            pts = smooth[sec]
            ax.plot(pts[:, 0], pts[:, 1], color='orange', lw=2, alpha=0.4)

        # Highlight valid
        for i, crypt in enumerate(valid):
            pts = smooth[crypt['indices']]
            ax.fill(pts[:, 0], pts[:, 1], color='lime', alpha=0.5)
            ax.plot(pts[:, 0], pts[:, 1], color='darkgreen', lw=2)

        # Section metrics summary
        info_lines = []
        for i, sec in enumerate(sections):
            sa = compute_section_area(smooth, sec) / total_area if total_area > 0 else 0
            sl = compute_arc_length(smooth[sec]) / total_arc if total_arc > 0 else 0
            pass_a = params['min_area'] <= sa < params['max_area']
            pass_l = sl >= params['min_arc_length']
            marker = "✓" if (pass_a and pass_l) else "✗"
            info_lines.append(f'{marker} S{i+1}: A={sa:.4f} L={sl:.4f}')

        info = '\n'.join(info_lines[:12])
        ax.text(0.02, 0.98, info, transform=ax.transAxes, fontsize=6,
                va='top', ha='left', family='monospace',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))

        ax.set_title(f'{label}\n{len(sections)} sections → {len(valid)} crypts')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)

    fig.suptitle(f'{title} — Parameter Sensitivity (Direct Fourier)', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    path = os.path.join(OUTPUT_DIR, f'{title.replace(" ", "_")}_params.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")


# ===========================================================================
# Main
# ===========================================================================
def main():
    print("=" * 60)
    print("Crypt Counting Debug Analysis")
    print("=" * 60)

    # --- Node2d: pick a low-stiffness (should be circular) and high-stiffness (should have budding) ---
    print("\n--- Node2d Diagnostics ---")

    test_cases_node2d = [
        (0.5, 0, "node2d_s0.5_r0"),   # Low stiffness, expect circular
        (100.0, 1, "node2d_s100_r1"),  # High stiffness, expect budding (circularity=0.874)
        (50.0, 0, "node2d_s50_r0"),    # Medium-high stiffness
    ]

    for stiffness, run, label in test_cases_node2d:
        print(f"\n  [{label}] stiffness={stiffness}, run={run}")
        try:
            boundary, cell_types = load_node2d_outline(stiffness, run)
            print(f"    Loaded boundary: {len(boundary)} points, range x=[{boundary[:,0].min():.2f}, {boundary[:,0].max():.2f}], y=[{boundary[:,1].min():.2f}, {boundary[:,1].max():.2f}]")

            diagnose_chain_code(boundary, label)
            diagnose_fourier(boundary, label)
            diagnose_full_pipeline(boundary, label)
            diagnose_parameter_sensitivity(boundary, label)

        except Exception as e:
            print(f"    ERROR: {e}")
            import traceback; traceback.print_exc()

    # --- Vertex2d: diagnose boundary extraction ---
    print("\n--- Vertex2d Diagnostics ---")

    test_cases_vertex2d = [
        (0.5, 0, "vertex2d_s0.5_r0"),
        (100.0, 0, "vertex2d_s100_r0"),
    ]

    for stiffness, run, label in test_cases_vertex2d:
        print(f"\n  [{label}] stiffness={stiffness}, run={run}")
        try:
            boundary, positions = diagnose_vertex2d_boundary(stiffness, run, label)

            if boundary is not None and len(boundary) > 10:
                diagnose_full_pipeline(boundary, label)
                diagnose_parameter_sensitivity(boundary, label)

            # Also try with alpha shape from centroids
            if positions is not None and len(positions) > 10:
                from simple_crypt_count import extract_boundary_alpha_shape
                alpha_boundary = extract_boundary_alpha_shape(positions)
                if len(alpha_boundary) > 10:
                    diagnose_full_pipeline(alpha_boundary, f"{label}_alpha")
                    diagnose_parameter_sensitivity(alpha_boundary, f"{label}_alpha")

        except Exception as e:
            print(f"    ERROR: {e}")
            import traceback; traceback.print_exc()

    # --- End-to-end verification: run count_crypts_simple_method with new defaults ---
    print("\n--- End-to-End Verification (SIMULATION_PARAMS + direct Fourier) ---")
    print(f"Parameters: {SIMULATION_PARAMS}")

    for model, stiffness_runs in [
        ("node2d", [(0.5, 0), (0.5, 5), (1.0, 0), (10.0, 0), (50.0, 0), (100.0, 0), (100.0, 1)]),
    ]:
        print(f"\n  {model}:")
        for stiffness, run in stiffness_runs:
            try:
                boundary, _ = load_node2d_outline(stiffness, run)
                result = count_crypts_simple_method(boundary, boundary_is_ordered=True)
                print(f"    s={stiffness:>5}, r={run}: {result.num_crypts:2d} crypts, circ={result.circularity:.4f}")
            except Exception as e:
                print(f"    s={stiffness:>5}, r={run}: ERROR {e}")

    print(f"\n  vertex2d (t=35h):")
    for stiffness, run in [(0.5, 0), (0.5, 5), (10.0, 0), (50.0, 0), (100.0, 0)]:
        try:
            boundary, _ = load_vertex2d_boundary(stiffness, run, target_time=35.0)
            result = count_crypts_simple_method(boundary, boundary_is_ordered=True)
            n_cross = count_self_intersections(boundary) if len(boundary) < 500 else -1
            print(f"    s={stiffness:>5}, r={run}: {result.num_crypts:2d} crypts, circ={result.circularity:.4f}, boundary_pts={len(boundary)}, crossings={n_cross}")
        except Exception as e:
            print(f"    s={stiffness:>5}, r={run}: ERROR {e}")

    # Generate a final "fixed" debug image for key cases
    print("\n--- Generating fixed debug images ---")
    from simple_crypt_count import plot_crypt_outline

    for stiffness, run in [(0.5, 0), (100.0, 1)]:
        boundary, _ = load_node2d_outline(stiffness, run)
        result = count_crypts_simple_method(boundary, boundary_is_ordered=True)
        path = os.path.join(OUTPUT_DIR, f'FIXED_node2d_s{stiffness}_r{run}.png')
        title = f'FIXED node2d - ECM Stiffness = {stiffness} - run_{run}'
        plot_crypt_outline(result, output_path=path, title=title)

    for stiffness, run in [(0.5, 0), (100.0, 0)]:
        try:
            boundary, _ = load_vertex2d_boundary(stiffness, run, target_time=35.0)
            result = count_crypts_simple_method(boundary, boundary_is_ordered=True)
            path = os.path.join(OUTPUT_DIR, f'FIXED_vertex2d_s{stiffness}_r{run}_t35.png')
            title = f'FIXED vertex2d - ECM Stiffness = {stiffness} - run_{run} (t=35h)'
            plot_crypt_outline(result, output_path=path, title=title)
        except Exception as e:
            print(f"  vertex2d s={stiffness} r={run}: ERROR {e}")

    print(f"\n{'=' * 60}")
    print(f"Debug images saved to: {OUTPUT_DIR}")
    print(f"{'=' * 60}")


if __name__ == '__main__':
    main()
