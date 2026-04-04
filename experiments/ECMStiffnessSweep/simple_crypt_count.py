#!/usr/bin/env python3
"""SimpleCryptCount algorithm for Chaste cell-based simulation outputs."""

import numpy as np
import os
import glob
import xml.etree.ElementTree as ET
from collections import namedtuple, defaultdict

from scipy.spatial import ConvexHull
from scipy.ndimage import uniform_filter1d
from matplotlib.path import Path
import matplotlib.pyplot as plt
import matplotlib.cm as cm

PAPER_PARAMS = {
    'fourier_harmonics': 25,
    'n_points': 1000,
    'min_area': 0.0666,
    'max_area': 0.2736,
    'min_arc_length': 0.1466,
    'circularity_threshold': 0.98
}

SIMULATION_PARAMS = {
    'fourier_harmonics': 24,
    'n_points': 1000,
    'min_area': 0.0068,
    'max_area': 0.0474,
    'min_arc_length': 0.0357,
    'circularity_threshold': 0.98
}

DEFAULT_PARAMS = SIMULATION_PARAMS.copy()


def load_outline_from_vtp(vtp_path):
    """Load ordered boundary points from a Chaste VTP outline file."""
    tree = ET.parse(vtp_path)
    root = tree.getroot()

    for piece in root.iter('Piece'):
        for points in piece.iter('Points'):
            for da in points.iter('DataArray'):
                text = da.text.strip()
                vals = [float(v) for v in text.split()]
                n_points = len(vals) // 3
                boundary = np.array([[vals[3*i], vals[3*i + 1]]
                                     for i in range(n_points)])

        cell_types = None
        for pdata in piece.iter('PointData'):
            for da in pdata.iter('DataArray'):
                if da.get('Name') == 'cell_type':
                    text = da.text.strip()
                    cell_types = np.array([float(v) for v in text.split()])

    return boundary, cell_types


def load_final_outline(data_dir):
    """Load the final (latest timestep) outline from a simulation directory."""
    vtp_files = sorted(glob.glob(os.path.join(data_dir, 'outline_*.vtp')))
    if not vtp_files:
        raise FileNotFoundError(f"No outline_*.vtp files in {data_dir}")

    def get_timestep(path):
        return int(os.path.basename(path).replace('outline_', '').replace('.vtp', ''))

    vtp_files.sort(key=get_timestep)
    return load_outline_from_vtp(vtp_files[-1])


def load_outline_at_time(data_dir, target_time, dt):
    """Load the outline VTP file closest to the given simulation time."""
    vtp_files = sorted(glob.glob(os.path.join(data_dir, 'outline_*.vtp')))
    if not vtp_files:
        raise FileNotFoundError(f"No outline_*.vtp files in {data_dir}")

    target_ts = target_time / dt

    def get_timestep(path):
        return int(os.path.basename(path).replace('outline_', '').replace('.vtp', ''))

    best = min(vtp_files, key=lambda p: abs(get_timestep(p) - target_ts))
    return load_outline_from_vtp(best)


def load_boundary_from_vtu(vtu_path):
    """Load ordered outer boundary from a Chaste vertex mesh VTU file."""
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_path)
    reader.Update()

    mesh = reader.GetOutput()
    vtk_points = mesh.GetPoints()
    points_array = vtk_to_numpy(vtk_points.GetData())
    points = points_array[:, :2]

    n_cells = mesh.GetNumberOfCells()
    cells = []
    for i in range(n_cells):
        cell = mesh.GetCell(i)
        n_pts = cell.GetNumberOfPoints()
        cell_verts = [cell.GetPointId(j) for j in range(n_pts)]
        cells.append(np.array(cell_verts))

    cell_type_data = None
    cell_data = mesh.GetCellData()
    for i in range(cell_data.GetNumberOfArrays()):
        name = cell_data.GetArrayName(i)
        if name and 'cell_type' in name.lower():
            cell_type_data = vtk_to_numpy(cell_data.GetArray(i))
            break

    boundary_pts, boundary_order = _extract_outer_boundary(points, cells)

    boundary_cell_types = None
    if cell_type_data is not None:
        boundary_cell_types = _get_boundary_cell_types(
            boundary_pts, points, cells, cell_type_data)

    return boundary_pts, boundary_cell_types


def _extract_outer_boundary(points, cells):
    """Extract the outer boundary from a vertex mesh using directed edge walking.

    Selects the loop enclosing the largest area (basal boundary).
    Falls back to angular sorting if directed walking fails.
    """
    directed_edges = set()
    for cell_verts in cells:
        n = len(cell_verts)
        for i in range(n):
            v1, v2 = int(cell_verts[i]), int(cell_verts[(i + 1) % n])
            directed_edges.add((v1, v2))

    boundary_directed = []
    for (v1, v2) in directed_edges:
        if (v2, v1) not in directed_edges:
            boundary_directed.append((v1, v2))

    if not boundary_directed:
        raise ValueError("No boundary edges found - mesh may not have an outer boundary")

    next_map = defaultdict(list)
    for v1, v2 in boundary_directed:
        next_map[v1].append(v2)

    visited_edges = set()
    loops = []

    for start_v1, start_v2 in boundary_directed:
        if (start_v1, start_v2) in visited_edges:
            continue

        loop = [start_v1]
        current = start_v1

        while True:
            candidates = next_map.get(current, [])
            chosen = None
            for nxt in candidates:
                if (current, nxt) not in visited_edges:
                    chosen = nxt
                    break

            if chosen is None:
                break

            visited_edges.add((current, chosen))
            if chosen == loop[0] and len(loop) > 2:
                break

            loop.append(chosen)
            current = chosen

        if len(loop) > 2:
            loops.append(loop)

    if not loops:
        boundary_verts = _angular_sort_boundary(points, boundary_directed)
        return points[boundary_verts], boundary_verts

    def _loop_enclosed_area(loop):
        pts = points[np.array(loop)]
        x, y = pts[:, 0], pts[:, 1]
        return 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))

    outer = max(loops, key=_loop_enclosed_area)
    ordered_verts = np.array(outer)
    return points[ordered_verts], ordered_verts


def _angular_sort_boundary(points, boundary_directed):
    """Fallback: sort boundary vertices by angle from centroid."""
    boundary_vert_set = set()
    for v1, v2 in boundary_directed:
        boundary_vert_set.add(v1)
        boundary_vert_set.add(v2)

    boundary_verts = np.array(sorted(boundary_vert_set))
    boundary_pts = points[boundary_verts]
    centroid = boundary_pts.mean(axis=0)
    angles = np.arctan2(boundary_pts[:, 1] - centroid[1],
                        boundary_pts[:, 0] - centroid[0])
    order = np.argsort(angles)
    return boundary_verts[order]


def _get_boundary_cell_types(boundary_pts, all_points, cells, cell_type_data):
    """Get cell type for each boundary point by finding nearest cell."""
    boundary_cell_types = np.zeros(len(boundary_pts))

    for i, pt in enumerate(boundary_pts):
        min_dist = float('inf')
        best_cell = 0
        for cell_idx, cell_verts in enumerate(cells):
            cell_pts = all_points[cell_verts]
            dists = np.linalg.norm(cell_pts - pt, axis=1)
            if dists.min() < min_dist:
                min_dist = dists.min()
                best_cell = cell_idx
        boundary_cell_types[i] = cell_type_data[best_cell] if cell_type_data is not None else 0

    return boundary_cell_types


def load_final_vertex_boundary(data_dir, target_time=None, dt=None):
    """Load boundary from a vertex mesh simulation."""
    vtu_files = [f for f in glob.glob(os.path.join(data_dir, 'results_*.vtu'))
                 if 'ecm_grid' not in f]
    if not vtu_files:
        raise FileNotFoundError(f"No results_*.vtu files in {data_dir}")

    def get_timestep(path):
        return int(os.path.basename(path).replace('results_', '').replace('.vtu', ''))

    if target_time is not None and dt is not None:
        target_ts = target_time / dt
        chosen_vtu = min(vtu_files, key=lambda p: abs(get_timestep(p) - target_ts))
    else:
        vtu_files.sort(key=get_timestep)
        chosen_vtu = vtu_files[-1]

    print(f"Loading vertex boundary from: {os.path.basename(chosen_vtu)}")
    return load_boundary_from_vtu(chosen_vtu)


def extract_boundary_convex_hull(positions):
    """Extract ordered boundary points using convex hull."""
    if len(positions) < 3:
        return positions
    hull = ConvexHull(positions)
    return positions[hull.vertices]


def extract_boundary_alpha_shape(positions, alpha=None):
    """Extract boundary using alpha shape (concave hull)."""
    from scipy.spatial import Delaunay

    if len(positions) < 4:
        return positions

    tri = Delaunay(positions)

    if alpha is None:
        edges = []
        for simplex in tri.simplices:
            for i in range(3):
                for j in range(i + 1, 3):
                    edge_len = np.linalg.norm(positions[simplex[i]] - positions[simplex[j]])
                    edges.append(edge_len)
        alpha = np.mean(edges) * 2.0

    edge_count = {}
    for simplex in tri.simplices:
        pts = positions[simplex]
        a = np.linalg.norm(pts[0] - pts[1])
        b = np.linalg.norm(pts[1] - pts[2])
        c = np.linalg.norm(pts[2] - pts[0])
        s = (a + b + c) / 2
        area = np.sqrt(max(0, s * (s-a) * (s-b) * (s-c)))
        if area > 1e-10:
            circumradius = (a * b * c) / (4 * area)
        else:
            circumradius = float('inf')

        if circumradius < alpha:
            for i in range(3):
                edge = tuple(sorted([simplex[i], simplex[(i+1) % 3]]))
                edge_count[edge] = edge_count.get(edge, 0) + 1

    boundary_edges = [e for e, count in edge_count.items() if count == 1]

    if not boundary_edges:
        return extract_boundary_convex_hull(positions)

    return order_boundary_edges(positions, boundary_edges)


def order_boundary_edges(positions, edges):
    """Order boundary edges into a continuous path."""
    if not edges:
        return np.array([])

    adj = {}
    for e in edges:
        adj.setdefault(e[0], []).append(e[1])
        adj.setdefault(e[1], []).append(e[0])

    ordered = [edges[0][0]]
    visited = {edges[0][0]}

    while True:
        current = ordered[-1]
        neighbors = adj.get(current, [])
        next_point = None
        for n in neighbors:
            if n not in visited:
                next_point = n
                break
        if next_point is None:
            break
        ordered.append(next_point)
        visited.add(next_point)

    return positions[ordered]


def direct_fourier_reconstruct(boundary, n_harmonics=25, n_points=1000):
    """Fourier reconstruction working directly on continuous boundary coordinates,
    parameterized by cumulative arc length. Caps harmonics at Nyquist limit.
    """
    n = len(boundary)
    if n < 4:
        return boundary

    n_harmonics = min(n_harmonics, n // 2)
    if n_harmonics < 1:
        return boundary

    closed = np.vstack([boundary, boundary[0:1]])
    diffs = np.diff(closed, axis=0)
    seg_lengths = np.linalg.norm(diffs, axis=1)
    cum_length = np.concatenate([[0], np.cumsum(seg_lengths)])
    total_length = cum_length[-1]

    if total_length < 1e-10:
        return boundary

    t = cum_length[:-1] / total_length * 2 * np.pi

    x = boundary[:, 0]
    y = boundary[:, 1]

    a0_x = np.mean(x)
    a0_y = np.mean(y)

    a_x = np.zeros(n_harmonics)
    b_x = np.zeros(n_harmonics)
    a_y = np.zeros(n_harmonics)
    b_y = np.zeros(n_harmonics)

    for k in range(1, n_harmonics + 1):
        cos_kt = np.cos(k * t)
        sin_kt = np.sin(k * t)
        a_x[k-1] = 2 * np.mean(x * cos_kt)
        b_x[k-1] = 2 * np.mean(x * sin_kt)
        a_y[k-1] = 2 * np.mean(y * cos_kt)
        b_y[k-1] = 2 * np.mean(y * sin_kt)

    t_new = np.linspace(0, 2 * np.pi, n_points, endpoint=False)

    x_new = np.full(n_points, a0_x)
    y_new = np.full(n_points, a0_y)

    for k in range(1, n_harmonics + 1):
        cos_kt = np.cos(k * t_new)
        sin_kt = np.sin(k * t_new)
        x_new += a_x[k-1] * cos_kt + b_x[k-1] * sin_kt
        y_new += a_y[k-1] * cos_kt + b_y[k-1] * sin_kt

    return np.column_stack([x_new, y_new])


def compute_curvature(boundary):
    """Compute signed curvature at each point along a 2D boundary."""
    n = len(boundary)
    if n < 3:
        return np.zeros(n)

    prev_idx = np.roll(np.arange(n), 1)
    next_idx = np.roll(np.arange(n), -1)

    p_prev = boundary[prev_idx]
    p_curr = boundary
    p_next = boundary[next_idx]

    v1 = p_curr - p_prev
    v2 = p_next - p_curr

    d1 = np.maximum(np.linalg.norm(v1, axis=1), 1e-10)
    d2 = np.maximum(np.linalg.norm(v2, axis=1), 1e-10)

    cross = v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0]
    denom = np.maximum(d1 * d2 * np.linalg.norm(v1 + v2, axis=1), 1e-10)

    return 2 * cross / denom


def compute_normals(boundary):
    """Compute outward-pointing unit normals at each boundary point."""
    n = len(boundary)
    if n < 2:
        return np.zeros((n, 2))

    next_idx = np.roll(np.arange(n), -1)
    prev_idx = np.roll(np.arange(n), 1)

    tangent = boundary[next_idx] - boundary[prev_idx]
    tangent_norm = np.maximum(np.linalg.norm(tangent, axis=1, keepdims=True), 1e-10)
    tangent = tangent / tangent_norm

    normals = np.column_stack([-tangent[:, 1], tangent[:, 0]])

    centroid = boundary.mean(axis=0)
    to_centroid = centroid - boundary
    dots = np.sum(normals * to_centroid, axis=1)
    normals[dots > 0] *= -1

    return normals


def find_crypt_sections(boundary, curvature, normals, params=None, return_mask=False):
    """Identify crypt sections where curvature-scaled normals point inward."""
    if params is None:
        params = DEFAULT_PARAMS

    n = len(boundary)
    if n < 10:
        return ([], np.array([])) if return_mask else []

    scaled_normal_endpoints = boundary + curvature[:, np.newaxis] * normals
    polygon_path = Path(boundary)
    inside = polygon_path.contains_points(scaled_normal_endpoints)

    sections = []
    current_section = []

    for i in range(n):
        if inside[i]:
            current_section.append(i)
        else:
            if len(current_section) >= 3:
                sections.append(current_section)
            current_section = []

    if len(current_section) >= 3:
        if sections and sections[0][0] == 0:
            sections[0] = current_section + sections[0]
        else:
            sections.append(current_section)

    if return_mask:
        return sections, inside
    return sections


def compute_arc_length(points):
    """Compute total arc length of a path."""
    if len(points) < 2:
        return 0
    diff = np.diff(points, axis=0)
    return np.sum(np.linalg.norm(diff, axis=1))


def compute_section_area(boundary, section_indices):
    """Compute area of a section polygon (shoelace formula)."""
    if len(section_indices) < 3:
        return 0

    section_points = boundary[section_indices]
    n = len(section_points)
    area = 0
    for i in range(n):
        j = (i + 1) % n
        area += section_points[i, 0] * section_points[j, 1]
        area -= section_points[j, 0] * section_points[i, 1]
    return abs(area) / 2


def filter_crypts(boundary, sections, params=None):
    """Filter detected crypt sections by area and arc length criteria."""
    if params is None:
        params = DEFAULT_PARAMS

    total_area = compute_section_area(boundary, list(range(len(boundary))))
    total_arc_length = compute_arc_length(np.vstack([boundary, boundary[0:1]]))

    if total_area < 1e-10 or total_arc_length < 1e-10:
        return []

    valid_crypts = []

    for section in sections:
        section_area = compute_section_area(boundary, section)
        norm_area = section_area / total_area

        section_points = boundary[section]
        section_arc_length = compute_arc_length(section_points)
        norm_arc_length = section_arc_length / total_arc_length

        if (params['min_area'] <= norm_area < params['max_area'] and
            norm_arc_length >= params['min_arc_length']):
            valid_crypts.append({
                'indices': section,
                'norm_area': norm_area,
                'norm_arc_length': norm_arc_length
            })

    return valid_crypts


def compute_circularity(boundary):
    """Compute circularity = 4*pi*Area/Perimeter^2. Returns 1.0 for a perfect circle."""
    if len(boundary) < 3:
        return np.nan

    n = len(boundary)
    area = 0
    for i in range(n):
        j = (i + 1) % n
        area += boundary[i, 0] * boundary[j, 1]
        area -= boundary[j, 0] * boundary[i, 1]
    area = abs(area) / 2

    perimeter = compute_arc_length(np.vstack([boundary, boundary[0:1]]))

    if perimeter < 1e-10:
        return np.nan

    return 4 * np.pi * area / (perimeter ** 2)


CryptCountResult = namedtuple('CryptCountResult', [
    'num_crypts', 'circularity', 'crypts', 'boundary', 'smooth_boundary',
    'curvature', 'normals', 'all_sections', 'inside_mask'
])


def count_crypts_simple_method(positions, params=None, use_alpha_shape=True,
                                boundary_is_ordered=False):
    """Count crypts using the SimpleCryptCount method.

    Args:
        positions: Cell positions (N x 2) or pre-ordered boundary points
        params: Detection parameters (see DEFAULT_PARAMS)
        use_alpha_shape: Use alpha shape (True) or convex hull for boundary
        boundary_is_ordered: If True, treat positions as ordered boundary

    Returns:
        CryptCountResult namedtuple
    """
    if params is None:
        params = DEFAULT_PARAMS.copy()

    if len(positions) < 10:
        return CryptCountResult(0, np.nan, [], np.array([]), np.array([]),
                                np.array([]), np.array([]), [], np.array([]))

    if boundary_is_ordered:
        boundary = np.asarray(positions)
    elif use_alpha_shape:
        boundary = extract_boundary_alpha_shape(positions)
    else:
        boundary = extract_boundary_convex_hull(positions)

    if len(boundary) < 10:
        return CryptCountResult(0, np.nan, [], boundary, np.array([]),
                                np.array([]), np.array([]), [], np.array([]))

    smooth_boundary = direct_fourier_reconstruct(
        boundary,
        n_harmonics=params.get('fourier_harmonics', 25),
        n_points=params.get('n_points', 1000),
    )

    circularity = compute_circularity(boundary)

    curvature = compute_curvature(smooth_boundary)
    normals = compute_normals(smooth_boundary)

    sections, inside_mask = find_crypt_sections(smooth_boundary, curvature, normals,
                                                 params, return_mask=True)

    if circularity > params.get('circularity_threshold', 0.98):
        return CryptCountResult(0, circularity, [], boundary, smooth_boundary,
                                curvature, normals, sections, inside_mask)

    valid_crypts = filter_crypts(smooth_boundary, sections, params)

    return CryptCountResult(
        num_crypts=len(valid_crypts),
        circularity=circularity,
        crypts=valid_crypts,
        boundary=boundary,
        smooth_boundary=smooth_boundary,
        curvature=curvature,
        normals=normals,
        all_sections=sections,
        inside_mask=inside_mask
    )


def plot_crypt_analysis(result, positions=None, output_path=None):
    """Visualize the crypt detection results (3-panel)."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    ax1 = axes[0]
    if positions is not None:
        ax1.scatter(positions[:, 0], positions[:, 1], s=5, alpha=0.5, label='Cells')
    if len(result.boundary) > 0:
        boundary_closed = np.vstack([result.boundary, result.boundary[0:1]])
        ax1.plot(boundary_closed[:, 0], boundary_closed[:, 1], 'b-', lw=2, label='Boundary')
    ax1.set_title('Cell Positions & Boundary')
    ax1.set_aspect('equal')
    ax1.legend()

    ax2 = axes[1]
    if len(result.smooth_boundary) > 0:
        smooth_closed = np.vstack([result.smooth_boundary, result.smooth_boundary[0:1]])
        ax2.plot(smooth_closed[:, 0], smooth_closed[:, 1], 'k-', lw=1)
        sc = ax2.scatter(result.smooth_boundary[:, 0], result.smooth_boundary[:, 1],
                        c=result.curvature, cmap='RdBu', s=2)
        plt.colorbar(sc, ax=ax2, label='Curvature')
    ax2.set_title(f'Smooth Boundary (Circularity: {result.circularity:.3f})')
    ax2.set_aspect('equal')

    ax3 = axes[2]
    if len(result.smooth_boundary) > 0:
        smooth_closed = np.vstack([result.smooth_boundary, result.smooth_boundary[0:1]])
        ax3.plot(smooth_closed[:, 0], smooth_closed[:, 1], 'k-', lw=1)
        colors = plt.cm.Set1(np.linspace(0, 1, max(len(result.crypts), 1)))
        for i, crypt in enumerate(result.crypts):
            indices = crypt['indices']
            ax3.plot(result.smooth_boundary[indices, 0],
                    result.smooth_boundary[indices, 1],
                    color=colors[i], lw=3)
    ax3.set_title(f'Detected Crypts: {result.num_crypts}')
    ax3.set_aspect('equal')

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_path}")
    else:
        plt.show()
    plt.close()


def plot_debug_analysis(result, positions=None, cell_types=None, output_path=None,
                        params=None, title_prefix=''):
    """Comprehensive 2x3 debug visualization of each SimpleCryptCount step."""
    if params is None:
        params = DEFAULT_PARAMS

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    # Panel 1: Original boundary
    ax1 = axes[0, 0]
    if positions is not None and not np.array_equal(positions, result.boundary):
        ax1.scatter(positions[:, 0], positions[:, 1], s=10, alpha=0.3,
                   c='gray', label='Cell positions')

    if len(result.boundary) > 0:
        boundary_closed = np.vstack([result.boundary, result.boundary[0:1]])
        if cell_types is not None:
            colors = ['blue' if ct == 0 else 'red' for ct in cell_types]
            for i in range(len(result.boundary)):
                ax1.plot([result.boundary[i, 0], result.boundary[(i+1) % len(result.boundary), 0]],
                        [result.boundary[i, 1], result.boundary[(i+1) % len(result.boundary), 1]],
                        c=colors[i], lw=2)
            ax1.scatter([], [], c='blue', s=50, label='Stem cells (type 0)')
            ax1.scatter([], [], c='red', s=50, label='Other cells (type 1)')
        else:
            ax1.plot(boundary_closed[:, 0], boundary_closed[:, 1], 'b-', lw=2,
                    label='Input boundary')
        ax1.scatter(result.boundary[0, 0], result.boundary[0, 1],
                   s=100, c='green', marker='*', zorder=5, label='Start point')

    ax1.set_title(f'1. Input Boundary\n({len(result.boundary)} points)')
    ax1.set_aspect('equal')
    ax1.legend(loc='upper right', fontsize=8)
    ax1.grid(True, alpha=0.3)

    # Panel 2: Fourier-smoothed boundary
    ax2 = axes[0, 1]
    if len(result.boundary) > 0:
        boundary_closed = np.vstack([result.boundary, result.boundary[0:1]])
        ax2.plot(boundary_closed[:, 0], boundary_closed[:, 1], 'b--', lw=1,
                alpha=0.5, label='Original')
    if len(result.smooth_boundary) > 0:
        smooth_closed = np.vstack([result.smooth_boundary, result.smooth_boundary[0:1]])
        ax2.plot(smooth_closed[:, 0], smooth_closed[:, 1], 'r-', lw=2,
                label=f'Fourier ({params.get("fourier_harmonics", 25)} harmonics)')
    ax2.set_title(f'2. Fourier Smoothing\n(Circularity: {result.circularity:.4f})')
    ax2.set_aspect('equal')
    ax2.legend(loc='upper right', fontsize=8)
    ax2.grid(True, alpha=0.3)

    # Panel 3: Curvature along boundary
    ax3 = axes[0, 2]
    if len(result.smooth_boundary) > 0 and len(result.curvature) > 0:
        smooth_closed = np.vstack([result.smooth_boundary, result.smooth_boundary[0:1]])
        ax3.plot(smooth_closed[:, 0], smooth_closed[:, 1], 'k-', lw=0.5, alpha=0.3)
        vmax = max(abs(result.curvature.min()), abs(result.curvature.max()))
        sc = ax3.scatter(result.smooth_boundary[:, 0], result.smooth_boundary[:, 1],
                        c=result.curvature, cmap='RdBu_r', s=5, vmin=-vmax, vmax=vmax)
        cbar = plt.colorbar(sc, ax=ax3)
        cbar.set_label('Curvature\n(+ convex, - concave)')
        ax3.text(0.02, 0.98, f'Min: {result.curvature.min():.3f}\nMax: {result.curvature.max():.3f}',
                transform=ax3.transAxes, fontsize=8, va='top', ha='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax3.set_title('3. Curvature (LineCurvature2D)')
    ax3.set_aspect('equal')
    ax3.grid(True, alpha=0.3)

    # Panel 4: Curvature-scaled normals
    ax4 = axes[1, 0]
    if len(result.smooth_boundary) > 0 and len(result.normals) > 0:
        smooth_closed = np.vstack([result.smooth_boundary, result.smooth_boundary[0:1]])
        ax4.plot(smooth_closed[:, 0], smooth_closed[:, 1], 'k-', lw=1)
        step = max(1, len(result.smooth_boundary) // 50)
        for i in range(0, len(result.smooth_boundary), step):
            pt = result.smooth_boundary[i]
            normal_end = pt + result.curvature[i] * result.normals[i]
            is_inside = result.inside_mask[i] if len(result.inside_mask) > 0 else False
            color = 'green' if is_inside else 'red'
            ax4.plot([pt[0], normal_end[0]], [pt[1], normal_end[1]],
                    color=color, lw=1, alpha=0.7)
            ax4.scatter(normal_end[0], normal_end[1], c=color, s=10, alpha=0.7)
        ax4.scatter([], [], c='green', s=30, label='Normal tip INSIDE (crypt candidate)')
        ax4.scatter([], [], c='red', s=30, label='Normal tip OUTSIDE (organoid body)')
        n_inside = np.sum(result.inside_mask) if len(result.inside_mask) > 0 else 0
        n_outside = len(result.inside_mask) - n_inside if len(result.inside_mask) > 0 else 0
        ax4.text(0.02, 0.98, f'Inside: {n_inside}\nOutside: {n_outside}',
                transform=ax4.transAxes, fontsize=8, va='top', ha='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax4.set_title('4. Curvature-Scaled Normals\n(inpolygon test)')
    ax4.set_aspect('equal')
    ax4.legend(loc='upper right', fontsize=7)
    ax4.grid(True, alpha=0.3)

    # Panel 5: All detected sections (before filtering)
    ax5 = axes[1, 1]
    if len(result.smooth_boundary) > 0:
        smooth_closed = np.vstack([result.smooth_boundary, result.smooth_boundary[0:1]])
        ax5.plot(smooth_closed[:, 0], smooth_closed[:, 1], 'k-', lw=1, alpha=0.5)
        total_area = compute_section_area(result.smooth_boundary,
                                          list(range(len(result.smooth_boundary))))
        total_arc = compute_arc_length(smooth_closed)
        colors = plt.cm.tab10(np.linspace(0, 1, max(len(result.all_sections), 1)))
        section_info = []
        for i, section in enumerate(result.all_sections):
            section_pts = result.smooth_boundary[section]
            ax5.plot(section_pts[:, 0], section_pts[:, 1],
                    color=colors[i], lw=4, alpha=0.8)
            sec_area = compute_section_area(result.smooth_boundary, section)
            sec_arc = compute_arc_length(section_pts)
            norm_area = sec_area / total_area if total_area > 0 else 0
            norm_arc = sec_arc / total_arc if total_arc > 0 else 0
            mid_idx = len(section) // 2
            mid_pt = section_pts[mid_idx]
            ax5.annotate(f'{i+1}', mid_pt, fontsize=8, ha='center', va='center',
                        bbox=dict(boxstyle='circle', facecolor=colors[i], alpha=0.8))
            section_info.append(f'S{i+1}: A={norm_area:.3f}, L={norm_arc:.3f}')
        info_text = '\n'.join(section_info[:10])
        if len(section_info) > 10:
            info_text += f'\n... +{len(section_info)-10} more'
        ax5.text(0.02, 0.98, info_text, transform=ax5.transAxes, fontsize=7,
                va='top', ha='left', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                family='monospace')
    ax5.set_title(f'5. All Sections (before filtering)\n({len(result.all_sections)} sections)')
    ax5.set_aspect('equal')
    ax5.grid(True, alpha=0.3)

    # Panel 6: Final filtered crypts
    ax6 = axes[1, 2]
    if len(result.smooth_boundary) > 0:
        smooth_closed = np.vstack([result.smooth_boundary, result.smooth_boundary[0:1]])
        ax6.fill(smooth_closed[:, 0], smooth_closed[:, 1],
                color='lightgray', alpha=0.5, label='Organoid body')
        ax6.plot(smooth_closed[:, 0], smooth_closed[:, 1], 'k-', lw=1)
        colors = plt.cm.Set1(np.linspace(0, 1, max(len(result.crypts), 1)))
        for i, crypt in enumerate(result.crypts):
            indices = crypt['indices']
            crypt_pts = result.smooth_boundary[indices]
            ax6.fill(crypt_pts[:, 0], crypt_pts[:, 1], color=colors[i], alpha=0.6)
            ax6.plot(crypt_pts[:, 0], crypt_pts[:, 1], color=colors[i], lw=3)
            mid_idx = len(indices) // 2
            mid_pt = crypt_pts[mid_idx]
            ax6.annotate(f"Crypt {i+1}\nA={crypt['norm_area']:.3f}\nL={crypt['norm_arc_length']:.3f}",
                        mid_pt, fontsize=7, ha='center', va='center',
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
        filter_text = (f"Filters:\n"
                      f"  Area: [{params.get('min_area', 0):.3f}, {params.get('max_area', 1):.3f})\n"
                      f"  Arc: >= {params.get('min_arc_length', 0):.3f}\n"
                      f"  Circ: < {params.get('circularity_threshold', 0.98):.2f}")
        ax6.text(0.02, 0.02, filter_text, transform=ax6.transAxes, fontsize=7,
                va='bottom', ha='left', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                family='monospace')
    ax6.set_title(f'6. Final Crypts (after filtering)\n({result.num_crypts} crypts detected)')
    ax6.set_aspect('equal')
    ax6.grid(True, alpha=0.3)

    fig.suptitle(f'{title_prefix}SimpleCryptCount Analysis\n'
                 f'Circularity={result.circularity:.4f}, '
                 f'Sections={len(result.all_sections)}, '
                 f'Crypts={result.num_crypts}',
                 fontsize=14, fontweight='bold')

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved debug plot: {output_path}")
    else:
        plt.show()
    plt.close()
    return fig


def plot_crypt_outline(result, output_path=None, title=None, show_metrics=True):
    """Single-panel visualization: organoid outline with crypt regions shaded."""
    fig, ax = plt.subplots(figsize=(10, 10))

    if len(result.smooth_boundary) > 0:
        smooth_closed = np.vstack([result.smooth_boundary, result.smooth_boundary[0:1]])
        ax.fill(smooth_closed[:, 0], smooth_closed[:, 1],
                color='#E8E8E8', alpha=0.8, label='Organoid body')
        ax.plot(smooth_closed[:, 0], smooth_closed[:, 1], 'k-', lw=2, label='Boundary')

        if len(result.crypts) > 0:
            crypt_colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4',
                           '#FFEAA7', '#DDA0DD', '#98D8C8', '#F7DC6F']
            centroid = result.smooth_boundary.mean(axis=0)
            for i, crypt in enumerate(result.crypts):
                indices = crypt['indices']
                crypt_pts = result.smooth_boundary[indices]
                color = crypt_colors[i % len(crypt_colors)]
                ax.fill(crypt_pts[:, 0], crypt_pts[:, 1],
                       color=color, alpha=0.6, edgecolor=color, lw=2)
                mid_idx = len(indices) // 2
                mid_pt = crypt_pts[mid_idx]
                ax.annotate(f'{i+1}', mid_pt, fontsize=14, fontweight='bold',
                           ha='center', va='center', color='white',
                           bbox=dict(boxstyle='circle', facecolor=color,
                                    edgecolor='white', alpha=0.9))

    if show_metrics:
        metrics_text = (f"Crypts: {result.num_crypts}\n"
                       f"Circularity: {result.circularity:.3f}")
        ax.text(0.02, 0.98, metrics_text, transform=ax.transAxes,
                fontsize=12, va='top', ha='left',
                bbox=dict(boxstyle='round', facecolor='white',
                         edgecolor='gray', alpha=0.9))

    ax.set_aspect('equal')
    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_title(title or f'Organoid Outline with {result.num_crypts} Detected Crypts',
                fontsize=14, fontweight='bold')
    ax.grid(False)
    ax.set_facecolor('#FAFAFA')

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        print(f"Saved crypt outline plot: {output_path}")
    else:
        plt.show()
    plt.close()
    return fig


if __name__ == '__main__':
    import argparse
    import sys

    sys.path.insert(0, '.')

    parser = argparse.ArgumentParser(
        description='SimpleCryptCount method for Chaste simulation output')
    parser.add_argument('--data-dir', type=str, required=True)
    parser.add_argument('--use-outline', action='store_true')
    parser.add_argument('--use-vertex-mesh', action='store_true')
    parser.add_argument('--output', type=str, default=None)
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--fourier-harmonics', type=int, default=25)
    parser.add_argument('--min-area', type=float, default=0.0666)
    parser.add_argument('--max-area', type=float, default=0.2736)
    parser.add_argument('--min-arc-length', type=float, default=0.1466)

    args = parser.parse_args()

    params = DEFAULT_PARAMS.copy()
    params['fourier_harmonics'] = args.fourier_harmonics
    params['min_area'] = args.min_area
    params['max_area'] = args.max_area
    params['min_arc_length'] = args.min_arc_length

    positions = None
    cell_types = None
    boundary_is_ordered = False

    if args.use_outline:
        boundary, cell_types = load_final_outline(args.data_dir)
        positions = boundary
        boundary_is_ordered = True
        print(f"Loaded outline with {len(boundary)} points from VTP file")
    elif args.use_vertex_mesh:
        boundary, cell_types = load_final_vertex_boundary(args.data_dir)
        positions = boundary
        boundary_is_ordered = True
        print(f"Loaded vertex mesh boundary with {len(boundary)} points from VTU file")
    else:
        from analyse_crypt_budding import load_final_positions
        positions = load_final_positions(args.data_dir, dim=2)
        print(f"Loaded {len(positions)} cell positions")

    print(f"\nRunning SimpleCryptCount analysis...")
    print(f"  Fourier harmonics: {params['fourier_harmonics']}")
    print(f"  Area range: [{params['min_area']:.4f}, {params['max_area']:.4f})")
    print(f"  Min arc length: {params['min_arc_length']:.4f}")

    result = count_crypts_simple_method(positions, params,
                                        boundary_is_ordered=boundary_is_ordered)

    print(f"\nResults:")
    print(f"  Circularity: {result.circularity:.4f}")
    print(f"  All sections detected: {len(result.all_sections)}")
    print(f"  Final crypts (after filtering): {result.num_crypts}")

    for i, crypt in enumerate(result.crypts):
        print(f"    Crypt {i+1}: norm_area={crypt['norm_area']:.4f}, "
              f"norm_arc_length={crypt['norm_arc_length']:.4f}")

    if args.debug or args.output:
        output_path = args.output or 'simple_crypt_debug.png'
        if args.debug:
            plot_debug_analysis(result,
                               positions if not boundary_is_ordered else None,
                               cell_types=cell_types,
                               output_path=output_path,
                               params=params)
        else:
            plot_crypt_analysis(result,
                               positions if not boundary_is_ordered else None,
                               output_path)
