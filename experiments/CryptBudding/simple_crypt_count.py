#!/usr/bin/env python3
"""
simple_crypt_count.py

Python implementation of the SimpleCryptCount algorithm from:
  "In-silico and in-vitro morphometric analysis of intestinal organoids"
  (Montes-Olivas et al.)

This module adapts the MATLAB crypt counting method for use with Chaste
cell-based simulation outputs. The algorithm:
  1. Extracts boundary from cell positions OR loads from VTP outline files
  2. Computes elliptic Fourier descriptors for smooth boundary reconstruction
  3. Calculates curvature and normals along the boundary
  4. Identifies crypts as inward-pointing (concave) regions
  5. Filters by normalized area and arc length thresholds

Usage:
  from simple_crypt_count import count_crypts_simple_method
  n_crypts, circularity, debug_info = count_crypts_simple_method(positions)
  
  # Or load directly from VTP outline:
  boundary = load_outline_from_vtp('/path/to/outline_6000.vtp')
  result = count_crypts_simple_method(boundary, boundary_is_ordered=True)
"""

import numpy as np
import os
import glob
import xml.etree.ElementTree as ET
from collections import namedtuple

try:
    from scipy.spatial import ConvexHull
    from scipy.ndimage import uniform_filter1d
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("WARNING: scipy not available. SimpleCryptCount requires scipy.")

# Default parameters (from SimpleCryptCount paper, Table 2 - In silico)
DEFAULT_PARAMS = {
    'fourier_harmonics': 25,      # Number of Fourier harmonics
    'n_points': 1000,             # Points for reconstructed boundary
    'min_area': 0.0666,           # Minimum normalized crypt area
    'max_area': 0.2736,           # Maximum normalized crypt area
    'min_arc_length': 0.1466,     # Minimum normalized arc length
    'circularity_threshold': 0.98 # Skip crypt counting if too circular
}


# =====================================================================
# VTP Outline file loading
# =====================================================================

def load_outline_from_vtp(vtp_path):
    """
    Load ordered boundary points from a Chaste VTP outline file.
    
    Args:
        vtp_path: Path to outline_*.vtp file
    
    Returns:
        boundary: Ordered boundary points (N x 2 array)
        cell_types: Cell type values at each boundary point (optional)
    """
    tree = ET.parse(vtp_path)
    root = tree.getroot()
    
    # Find Points data
    for piece in root.iter('Piece'):
        for points in piece.iter('Points'):
            for da in points.iter('DataArray'):
                text = da.text.strip()
                vals = [float(v) for v in text.split()]
                n_points = len(vals) // 3
                # Extract x, y coordinates (ignore z)
                boundary = np.array([[vals[3*i], vals[3*i + 1]] 
                                     for i in range(n_points)])
        
        # Try to get cell types
        cell_types = None
        for pdata in piece.iter('PointData'):
            for da in pdata.iter('DataArray'):
                if da.get('Name') == 'cell_type':
                    text = da.text.strip()
                    cell_types = np.array([float(v) for v in text.split()])
    
    return boundary, cell_types


def load_final_outline(data_dir):
    """
    Load the final (latest timestep) outline from a simulation directory.
    
    Args:
        data_dir: Path to results_from_time_* directory
    
    Returns:
        boundary: Ordered boundary points (N x 2)
        cell_types: Cell type at each point (or None)
    """
    # Find all outline VTP files
    vtp_files = sorted(glob.glob(os.path.join(data_dir, 'outline_*.vtp')))
    
    if not vtp_files:
        raise FileNotFoundError(f"No outline_*.vtp files in {data_dir}")
    
    # Get the one with the largest timestep
    def get_timestep(path):
        basename = os.path.basename(path)
        # outline_6000.vtp -> 6000
        try:
            return int(basename.replace('outline_', '').replace('.vtp', ''))
        except ValueError:
            return -1
    
    vtp_files.sort(key=get_timestep)
    latest_vtp = vtp_files[-1]
    
    return load_outline_from_vtp(latest_vtp)


# =====================================================================
# VTU Vertex Mesh loading (for vertex2d model)
# =====================================================================

def load_boundary_from_vtu(vtu_path):
    """
    Load ordered outer boundary from a Chaste vertex mesh VTU file.
    
    For vertex models, the mesh is stored as polygons. The outer boundary
    consists of edges that appear in only one cell (not shared).
    
    Args:
        vtu_path: Path to results_*.vtu file
    
    Returns:
        boundary: Ordered boundary points (N x 2 array)
        cell_types: Cell type values at each boundary point (or None)
    """
    # Try VTK library first (most robust)
    try:
        return _load_vtu_with_vtk(vtu_path)
    except ImportError:
        pass  # Fall back to manual parsing
    except Exception as e:
        print(f"VTK loading failed: {e}, trying manual parsing...")
    
    # Fall back to manual XML parsing
    return _load_vtu_manual(vtu_path)


def _load_vtu_with_vtk(vtu_path):
    """Load VTU file using VTK library."""
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_path)
    reader.Update()
    
    mesh = reader.GetOutput()
    
    # Get points
    vtk_points = mesh.GetPoints()
    points_array = vtk_to_numpy(vtk_points.GetData())
    points = points_array[:, :2]  # x, y only
    
    # Get cells
    n_cells = mesh.GetNumberOfCells()
    cells = []
    for i in range(n_cells):
        cell = mesh.GetCell(i)
        n_pts = cell.GetNumberOfPoints()
        cell_verts = [cell.GetPointId(j) for j in range(n_pts)]
        cells.append(np.array(cell_verts))
    
    # Get cell type data if available
    cell_type_data = None
    cell_data = mesh.GetCellData()
    for i in range(cell_data.GetNumberOfArrays()):
        name = cell_data.GetArrayName(i)
        if name and 'cell_type' in name.lower():
            arr = cell_data.GetArray(i)
            cell_type_data = vtk_to_numpy(arr)
            break
    
    # Extract boundary
    boundary_pts, boundary_order = _extract_outer_boundary(points, cells)
    
    # Get cell types for boundary points
    boundary_cell_types = None
    if cell_type_data is not None:
        boundary_cell_types = _get_boundary_cell_types(boundary_pts, points, cells, cell_type_data)
    
    return boundary_pts, boundary_cell_types


def _load_vtu_manual(vtu_path):
    """Load VTU file using manual XML parsing (fallback)."""
    import base64
    import struct
    
    tree = ET.parse(vtu_path)
    root = tree.getroot()
    
    points = None
    connectivity = None
    offsets = None
    cell_type_data = None
    
    for piece in root.iter('Piece'):
        n_points = int(piece.get('NumberOfPoints', 0))
        n_cells = int(piece.get('NumberOfCells', 0))
        
        # Get Points
        for pts in piece.iter('Points'):
            for da in pts.iter('DataArray'):
                if da.get('format') == 'appended':
                    # Need to decode from appended data
                    pass  # Fall back to ASCII parsing below
                else:
                    text = da.text.strip() if da.text else ""
                    if text:
                        vals = [float(v) for v in text.split()]
                        points = np.array(vals).reshape(-1, 3)[:, :2]  # x, y only
        
        # Get Cells connectivity and offsets
        for cells in piece.iter('Cells'):
            for da in cells.iter('DataArray'):
                name = da.get('Name', '')
                if name == 'connectivity':
                    if da.get('format') != 'appended' and da.text:
                        connectivity = np.array([int(v) for v in da.text.strip().split()])
                elif name == 'offsets':
                    if da.get('format') != 'appended' and da.text:
                        offsets = np.array([int(v) for v in da.text.strip().split()])
        
        # Get cell type data if available
        for cdata in piece.iter('CellData'):
            for da in cdata.iter('DataArray'):
                name = da.get('Name', '')
                if 'cell_type' in name.lower() and da.get('format') != 'appended' and da.text:
                    cell_type_data = np.array([float(v) for v in da.text.strip().split()])
    
    # If data was in appended format, try to parse it
    if points is None or connectivity is None or offsets is None:
        points, connectivity, offsets, cell_type_data = _parse_vtu_appended(root, vtu_path)
    
    if points is None or connectivity is None or offsets is None:
        raise ValueError(f"Could not parse VTU file: {vtu_path}")
    
    # Build cells as lists of vertex indices
    cells = []
    prev_offset = 0
    for offset in offsets:
        cell_verts = connectivity[prev_offset:offset]
        cells.append(cell_verts)
        prev_offset = offset
    
    # Extract boundary edges
    boundary_pts, boundary_order = _extract_outer_boundary(points, cells)
    
    # Get cell types for boundary points if available
    boundary_cell_types = None
    if cell_type_data is not None:
        # For vertex model, we need to map boundary vertices to cells
        # This is approximate - use nearest cell
        boundary_cell_types = _get_boundary_cell_types(boundary_pts, points, cells, cell_type_data)
    
    return boundary_pts, boundary_cell_types


def _parse_vtu_appended(root, vtu_path):
    """Parse VTU with base64 appended data format."""
    import base64
    import struct
    
    # Find AppendedData
    appended_data = None
    for ad in root.iter('AppendedData'):
        encoding = ad.get('encoding', 'base64')
        if encoding != 'base64':
            raise ValueError(f"Unsupported encoding: {encoding}")
        text = ad.text
        if text and text.startswith('_'):
            appended_data = text[1:]  # Skip leading underscore
    
    if not appended_data:
        return None, None, None, None
    
    points = None
    connectivity = None  
    offsets = None
    cell_type_data = None
    
    for piece in root.iter('Piece'):
        n_points = int(piece.get('NumberOfPoints', 0))
        n_cells = int(piece.get('NumberOfCells', 0))
        
        # Parse Points
        for pts in piece.iter('Points'):
            for da in pts.iter('DataArray'):
                if da.get('format') == 'appended':
                    offset = int(da.get('offset', 0))
                    n_components = int(da.get('NumberOfComponents', 3))
                    dtype = da.get('type', 'Float64')
                    
                    data = _decode_appended_array(appended_data, offset, n_points * n_components, dtype)
                    if data is not None:
                        points = data.reshape(-1, n_components)[:, :2]
        
        # Parse Cells
        for cells in piece.iter('Cells'):
            for da in cells.iter('DataArray'):
                name = da.get('Name', '')
                if da.get('format') == 'appended':
                    offset = int(da.get('offset', 0))
                    dtype = da.get('type', 'Int64')
                    
                    if name == 'connectivity':
                        # Connectivity length = total vertices in all cells
                        # We don't know exactly, so decode all and trim later
                        data = _decode_appended_array(appended_data, offset, n_cells * 10, dtype)  # Estimate
                        if data is not None:
                            connectivity = data
                    elif name == 'offsets':
                        data = _decode_appended_array(appended_data, offset, n_cells, dtype)
                        if data is not None:
                            offsets = data
        
        # Parse CellData for cell types
        for cdata in piece.iter('CellData'):
            for da in cdata.iter('DataArray'):
                name = da.get('Name', '')
                if 'cell_type' in name.lower() and da.get('format') == 'appended':
                    offset = int(da.get('offset', 0))
                    dtype = da.get('type', 'Float64')
                    data = _decode_appended_array(appended_data, offset, n_cells, dtype)
                    if data is not None:
                        cell_type_data = data
    
    # Trim connectivity to actual length based on offsets
    if connectivity is not None and offsets is not None:
        connectivity = connectivity[:int(offsets[-1])]
    
    return points, connectivity, offsets, cell_type_data


def _decode_appended_array(appended_data, offset, expected_count, dtype):
    """Decode a base64-encoded array from VTK appended data."""
    import base64
    import struct
    
    # Map VTK types to struct format
    type_map = {
        'Float64': ('d', 8),
        'Float32': ('f', 4),
        'Int64': ('q', 8),
        'Int32': ('i', 4),
        'UInt8': ('B', 1),
    }
    
    if dtype not in type_map:
        return None
    
    fmt, size = type_map[dtype]
    
    # The appended data has a 4-byte header (UInt32) giving the size in bytes
    # Then the actual data follows
    
    # Base64 encodes 3 bytes into 4 characters
    # Offset is in decoded bytes, need to find corresponding position in encoded string
    
    # VTK appended data: each array is prefixed with a 4-byte size header
    # The offset points to the start of THIS array's header
    
    try:
        # Decode the header first (4 bytes = ~6 base64 chars, but need alignment)
        # This is complex because base64 encoding chunks data
        # Simpler: decode everything and slice
        
        decoded = base64.b64decode(appended_data)
        
        # Read the size header at the given offset
        array_size = struct.unpack('<I', decoded[offset:offset+4])[0]
        
        # Read the actual array data
        data_start = offset + 4
        data_end = data_start + array_size
        
        n_elements = array_size // size
        arr = np.frombuffer(decoded[data_start:data_end], dtype=np.dtype(f'<{fmt}'))
        
        return arr
    except Exception as e:
        print(f"Warning: Could not decode appended array: {e}")
        return None


def _extract_outer_boundary(points, cells):
    """
    Extract the outer boundary from a vertex mesh.
    
    An outer boundary edge is one that belongs to only ONE cell.
    
    Args:
        points: Nx2 array of vertex positions
        cells: List of cells, each cell is a list of vertex indices
    
    Returns:
        boundary_pts: Ordered boundary points (Mx2)
        boundary_indices: Indices of boundary vertices
    """
    from collections import defaultdict
    
    # Count how many times each edge appears
    edge_count = defaultdict(int)
    edge_to_cell = {}  # Track which cell each edge belongs to
    
    for cell_idx, cell_verts in enumerate(cells):
        n = len(cell_verts)
        for i in range(n):
            v1, v2 = cell_verts[i], cell_verts[(i + 1) % n]
            # Canonical edge representation (sorted)
            edge = tuple(sorted([v1, v2]))
            edge_count[edge] += 1
            edge_to_cell[edge] = cell_idx
    
    # Boundary edges appear exactly once
    boundary_edges = [edge for edge, count in edge_count.items() if count == 1]
    
    if not boundary_edges:
        raise ValueError("No boundary edges found - mesh may not have an outer boundary")
    
    # Order the boundary edges into a continuous loop
    ordered_verts = _order_boundary_edges(boundary_edges)
    
    # Get the boundary points
    boundary_pts = points[ordered_verts]
    
    return boundary_pts, ordered_verts


def _order_boundary_edges(edges):
    """
    Order boundary edges into a continuous loop.
    
    Args:
        edges: List of (v1, v2) tuples
    
    Returns:
        ordered_verts: List of vertex indices in order around the boundary
    """
    from collections import defaultdict
    
    # Build adjacency: for each vertex, which vertices are connected
    adj = defaultdict(list)
    for v1, v2 in edges:
        adj[v1].append(v2)
        adj[v2].append(v1)
    
    # Start from any boundary vertex
    start = edges[0][0]
    ordered = [start]
    visited = {start}
    
    current = start
    while True:
        # Find next unvisited neighbor
        found_next = False
        for neighbor in adj[current]:
            if neighbor not in visited:
                ordered.append(neighbor)
                visited.add(neighbor)
                current = neighbor
                found_next = True
                break
        
        if not found_next:
            # Check if we've completed the loop
            if start in adj[current] and len(ordered) > 2:
                break  # Complete loop
            else:
                # Dead end - shouldn't happen with valid boundary
                break
    
    return np.array(ordered)


def _get_boundary_cell_types(boundary_pts, all_points, cells, cell_type_data):
    """
    Get cell type for each boundary point by finding which cell it belongs to.
    """
    boundary_cell_types = np.zeros(len(boundary_pts))
    
    for i, pt in enumerate(boundary_pts):
        # Find which cell contains this point as a vertex
        min_dist = float('inf')
        best_cell = 0
        
        for cell_idx, cell_verts in enumerate(cells):
            # Check if point matches any cell vertex
            cell_pts = all_points[cell_verts]
            dists = np.linalg.norm(cell_pts - pt, axis=1)
            if dists.min() < min_dist:
                min_dist = dists.min()
                best_cell = cell_idx
        
        boundary_cell_types[i] = cell_type_data[best_cell] if cell_type_data is not None else 0
    
    return boundary_cell_types


def load_final_vertex_boundary(data_dir):
    """
    Load the final (latest timestep) boundary from a vertex mesh simulation.
    
    Args:
        data_dir: Path to results_from_time_* directory
    
    Returns:
        boundary: Ordered boundary points (N x 2)
        cell_types: Cell type at each point (or None)
    """
    # Find all results VTU files (not ecm_grid files)
    vtu_files = [f for f in glob.glob(os.path.join(data_dir, 'results_*.vtu'))
                 if 'ecm_grid' not in f]
    
    if not vtu_files:
        raise FileNotFoundError(f"No results_*.vtu files in {data_dir}")
    
    # Get the one with the largest timestep
    def get_timestep(path):
        basename = os.path.basename(path)
        # results_6000.vtu -> 6000
        try:
            return int(basename.replace('results_', '').replace('.vtu', ''))
        except ValueError:
            return -1
    
    vtu_files.sort(key=get_timestep)
    latest_vtu = vtu_files[-1]
    
    print(f"Loading vertex boundary from: {os.path.basename(latest_vtu)}")
    return load_boundary_from_vtu(latest_vtu)


# =====================================================================
# Boundary extraction from cell positions
# =====================================================================

def extract_boundary_convex_hull(positions):
    """
    Extract ordered boundary points using convex hull.
    For point clouds from cell simulations, this gives the outer envelope.
    """
    if len(positions) < 3:
        return positions
    
    hull = ConvexHull(positions)
    # Get hull vertices in order
    boundary = positions[hull.vertices]
    return boundary


def extract_boundary_alpha_shape(positions, alpha=None):
    """
    Extract boundary using alpha shape (concave hull).
    Better than convex hull for capturing indentations.
    """
    try:
        from scipy.spatial import Delaunay
    except ImportError:
        return extract_boundary_convex_hull(positions)
    
    if len(positions) < 4:
        return positions
    
    # Compute Delaunay triangulation
    tri = Delaunay(positions)
    
    # Auto-compute alpha if not provided
    if alpha is None:
        # Use mean edge length as alpha
        edges = []
        for simplex in tri.simplices:
            for i in range(3):
                for j in range(i + 1, 3):
                    edge_len = np.linalg.norm(positions[simplex[i]] - positions[simplex[j]])
                    edges.append(edge_len)
        alpha = np.mean(edges) * 2.0
    
    # Find boundary edges (edges that appear in only one triangle)
    edge_count = {}
    for simplex in tri.simplices:
        # Check circumradius
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
            # Add edges from this triangle
            for i in range(3):
                edge = tuple(sorted([simplex[i], simplex[(i+1) % 3]]))
                edge_count[edge] = edge_count.get(edge, 0) + 1
    
    # Boundary edges appear exactly once
    boundary_edges = [e for e, count in edge_count.items() if count == 1]
    
    if not boundary_edges:
        return extract_boundary_convex_hull(positions)
    
    # Order boundary points
    boundary_points = order_boundary_edges(positions, boundary_edges)
    return boundary_points


def order_boundary_edges(positions, edges):
    """Order boundary edges into a continuous path."""
    if not edges:
        return np.array([])
    
    # Build adjacency
    adj = {}
    for e in edges:
        adj.setdefault(e[0], []).append(e[1])
        adj.setdefault(e[1], []).append(e[0])
    
    # Walk the boundary
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


# =====================================================================
# Chain code and Elliptic Fourier Descriptors
# =====================================================================

def compute_chain_code(boundary):
    """
    Compute Freeman 8-connected chain code from boundary points.
    The chain code describes direction changes along the boundary.
    """
    if len(boundary) < 2:
        return np.array([])
    
    # Direction to code mapping (8-connected)
    # Code: 0=right, 1=up-right, 2=up, 3=up-left, 4=left, 5=down-left, 6=down, 7=down-right
    
    # Compute relative displacements (circular - last point back to first)
    diff = np.vstack([np.diff(boundary, axis=0), boundary[0] - boundary[-1]])
    
    # Convert to discrete chain codes
    angles = np.arctan2(diff[:, 1], diff[:, 0])  # angle in radians
    # Map angle to chain code (0-7)
    codes = np.round((angles / (np.pi / 4)) % 8).astype(int)
    
    return codes


def calc_harmonic_coefficients(chain_code, n):
    """
    Calculate the nth harmonic coefficients (a_n, b_n, c_n, d_n).
    Based on elliptic Fourier descriptor formulation.
    """
    K = len(chain_code)
    if K == 0:
        return 0, 0, 0, 0
    
    # Compute traversal distances
    dt = np.ones(K)  # Assume unit steps for simplicity
    t = np.cumsum(dt)
    T = t[-1]
    
    # Delta x and delta y from chain code
    dx = np.cos(chain_code * np.pi / 4)
    dy = np.sin(chain_code * np.pi / 4)
    
    # Compute coefficients
    a_n = 0
    b_n = 0
    c_n = 0
    d_n = 0
    
    t_prev = np.concatenate([[0], t[:-1]])
    
    for i in range(K):
        if dt[i] < 1e-10:
            continue
        cos_term_curr = np.cos(2 * n * np.pi * t[i] / T)
        cos_term_prev = np.cos(2 * n * np.pi * t_prev[i] / T)
        sin_term_curr = np.sin(2 * n * np.pi * t[i] / T)
        sin_term_prev = np.sin(2 * n * np.pi * t_prev[i] / T)
        
        a_n += (dx[i] / dt[i]) * (cos_term_curr - cos_term_prev)
        b_n += (dx[i] / dt[i]) * (sin_term_curr - sin_term_prev)
        c_n += (dy[i] / dt[i]) * (cos_term_curr - cos_term_prev)
        d_n += (dy[i] / dt[i]) * (sin_term_curr - sin_term_prev)
    
    factor = T / (2 * n * n * np.pi * np.pi)
    return a_n * factor, b_n * factor, c_n * factor, d_n * factor


def calc_dc_components(chain_code):
    """Calculate DC components (A0, C0) for Fourier reconstruction."""
    K = len(chain_code)
    if K == 0:
        return 0, 0
    
    dt = np.ones(K)
    t = np.cumsum(dt)
    T = t[-1]
    
    dx = np.cos(chain_code * np.pi / 4)
    dy = np.sin(chain_code * np.pi / 4)
    
    # Cumulative sums
    xi = np.cumsum(dx)
    eta = np.cumsum(dy)
    
    A0 = 0
    C0 = 0
    
    t_prev = np.concatenate([[0], t[:-1]])
    xi_prev = np.concatenate([[0], xi[:-1]])
    eta_prev = np.concatenate([[0], eta[:-1]])
    
    for i in range(K):
        if dt[i] < 1e-10:
            continue
        A0 += (dx[i] / (2 * dt[i])) * (t[i]**2 - t_prev[i]**2) + \
              (xi_prev[i] - (dx[i] / dt[i]) * t_prev[i]) * dt[i]
        C0 += (dy[i] / (2 * dt[i])) * (t[i]**2 - t_prev[i]**2) + \
              (eta_prev[i] - (dy[i] / dt[i]) * t_prev[i]) * dt[i]
    
    return A0 / T, C0 / T


def fourier_reconstruct(boundary, n_harmonics=25, n_points=1000, normalize=True):
    """
    Reconstruct boundary using elliptic Fourier descriptors.
    
    Args:
        boundary: Original boundary points (N x 2)
        n_harmonics: Number of Fourier harmonics to use
        n_points: Number of points in reconstructed boundary
        normalize: Whether to normalize the descriptors
    
    Returns:
        Reconstructed boundary points (n_points x 2)
    """
    if len(boundary) < 4:
        return boundary
    
    # Compute chain code
    chain_code = compute_chain_code(boundary)
    if len(chain_code) == 0:
        return boundary
    
    # Compute harmonic coefficients
    a = np.zeros(n_harmonics)
    b = np.zeros(n_harmonics)
    c = np.zeros(n_harmonics)
    d = np.zeros(n_harmonics)
    
    for i in range(n_harmonics):
        a[i], b[i], c[i], d[i] = calc_harmonic_coefficients(chain_code, i + 1)
    
    A0, C0 = calc_dc_components(chain_code)
    
    # Normalization (optional)
    if normalize:
        A0, C0 = 0, 0
        
        # Compute normalization parameters from first harmonic
        if abs(a[0]) > 1e-10 or abs(c[0]) > 1e-10:
            theta1 = 0.5 * np.arctan2(2 * (a[0] * b[0] + c[0] * d[0]),
                                       a[0]**2 + c[0]**2 - b[0]**2 - d[0]**2)
            
            cos_th = np.cos(theta1)
            sin_th = np.sin(theta1)
            
            a_star = cos_th * a[0] + sin_th * b[0]
            c_star = cos_th * c[0] + sin_th * d[0]
            
            psi1 = np.arctan2(c_star, a_star)
            E = np.sqrt(a_star**2 + c_star**2)
            
            if E > 1e-10:
                cos_psi = np.cos(psi1)
                sin_psi = np.sin(psi1)
                
                for i in range(n_harmonics):
                    rot_angle = theta1 * (i + 1)
                    cos_rot = np.cos(rot_angle)
                    sin_rot = np.sin(rot_angle)
                    
                    # Apply rotation and normalization
                    mat1 = np.array([[cos_psi, sin_psi], [-sin_psi, cos_psi]])
                    mat2 = np.array([[a[i], b[i]], [c[i], d[i]]])
                    mat3 = np.array([[cos_rot, -sin_rot], [sin_rot, cos_rot]])
                    
                    result = mat1 @ mat2 @ mat3 / E
                    a[i], b[i] = result[0, 0], result[0, 1]
                    c[i], d[i] = result[1, 0], result[1, 1]
    
    # Reconstruct boundary
    output = np.zeros((n_points, 2))
    t_vals = np.arange(1, n_points + 1)
    
    for i in range(n_harmonics):
        n = i + 1
        cos_term = np.cos(2 * n * np.pi * t_vals / n_points)
        sin_term = np.sin(2 * n * np.pi * t_vals / n_points)
        output[:, 0] += a[i] * cos_term + b[i] * sin_term
        output[:, 1] += c[i] * cos_term + d[i] * sin_term
    
    output[:, 0] += A0
    output[:, 1] += C0
    
    # Scale to match original size
    original_scale = np.max(np.ptp(boundary, axis=0))
    reconstructed_scale = np.max(np.ptp(output, axis=0))
    if reconstructed_scale > 1e-10:
        output *= original_scale / reconstructed_scale
    
    # Center to match original
    output += boundary.mean(axis=0) - output.mean(axis=0)
    
    return output


# =====================================================================
# Curvature and normal computation
# =====================================================================

def compute_curvature(boundary):
    """
    Compute signed curvature at each point along a 2D boundary.
    Uses polynomial fitting similar to LineCurvature2D.m.
    
    Positive curvature = convex (bulging outward)
    Negative curvature = concave (crypt-like indentation)
    """
    n = len(boundary)
    if n < 3:
        return np.zeros(n)
    
    # Get neighbors (circular)
    prev_idx = np.roll(np.arange(n), 1)
    next_idx = np.roll(np.arange(n), -1)
    
    # Neighbors
    p_prev = boundary[prev_idx]
    p_curr = boundary
    p_next = boundary[next_idx]
    
    # Vectors
    v1 = p_curr - p_prev
    v2 = p_next - p_curr
    
    # Distances
    d1 = np.linalg.norm(v1, axis=1)
    d2 = np.linalg.norm(v2, axis=1)
    
    # Handle zero distances
    d1 = np.maximum(d1, 1e-10)
    d2 = np.maximum(d2, 1e-10)
    
    # Fit parabola and compute curvature
    # k = 2 * (v1 x v2) / (|v1| * |v2| * |v1 + v2|)
    cross = v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0]  # 2D cross product (scalar)
    denom = d1 * d2 * np.linalg.norm(v1 + v2, axis=1)
    denom = np.maximum(denom, 1e-10)
    
    curvature = 2 * cross / denom
    
    return curvature


def compute_normals(boundary):
    """
    Compute outward-pointing unit normals at each boundary point.
    """
    n = len(boundary)
    if n < 2:
        return np.zeros((n, 2))
    
    # Tangent vectors (forward difference, circular)
    next_idx = np.roll(np.arange(n), -1)
    prev_idx = np.roll(np.arange(n), 1)
    
    tangent = boundary[next_idx] - boundary[prev_idx]
    tangent_norm = np.linalg.norm(tangent, axis=1, keepdims=True)
    tangent_norm = np.maximum(tangent_norm, 1e-10)
    tangent = tangent / tangent_norm
    
    # Normal = perpendicular to tangent (rotate 90 degrees)
    # Assuming counter-clockwise boundary, outward normal points right
    normals = np.column_stack([-tangent[:, 1], tangent[:, 0]])
    
    # Check if normals point outward by testing against centroid
    centroid = boundary.mean(axis=0)
    to_centroid = centroid - boundary
    
    # If dot product with normal is positive, normal points inward - flip it
    dots = np.sum(normals * to_centroid, axis=1)
    flip_mask = dots > 0
    normals[flip_mask] *= -1
    
    return normals


# =====================================================================
# Crypt detection
# =====================================================================

def find_crypt_sections(boundary, curvature, normals, params=None, return_mask=False):
    """
    Identify crypt sections where curvature-scaled normals point inward.
    
    A crypt is a concave region where the normal (scaled by curvature)
    points into the interior of the organoid.
    
    Args:
        boundary: Smoothed boundary points (N x 2)
        curvature: Curvature at each point (N,)
        normals: Normal vectors at each point (N x 2)
        params: Parameters dict
        return_mask: If True, also return the inside mask
    
    Returns:
        sections: List of index lists for each section
        inside: (optional) Boolean mask of points with inward-pointing normals
    """
    if params is None:
        params = DEFAULT_PARAMS
    
    n = len(boundary)
    if n < 10:
        return ([], np.array([])) if return_mask else []
    
    # Compute endpoints of curvature-scaled normal vectors
    # When curvature is negative (concave), the scaled normal points inward
    scaled_normal_endpoints = boundary + curvature[:, np.newaxis] * normals
    
    # Check if endpoints are inside the polygon
    try:
        from matplotlib.path import Path
        polygon_path = Path(boundary)
        inside = polygon_path.contains_points(scaled_normal_endpoints)
    except ImportError:
        # Fallback: use cross product winding
        inside = points_in_polygon(scaled_normal_endpoints, boundary)
    
    # Find contiguous sections where normals point inward
    sections = []
    current_section = []
    
    for i in range(n):
        if inside[i]:
            current_section.append(i)
        else:
            if len(current_section) >= 3:  # Minimum size for a crypt
                sections.append(current_section)
            current_section = []
    
    # Handle wraparound
    if len(current_section) >= 3:
        # Check if it connects to the first section
        if sections and sections[0][0] == 0:
            sections[0] = current_section + sections[0]
        else:
            sections.append(current_section)
    
    if return_mask:
        return sections, inside
    return sections


def points_in_polygon(points, polygon):
    """Check if points are inside polygon using winding number."""
    n_poly = len(polygon)
    inside = np.zeros(len(points), dtype=bool)
    
    for i, pt in enumerate(points):
        winding = 0
        for j in range(n_poly):
            p1 = polygon[j]
            p2 = polygon[(j + 1) % n_poly]
            
            if p1[1] <= pt[1]:
                if p2[1] > pt[1]:
                    if ((p2[0] - p1[0]) * (pt[1] - p1[1]) - 
                        (pt[0] - p1[0]) * (p2[1] - p1[1])) > 0:
                        winding += 1
            else:
                if p2[1] <= pt[1]:
                    if ((p2[0] - p1[0]) * (pt[1] - p1[1]) - 
                        (pt[0] - p1[0]) * (p2[1] - p1[1])) < 0:
                        winding -= 1
        
        inside[i] = winding != 0
    
    return inside


def compute_arc_length(points):
    """Compute total arc length of a path."""
    if len(points) < 2:
        return 0
    diff = np.diff(points, axis=0)
    return np.sum(np.linalg.norm(diff, axis=1))


def compute_section_area(boundary, section_indices):
    """Compute area of a section polygon."""
    if len(section_indices) < 3:
        return 0
    
    section_points = boundary[section_indices]
    
    # Shoelace formula
    n = len(section_points)
    area = 0
    for i in range(n):
        j = (i + 1) % n
        area += section_points[i, 0] * section_points[j, 1]
        area -= section_points[j, 0] * section_points[i, 1]
    
    return abs(area) / 2


def filter_crypts(boundary, sections, params=None):
    """
    Filter detected crypt sections by area and arc length criteria.
    """
    if params is None:
        params = DEFAULT_PARAMS
    
    # Total organoid metrics
    total_area = compute_section_area(boundary, list(range(len(boundary))))
    total_arc_length = compute_arc_length(np.vstack([boundary, boundary[0:1]]))
    
    if total_area < 1e-10 or total_arc_length < 1e-10:
        return []
    
    valid_crypts = []
    
    for section in sections:
        # Compute normalized area
        section_area = compute_section_area(boundary, section)
        norm_area = section_area / total_area
        
        # Compute normalized arc length
        section_points = boundary[section]
        section_arc_length = compute_arc_length(section_points)
        norm_arc_length = section_arc_length / total_arc_length
        
        # Apply filters
        if (params['min_area'] <= norm_area < params['max_area'] and
            norm_arc_length >= params['min_arc_length']):
            valid_crypts.append({
                'indices': section,
                'norm_area': norm_area,
                'norm_arc_length': norm_arc_length
            })
    
    return valid_crypts


def compute_circularity(boundary):
    """
    Compute circularity of a 2D boundary.
    circularity = 4 * pi * Area / Perimeter^2
    Returns 1.0 for a perfect circle.
    """
    if len(boundary) < 3:
        return np.nan
    
    # Area using shoelace formula
    n = len(boundary)
    area = 0
    for i in range(n):
        j = (i + 1) % n
        area += boundary[i, 0] * boundary[j, 1]
        area -= boundary[j, 0] * boundary[i, 1]
    area = abs(area) / 2
    
    # Perimeter
    perimeter = compute_arc_length(np.vstack([boundary, boundary[0:1]]))
    
    if perimeter < 1e-10:
        return np.nan
    
    return 4 * np.pi * area / (perimeter ** 2)


# =====================================================================
# Main function
# =====================================================================

CryptCountResult = namedtuple('CryptCountResult', [
    'num_crypts', 'circularity', 'crypts', 'boundary', 'smooth_boundary',
    'curvature', 'normals', 'all_sections', 'inside_mask'
])


def count_crypts_simple_method(positions, params=None, use_alpha_shape=True, 
                                boundary_is_ordered=False):
    """
    Count crypts using the SimpleCryptCount method.
    
    Args:
        positions: Cell center positions (N x 2 array) OR pre-ordered boundary points
        params: Dict of parameters (see DEFAULT_PARAMS)
        use_alpha_shape: Use alpha shape (True) or convex hull (False) for boundary extraction
        boundary_is_ordered: If True, treat positions as already ordered boundary points
                            (e.g., from VTP outline file)
    
    Returns:
        CryptCountResult namedtuple with:
          - num_crypts: Number of detected crypts
          - circularity: Shape circularity (1.0 = circle)  
          - crypts: List of crypt info dicts (filtered)
          - boundary: Extracted/input boundary points
          - smooth_boundary: Fourier-smoothed boundary
          - curvature: Curvature at each smooth boundary point
          - normals: Normal vectors at each smooth boundary point
          - all_sections: All detected sections (before filtering)
          - inside_mask: Boolean mask indicating which points have inward normals
    """
    if params is None:
        params = DEFAULT_PARAMS.copy()
    
    if len(positions) < 10:
        return CryptCountResult(0, np.nan, [], np.array([]), np.array([]),
                                np.array([]), np.array([]), [], np.array([]))
    
    # Step 1: Get boundary
    if boundary_is_ordered:
        # Use positions directly as ordered boundary
        boundary = np.asarray(positions)
    else:
        # Extract boundary from cell positions
        if use_alpha_shape:
            boundary = extract_boundary_alpha_shape(positions)
        else:
            boundary = extract_boundary_convex_hull(positions)
    
    if len(boundary) < 10:
        return CryptCountResult(0, np.nan, [], boundary, np.array([]),
                                np.array([]), np.array([]), [], np.array([]))
    
    # Step 2: Smooth boundary using Fourier reconstruction
    smooth_boundary = fourier_reconstruct(
        boundary,
        n_harmonics=params.get('fourier_harmonics', 25),
        n_points=params.get('n_points', 1000),
        normalize=True
    )
    
    # Step 3: Compute circularity
    circularity = compute_circularity(smooth_boundary)
    
    # Step 4: Compute curvature and normals
    curvature = compute_curvature(smooth_boundary)
    normals = compute_normals(smooth_boundary)
    
    # Step 5: Find crypt sections (regions where curvature-scaled normals point inward)
    sections, inside_mask = find_crypt_sections(smooth_boundary, curvature, normals, 
                                                 params, return_mask=True)
    
    # If nearly circular, no crypts (but still return debug info)
    if circularity > params.get('circularity_threshold', 0.98):
        return CryptCountResult(0, circularity, [], boundary, smooth_boundary,
                                curvature, normals, sections, inside_mask)
    
    # Step 6: Filter by area and arc length
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


# =====================================================================
# Visualization (optional)
# =====================================================================

def plot_crypt_analysis(result, positions=None, output_path=None):
    """
    Visualize the crypt detection results.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available for plotting")
        return
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Plot 1: Original positions and boundary
    ax1 = axes[0]
    if positions is not None:
        ax1.scatter(positions[:, 0], positions[:, 1], s=5, alpha=0.5, label='Cells')
    if len(result.boundary) > 0:
        boundary_closed = np.vstack([result.boundary, result.boundary[0:1]])
        ax1.plot(boundary_closed[:, 0], boundary_closed[:, 1], 'b-', lw=2, label='Boundary')
    ax1.set_title('Cell Positions & Boundary')
    ax1.set_aspect('equal')
    ax1.legend()
    
    # Plot 2: Smooth boundary with curvature
    ax2 = axes[1]
    if len(result.smooth_boundary) > 0:
        smooth_closed = np.vstack([result.smooth_boundary, result.smooth_boundary[0:1]])
        ax2.plot(smooth_closed[:, 0], smooth_closed[:, 1], 'k-', lw=1)
        
        # Color by curvature
        sc = ax2.scatter(result.smooth_boundary[:, 0], result.smooth_boundary[:, 1],
                        c=result.curvature, cmap='RdBu', s=2)
        plt.colorbar(sc, ax=ax2, label='Curvature')
    ax2.set_title(f'Smooth Boundary (Circularity: {result.circularity:.3f})')
    ax2.set_aspect('equal')
    
    # Plot 3: Detected crypts
    ax3 = axes[2]
    if len(result.smooth_boundary) > 0:
        smooth_closed = np.vstack([result.smooth_boundary, result.smooth_boundary[0:1]])
        ax3.plot(smooth_closed[:, 0], smooth_closed[:, 1], 'k-', lw=1)
        
        # Highlight crypts
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
    """
    Create comprehensive debug visualization showing each step of Sandra's
    SimpleCryptCount method.
    
    Creates a 2x3 grid of subplots:
      1. Original boundary (from VTP or extracted)
      2. Fourier-smoothed boundary
      3. Curvature along boundary with color scale
      4. Curvature-scaled normals showing inside/outside detection
      5. All detected sections (before filtering)
      6. Final filtered crypts with area/arc-length annotations
    
    Args:
        result: CryptCountResult from count_crypts_simple_method
        positions: Original cell positions (optional)
        cell_types: Cell type values for boundary points (optional)
        output_path: Path to save figure (None shows interactive)
        params: Parameters used for detection
        title_prefix: Prefix for figure title (e.g., "stiffness_1.0/run_1")
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Polygon as MplPolygon
        from matplotlib.collections import LineCollection
        import matplotlib.cm as cm
    except ImportError:
        print("matplotlib not available for plotting")
        return
    
    if params is None:
        params = DEFAULT_PARAMS
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # =====================================================================
    # Plot 1: Original boundary
    # =====================================================================
    ax1 = axes[0, 0]
    if positions is not None and not np.array_equal(positions, result.boundary):
        ax1.scatter(positions[:, 0], positions[:, 1], s=10, alpha=0.3, 
                   c='gray', label='Cell positions')
    
    if len(result.boundary) > 0:
        boundary_closed = np.vstack([result.boundary, result.boundary[0:1]])
        
        if cell_types is not None:
            # Color boundary by cell type
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
        
        # Mark start point
        ax1.scatter(result.boundary[0, 0], result.boundary[0, 1], 
                   s=100, c='green', marker='*', zorder=5, label='Start point')
    
    ax1.set_title(f'1. Input Boundary\n({len(result.boundary)} points)')
    ax1.set_aspect('equal')
    ax1.legend(loc='upper right', fontsize=8)
    ax1.grid(True, alpha=0.3)
    
    # =====================================================================
    # Plot 2: Fourier-smoothed boundary
    # =====================================================================
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
    
    # =====================================================================
    # Plot 3: Curvature along boundary
    # =====================================================================
    ax3 = axes[0, 2]
    if len(result.smooth_boundary) > 0 and len(result.curvature) > 0:
        smooth_closed = np.vstack([result.smooth_boundary, result.smooth_boundary[0:1]])
        ax3.plot(smooth_closed[:, 0], smooth_closed[:, 1], 'k-', lw=0.5, alpha=0.3)
        
        # Color points by curvature
        vmax = max(abs(result.curvature.min()), abs(result.curvature.max()))
        sc = ax3.scatter(result.smooth_boundary[:, 0], result.smooth_boundary[:, 1],
                        c=result.curvature, cmap='RdBu_r', s=5, vmin=-vmax, vmax=vmax)
        cbar = plt.colorbar(sc, ax=ax3)
        cbar.set_label('Curvature\n(+ convex, - concave)')
        
        # Add annotations for curvature stats
        ax3.text(0.02, 0.98, f'Min: {result.curvature.min():.3f}\nMax: {result.curvature.max():.3f}',
                transform=ax3.transAxes, fontsize=8, va='top', ha='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax3.set_title('3. Curvature (LineCurvature2D)')
    ax3.set_aspect('equal')
    ax3.grid(True, alpha=0.3)
    
    # =====================================================================
    # Plot 4: Curvature-scaled normals (inside/outside detection)
    # =====================================================================
    ax4 = axes[1, 0]
    if len(result.smooth_boundary) > 0 and len(result.normals) > 0:
        smooth_closed = np.vstack([result.smooth_boundary, result.smooth_boundary[0:1]])
        ax4.plot(smooth_closed[:, 0], smooth_closed[:, 1], 'k-', lw=1)
        
        # Draw curvature-scaled normals (subsample for visibility)
        step = max(1, len(result.smooth_boundary) // 50)
        for i in range(0, len(result.smooth_boundary), step):
            pt = result.smooth_boundary[i]
            normal_end = pt + result.curvature[i] * result.normals[i]
            
            # Check if inside
            is_inside = result.inside_mask[i] if len(result.inside_mask) > 0 else False
            color = 'green' if is_inside else 'red'
            
            ax4.plot([pt[0], normal_end[0]], [pt[1], normal_end[1]], 
                    color=color, lw=1, alpha=0.7)
            ax4.scatter(normal_end[0], normal_end[1], c=color, s=10, alpha=0.7)
        
        # Legend
        ax4.scatter([], [], c='green', s=30, label='Normal tip INSIDE (crypt candidate)')
        ax4.scatter([], [], c='red', s=30, label='Normal tip OUTSIDE (organoid body)')
        
        # Count inside/outside
        n_inside = np.sum(result.inside_mask) if len(result.inside_mask) > 0 else 0
        n_outside = len(result.inside_mask) - n_inside if len(result.inside_mask) > 0 else 0
        ax4.text(0.02, 0.98, f'Inside: {n_inside}\nOutside: {n_outside}',
                transform=ax4.transAxes, fontsize=8, va='top', ha='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax4.set_title('4. Curvature-Scaled Normals\n(inpolygon test)')
    ax4.set_aspect('equal')
    ax4.legend(loc='upper right', fontsize=7)
    ax4.grid(True, alpha=0.3)
    
    # =====================================================================
    # Plot 5: All detected sections (before filtering)
    # =====================================================================
    ax5 = axes[1, 1]
    if len(result.smooth_boundary) > 0:
        smooth_closed = np.vstack([result.smooth_boundary, result.smooth_boundary[0:1]])
        ax5.plot(smooth_closed[:, 0], smooth_closed[:, 1], 'k-', lw=1, alpha=0.5)
        
        # Compute metrics for all sections
        total_area = compute_section_area(result.smooth_boundary, 
                                          list(range(len(result.smooth_boundary))))
        total_arc = compute_arc_length(smooth_closed)
        
        # Color each section
        colors = plt.cm.tab10(np.linspace(0, 1, max(len(result.all_sections), 1)))
        section_info = []
        
        for i, section in enumerate(result.all_sections):
            section_pts = result.smooth_boundary[section]
            ax5.plot(section_pts[:, 0], section_pts[:, 1], 
                    color=colors[i], lw=4, alpha=0.8)
            
            # Compute section metrics
            sec_area = compute_section_area(result.smooth_boundary, section)
            sec_arc = compute_arc_length(section_pts)
            norm_area = sec_area / total_area if total_area > 0 else 0
            norm_arc = sec_arc / total_arc if total_arc > 0 else 0
            
            # Mark midpoint with section number
            mid_idx = len(section) // 2
            mid_pt = section_pts[mid_idx]
            ax5.annotate(f'{i+1}', mid_pt, fontsize=8, ha='center', va='center',
                        bbox=dict(boxstyle='circle', facecolor=colors[i], alpha=0.8))
            
            section_info.append(f'S{i+1}: A={norm_area:.3f}, L={norm_arc:.3f}')
        
        # Add section info text box
        info_text = '\n'.join(section_info[:10])  # Limit to 10 sections
        if len(section_info) > 10:
            info_text += f'\n... +{len(section_info)-10} more'
        ax5.text(0.02, 0.98, info_text, transform=ax5.transAxes, fontsize=7, 
                va='top', ha='left', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                family='monospace')
    
    ax5.set_title(f'5. All Sections (before filtering)\n({len(result.all_sections)} sections)')
    ax5.set_aspect('equal')
    ax5.grid(True, alpha=0.3)
    
    # =====================================================================
    # Plot 6: Final filtered crypts
    # =====================================================================
    ax6 = axes[1, 2]
    if len(result.smooth_boundary) > 0:
        # Fill the organoid body in light gray
        smooth_closed = np.vstack([result.smooth_boundary, result.smooth_boundary[0:1]])
        ax6.fill(smooth_closed[:, 0], smooth_closed[:, 1], 
                color='lightgray', alpha=0.5, label='Organoid body')
        ax6.plot(smooth_closed[:, 0], smooth_closed[:, 1], 'k-', lw=1)
        
        # Highlight final crypts
        colors = plt.cm.Set1(np.linspace(0, 1, max(len(result.crypts), 1)))
        for i, crypt in enumerate(result.crypts):
            indices = crypt['indices']
            crypt_pts = result.smooth_boundary[indices]
            
            # Fill crypt region
            ax6.fill(crypt_pts[:, 0], crypt_pts[:, 1], 
                    color=colors[i], alpha=0.6)
            ax6.plot(crypt_pts[:, 0], crypt_pts[:, 1], 
                    color=colors[i], lw=3)
            
            # Annotate with metrics
            mid_idx = len(indices) // 2
            mid_pt = crypt_pts[mid_idx]
            ax6.annotate(f"Crypt {i+1}\nA={crypt['norm_area']:.3f}\nL={crypt['norm_arc_length']:.3f}",
                        mid_pt, fontsize=7, ha='center', va='center',
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
        
        # Add filter criteria
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
    
    # Overall title
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
    """
    Create a simple single-panel visualization showing the organoid outline
    with detected crypt regions shaded.
    
    This is a simplified version of the debug plot that just shows:
      - Organoid outline (black)
      - Organoid body (light gray fill)
      - Crypt regions (colored shading)
    
    Args:
        result: CryptCountResult from count_crypts_simple_method
        output_path: Path to save figure (None shows interactive)
        title: Custom title for the plot
        show_metrics: If True, show circularity and crypt count on plot
    
    Returns:
        matplotlib Figure object
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Polygon as MplPolygon
    except ImportError:
        print("matplotlib not available for plotting")
        return None
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    if len(result.smooth_boundary) > 0:
        # Close the boundary for plotting
        smooth_closed = np.vstack([result.smooth_boundary, result.smooth_boundary[0:1]])
        
        # Fill the organoid body in light gray
        ax.fill(smooth_closed[:, 0], smooth_closed[:, 1], 
                color='#E8E8E8', alpha=0.8, label='Organoid body')
        
        # Draw the outline
        ax.plot(smooth_closed[:, 0], smooth_closed[:, 1], 
                'k-', lw=2, label='Boundary')
        
        # Highlight detected crypts with distinct colors
        if len(result.crypts) > 0:
            crypt_colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', 
                           '#FFEAA7', '#DDA0DD', '#98D8C8', '#F7DC6F']
            
            for i, crypt in enumerate(result.crypts):
                indices = crypt['indices']
                crypt_pts = result.smooth_boundary[indices]
                
                # Create a closed polygon for the crypt region
                # Connect to centroid for proper fill
                centroid = result.smooth_boundary.mean(axis=0)
                
                # Fill crypt region (polygon from boundary to center)
                crypt_polygon = np.vstack([crypt_pts, centroid])
                color = crypt_colors[i % len(crypt_colors)]
                ax.fill(crypt_pts[:, 0], crypt_pts[:, 1], 
                       color=color, alpha=0.6, edgecolor=color, lw=2)
                
                # Add crypt label
                mid_idx = len(indices) // 2
                mid_pt = crypt_pts[mid_idx]
                
                # Offset label outward from centroid
                direction = mid_pt - centroid
                direction = direction / (np.linalg.norm(direction) + 1e-8)
                label_pt = mid_pt + direction * 0.3
                
                ax.annotate(f'{i+1}', mid_pt, fontsize=14, fontweight='bold',
                           ha='center', va='center', color='white',
                           bbox=dict(boxstyle='circle', facecolor=color, 
                                    edgecolor='white', alpha=0.9))
    
    # Add metrics text box
    if show_metrics:
        metrics_text = (f"Crypts: {result.num_crypts}\n"
                       f"Circularity: {result.circularity:.3f}")
        ax.text(0.02, 0.98, metrics_text, transform=ax.transAxes, 
                fontsize=12, va='top', ha='left',
                bbox=dict(boxstyle='round', facecolor='white', 
                         edgecolor='gray', alpha=0.9))
    
    # Set plot properties
    ax.set_aspect('equal')
    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    
    if title:
        ax.set_title(title, fontsize=14, fontweight='bold')
    else:
        ax.set_title(f'Organoid Outline with {result.num_crypts} Detected Crypts', 
                    fontsize=14, fontweight='bold')
    
    # Remove grid for cleaner look
    ax.grid(False)
    
    # Add subtle background
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


# =====================================================================
# CLI Entry Point
# =====================================================================

if __name__ == '__main__':
    import argparse
    import sys
    
    # Add parent directory to path for imports
    sys.path.insert(0, '.')
    
    parser = argparse.ArgumentParser(
        description='SimpleCryptCount method for Chaste simulation output',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyse using VTP outline file (recommended for node2d):
  python simple_crypt_count.py --data-dir /path/to/results_from_time_0 --use-outline --debug
  
  # Analyse vertex model data (vertex2d):
  python simple_crypt_count.py --data-dir /path/to/results_from_time_0 --use-vertex-mesh --debug
  
  # Analyse using cell positions (fallback):
  python simple_crypt_count.py --data-dir /path/to/results_from_time_0
  
  # Custom parameters:
  python simple_crypt_count.py --data-dir /path/to/results --use-outline --min-area 0.05 --max-area 0.3
""")
    
    parser.add_argument('--data-dir', type=str, required=True,
                        help='Directory with simulation output')
    parser.add_argument('--use-outline', action='store_true',
                        help='Load boundary from VTP outline files (node2d)')
    parser.add_argument('--use-vertex-mesh', action='store_true',
                        help='Load boundary from VTU vertex mesh files (vertex2d)')
    parser.add_argument('--output', type=str, default=None,
                        help='Output plot path (optional)')
    parser.add_argument('--debug', action='store_true',
                        help='Generate detailed debug visualization')
    parser.add_argument('--fourier-harmonics', type=int, default=25,
                        help='Number of Fourier harmonics (default: 25)')
    parser.add_argument('--min-area', type=float, default=0.0666,
                        help='Minimum normalized crypt area (default: 0.0666)')
    parser.add_argument('--max-area', type=float, default=0.2736,
                        help='Maximum normalized crypt area (default: 0.2736)')
    parser.add_argument('--min-arc-length', type=float, default=0.1466,
                        help='Minimum normalized arc length (default: 0.1466)')
    
    args = parser.parse_args()
    
    # Set up parameters
    params = DEFAULT_PARAMS.copy()
    params['fourier_harmonics'] = args.fourier_harmonics
    params['min_area'] = args.min_area
    params['max_area'] = args.max_area
    params['min_arc_length'] = args.min_arc_length
    
    positions = None
    cell_types = None
    boundary_is_ordered = False
    
    if args.use_outline:
        # Load from VTP outline file (node2d model)
        try:
            boundary, cell_types = load_final_outline(args.data_dir)
            positions = boundary
            boundary_is_ordered = True
            print(f"Loaded outline with {len(boundary)} points from VTP file")
        except FileNotFoundError as e:
            print(f"Error: {e}")
            print("No outline_*.vtp files found. Try --use-vertex-mesh for vertex2d.")
            sys.exit(1)
    elif args.use_vertex_mesh:
        # Load from VTU vertex mesh file (vertex2d model)
        try:
            boundary, cell_types = load_final_vertex_boundary(args.data_dir)
            positions = boundary
            boundary_is_ordered = True
            print(f"Loaded vertex mesh boundary with {len(boundary)} points from VTU file")
        except FileNotFoundError as e:
            print(f"Error: {e}")
            print("No results_*.vtu files found.")
            sys.exit(1)
        except Exception as e:
            print(f"Error loading vertex mesh: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)
    else:
        # Load cell positions
        try:
            from analyse_crypt_budding import load_final_positions
            positions = load_final_positions(args.data_dir, dim=2)
            print(f"Loaded {len(positions)} cell positions")
        except ImportError:
            print("Could not import load_final_positions from analyse_crypt_budding.py")
            sys.exit(1)
        except Exception as e:
            print(f"Error loading data: {e}")
            sys.exit(1)
    
    # Run analysis
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
    
    # Generate plots
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
