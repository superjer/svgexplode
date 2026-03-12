#!/usr/bin/env python3
"""
svgexplode - Break apart intersecting SVG paths into separate closed regions.

Given an SVG with intersecting paths, finds all intersection points,
splits paths at those points, builds a planar graph, and extracts every
enclosed face as a separate closed SVG path.
"""

import sys
import math
from collections import defaultdict
import numpy as np
import svgpathtools as spt

TOLERANCE = 0.5  # point clustering tolerance in SVG user units


def find_all_intersections(paths):
    """Find all intersection points between paths.

    Returns a list of (point, path_index, T_value) events.
    Handles proper crossings, shared endpoints, and T-junctions.
    """
    events = []

    # Proper crossings between path pairs
    for i in range(len(paths)):
        for j in range(i + 1, len(paths)):
            try:
                intersections = paths[i].intersect(paths[j], tol=1e-8)
            except Exception:
                intersections = []
            for (T1, seg1, t1), (T2, seg2, t2) in intersections:
                pt = paths[i].point(T1)
                events.append((pt, i, T1))
                events.append((pt, j, T2))

    # Collect all knot points (path endpoints AND internal segment joins)
    # Use t2T to get correct global T values at segment boundaries
    knot_points = []  # (point, path_index, T_value)
    for i, p in enumerate(paths):
        for k, seg in enumerate(p):
            T_start = p.t2T(k, 0)
            knot_points.append((seg.start, i, T_start))
        # Also add the final endpoint
        T_end = p.t2T(len(p) - 1, 1)
        knot_points.append((p[-1].end, i, T_end))

    # Check shared knot points between different paths
    for a in range(len(knot_points)):
        for b in range(a + 1, len(knot_points)):
            pt_a, idx_a, T_a = knot_points[a]
            pt_b, idx_b, T_b = knot_points[b]
            if idx_a != idx_b and abs(pt_a - pt_b) < TOLERANCE:
                events.append((pt_a, idx_a, T_a))
                events.append((pt_b, idx_b, T_b))

    # Check if any knot point lies on the interior of another path (T-junctions)
    for kp, kp_idx, kp_T in knot_points:
        for j, path in enumerate(paths):
            if j == kp_idx:
                continue
            # Sample the path to find closest point
            min_dist = float('inf')
            min_t = None
            for t in np.linspace(0, 1, 2000):
                dist = abs(kp - path.point(t))
                if dist < min_dist:
                    min_dist = dist
                    min_t = t
            if min_dist < TOLERANCE and 0.005 < min_t < 0.995:
                # Refine with finer sampling around min_t
                for t in np.linspace(max(0, min_t - 0.005), min(1, min_t + 0.005), 500):
                    dist = abs(kp - path.point(t))
                    if dist < min_dist:
                        min_dist = dist
                        min_t = t
                if min_dist < TOLERANCE:
                    events.append((kp, kp_idx, kp_T))
                    events.append((kp, j, min_t))

    # Deduplicate events
    seen = set()
    unique = []
    for pt, idx, T in events:
        key = (idx, round(T, 6))
        if key not in seen:
            seen.add(key)
            unique.append((pt, idx, T))

    return unique


def split_paths_at_events(paths, events):
    """Split each path at intersection T-values, returning sub-segments.

    Returns list of (sub_path, start_point, end_point).
    """
    # Group T-values by path index
    path_ts = defaultdict(set)
    for pt, idx, T in events:
        path_ts[idx].add(T)

    # Always include endpoints
    for i in range(len(paths)):
        path_ts[i].add(0.0)
        path_ts[i].add(1.0)

    segments = []
    for i, path in enumerate(paths):
        t_values = sorted(path_ts[i])
        for k in range(len(t_values) - 1):
            t0, t1 = t_values[k], t_values[k + 1]
            if t1 - t0 < 1e-8:
                continue
            try:
                sub = path.cropped(t0, t1)
                segments.append((sub, path.point(t0), path.point(t1)))
            except Exception as e:
                print(f"  Warning: crop failed on path {i} [{t0:.4f}, {t1:.4f}]: {e}",
                      file=sys.stderr)

    return segments


def cluster_point(pt, nodes, tol=TOLERANCE):
    """Find or create a node for a point. Returns node index."""
    for i, n in enumerate(nodes):
        if abs(pt - n) < tol:
            return i
    nodes.append(pt)
    return len(nodes) - 1


def build_planar_graph(segments):
    """Build a planar graph from sub-segments.

    Returns (nodes, edges) where:
      nodes: list of complex points
      edges: list of (node_a, node_b, segment_index)
    """
    nodes = []
    edges = []

    for seg_idx, (sub, start, end) in enumerate(segments):
        na = cluster_point(start, nodes)
        nb = cluster_point(end, nodes)
        if na != nb:
            edges.append((na, nb, seg_idx))

    return nodes, edges


def prune_graph(nodes, edges):
    """Iteratively remove degree-1 nodes (dangling edges can't form faces)."""
    edge_set = list(edges)
    changed = True
    while changed:
        changed = False
        degree = defaultdict(int)
        for a, b, si in edge_set:
            degree[a] += 1
            degree[b] += 1
        new_edges = []
        for a, b, si in edge_set:
            if degree[a] >= 2 and degree[b] >= 2:
                new_edges.append((a, b, si))
            else:
                changed = True
        edge_set = new_edges
    return edge_set


def outgoing_angle(nodes, from_node, segments, seg_idx, forward):
    """Compute the outgoing angle of an edge from a node.

    Samples slightly along the curve to handle curvature at the node.
    """
    sub = segments[seg_idx][0]
    if forward:
        pt0 = sub.point(0)
        pt1 = sub.point(0.02)
    else:
        pt0 = sub.point(1)
        pt1 = sub.point(0.98)
    dx = pt1.real - pt0.real
    dy = pt1.imag - pt0.imag
    return math.atan2(dy, dx)


def find_faces(nodes, edges, segments):
    """Find all minimal faces using planar face traversal.

    At each node, edges are sorted by outgoing angle. For each directed
    half-edge u->v, the next half-edge in the same face is found by:
    at node v, locate the twin edge v->u in the sorted list, then take
    the previous entry (next clockwise in SVG's y-down coords).
    """
    # Build adjacency: node -> sorted list of (to_node, seg_idx, forward, angle)
    adj = defaultdict(list)
    for na, nb, seg_idx in edges:
        ang_fwd = outgoing_angle(nodes, na, segments, seg_idx, True)
        adj[na].append((nb, seg_idx, True, ang_fwd))
        ang_bwd = outgoing_angle(nodes, nb, segments, seg_idx, False)
        adj[nb].append((na, seg_idx, False, ang_bwd))

    for node in adj:
        adj[node].sort(key=lambda x: x[3])

    # Build twin lookup: at each node, find the index of each half-edge
    # Key: (from_node, to_node, seg_idx, forward) -> index in adj[from_node]
    he_index = {}
    for node in adj:
        for i, (to, si, fwd, ang) in enumerate(adj[node]):
            he_index[(node, to, si, fwd)] = i

    def next_half_edge(u, v, si, fwd):
        """Given half-edge u->v, return the next half-edge in the face."""
        # Find the twin (v->u) in adj[v]
        twin_key = (v, u, si, not fwd)
        twin_idx = he_index.get(twin_key)
        if twin_idx is None:
            return None
        # Next clockwise = previous index in the CCW-sorted list
        nxt_idx = (twin_idx - 1) % len(adj[v])
        to, si2, fwd2, ang2 = adj[v][nxt_idx]
        return (v, to, si2, fwd2)

    # Traverse all faces
    visited = set()
    faces = []

    for na, nb, seg_idx in edges:
        for fwd in (True, False):
            start = (na if fwd else nb, nb if fwd else na, seg_idx, fwd)
            if start in visited:
                continue

            face = []
            current = start
            ok = True
            for _ in range(len(edges) * 2 + 2):  # safety bound
                if current in visited and current != start:
                    ok = False
                    break
                visited.add(current)
                face.append(current)
                nxt = next_half_edge(*current)
                if nxt is None:
                    ok = False
                    break
                if nxt == start:
                    break  # closed face
                current = nxt

            if ok and len(face) >= 2:
                faces.append(face)

    return faces


def face_to_path(face, segments):
    """Convert a face (list of half-edges) into a single closed svgpathtools.Path.

    Stitches segments end-to-start so there's only one M command,
    and appends a Z (close) at the end.
    """
    path_segs = []
    for u, v, si, fwd in face:
        sub = segments[si][0]
        if fwd:
            segs = list(sub)
        else:
            segs = [s.reversed() for s in reversed(sub)]

        for s in segs:
            if path_segs:
                # Force this segment's start to match the previous end
                # so svgpathtools doesn't insert a new M command
                prev_end = path_segs[-1].end
                if isinstance(s, spt.CubicBezier):
                    s = spt.CubicBezier(prev_end, s.control1, s.control2, s.end)
                elif isinstance(s, spt.QuadraticBezier):
                    s = spt.QuadraticBezier(prev_end, s.control, s.end)
                elif isinstance(s, spt.Line):
                    s = spt.Line(prev_end, s.end)
                elif isinstance(s, spt.Arc):
                    s = spt.Arc(prev_end, s.radius, s.rotation, s.large_arc,
                                s.sweep, s.end)
            path_segs.append(s)

    # Close: snap the last segment's end to the first segment's start
    if path_segs:
        first_start = path_segs[0].start
        last = path_segs[-1]
        if abs(last.end - first_start) > 1e-6:
            # Add a closing line if there's a real gap
            path_segs.append(spt.Line(last.end, first_start))
        else:
            # Snap the endpoint exactly to close cleanly
            if isinstance(last, spt.CubicBezier):
                path_segs[-1] = spt.CubicBezier(last.start, last.control1,
                                                  last.control2, first_start)
            elif isinstance(last, spt.QuadraticBezier):
                path_segs[-1] = spt.QuadraticBezier(last.start, last.control,
                                                      first_start)
            elif isinstance(last, spt.Line):
                path_segs[-1] = spt.Line(last.start, first_start)

    p = spt.Path(*path_segs)
    p.closed = True
    return p


def signed_area(path, n_samples=300):
    """Compute the signed area of a path using the shoelace formula on samples."""
    pts = [path.point(t / n_samples) for t in range(n_samples)]
    area = 0.0
    for i in range(n_samples):
        j = (i + 1) % n_samples
        area += pts[i].real * pts[j].imag
        area -= pts[j].real * pts[i].imag
    return area / 2.0


def rotate_path(path, angle_deg):
    """Rotate a path around its bounding box center by angle_deg degrees."""
    bbox = path.bbox()
    cx = (bbox[0] + bbox[1]) / 2
    cy = (bbox[2] + bbox[3]) / 2
    center = cx + cy * 1j
    angle_rad = math.radians(angle_deg)
    rot = math.cos(angle_rad) + 1j * math.sin(angle_rad)
    return path.translated(-center).scaled(rot).translated(center)


def path_to_polygon(path, n_samples=1000):
    """Sample a path into a polygon (list of (x, y) tuples)."""
    return [(path.point(t / n_samples).real, path.point(t / n_samples).imag)
            for t in range(n_samples)]


def rasterize_polygon(polygon, grid_res):
    """Rasterize a polygon to a binary mask relative to its bounding box.

    Returns (mask, x_offset, y_offset) where offsets are in grid coords.
    mask[row, col] = True means the polygon covers that cell.
    """
    xs = [p[0] for p in polygon]
    ys = [p[1] for p in polygon]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)

    cols = int(math.ceil((xmax - xmin) / grid_res)) + 1
    rows = int(math.ceil((ymax - ymin) / grid_res)) + 1
    mask = np.zeros((rows, cols), dtype=bool)

    # Scanline fill: for each row, find polygon edge crossings
    poly_pts = [(int((x - xmin) / grid_res), int((y - ymin) / grid_res))
                for x, y in polygon]
    n = len(poly_pts)

    for row in range(rows):
        crossings = []
        for i in range(n):
            x0, y0 = poly_pts[i]
            x1, y1 = poly_pts[(i + 1) % n]
            if y0 == y1:
                continue
            if min(y0, y1) <= row < max(y0, y1):
                # X at this row via linear interpolation
                t = (row - y0) / (y1 - y0)
                cx = x0 + t * (x1 - x0)
                crossings.append(cx)
        crossings.sort()
        # Fill between pairs of crossings
        for j in range(0, len(crossings) - 1, 2):
            c_start = max(0, int(math.floor(crossings[j])))
            c_end = min(cols, int(math.ceil(crossings[j + 1])))
            mask[row, c_start:c_end] = True

    return mask, xmin, ymin


def dilate_mask(mask, pixels):
    """Dilate a binary mask by the given number of pixels."""
    if pixels <= 0:
        return mask
    from scipy.ndimage import binary_dilation
    struct_size = 2 * pixels + 1
    struct = np.ones((struct_size, struct_size), dtype=bool)
    # Make it circular
    center = pixels
    for r in range(struct_size):
        for c in range(struct_size):
            if (r - center) ** 2 + (c - center) ** 2 > (pixels + 0.5) ** 2:
                struct[r, c] = False
    return binary_dilation(mask, structure=struct).astype(bool)


def prepare_piece(path, area, grid_res, spacing_px):
    """Prepare a piece for packing: rasterize and dilate."""
    poly = path_to_polygon(path)
    mask, ox, oy = rasterize_polygon(poly, grid_res)
    half_sp = max(1, spacing_px // 2)
    dmask = dilate_mask(mask, half_sp)
    bbox = path.bbox()
    return {
        'path': path,
        'mask': mask,
        'dmask': dmask,
        'bbox': bbox,
        'area': area,
    }


def find_best_placement(occupied, dmask, grid_h, grid_w, spacing_px, step):
    """Find the best (lowest, then leftmost) position for a piece."""
    mh, mw = dmask.shape
    best_pos = None
    best_y = grid_h
    best_x = grid_w

    for gy in range(spacing_px, grid_h - mh, step):
        if gy > best_y:
            break
        for gx in range(spacing_px, grid_w - mw, step):
            region = occupied[gy:gy + mh, gx:gx + mw]
            if not np.any(region & dmask):
                if gy < best_y or (gy == best_y and gx < best_x):
                    best_pos = (gx, gy)
                    best_y = gy
                    best_x = gx
                break
    return best_pos, best_y, best_x


def pack_regions(result, spacing, grid_res, rotations=(0, 90, 180, 270)):
    """Pack region paths tightly using raster-based collision detection.

    For each piece, tries multiple rotations and picks the one that
    packs tightest (lowest y, then leftmost x).
    Uses a bottom-left gravity heuristic with actual shape collision.
    """
    spacing_px = int(math.ceil(spacing / grid_res))

    # Pre-sort by area (largest first)
    result_sorted = sorted(result, key=lambda r: -abs(r[1]))

    # Master grid — 3:2 aspect ratio, sized from total piece area
    total_area = sum(abs(a) for _, a in result_sorted)
    # For a 3:2 rectangle: w = 3k, h = 2k, area = 6k²
    # We want ~4x the piece area for headroom
    k = math.sqrt(total_area * 4 / 6)
    grid_w_world = max(3 * k, 200)
    grid_h_world = grid_w_world * 2 / 3
    grid_w = int(grid_w_world / grid_res)
    grid_h = int(grid_h_world / grid_res)
    occupied = np.zeros((grid_h, grid_w), dtype=bool)

    placed_paths = []
    step = max(1, spacing_px // 3)  # finer scan step for tighter packing

    for pi, (path, area) in enumerate(result_sorted):
        best_overall = None  # (y, x, gx, gy, piece_data, rotated_path)

        for angle in rotations:
            if angle == 0:
                rpath = path
            else:
                rpath = rotate_path(path, angle)

            piece = prepare_piece(rpath, area, grid_res, spacing_px)
            dmask = piece['dmask']
            mh, mw = dmask.shape

            if mh >= grid_h - spacing_px or mw >= grid_w - spacing_px:
                continue  # piece doesn't fit at this rotation

            pos, by, bx = find_best_placement(
                occupied, dmask, grid_h, grid_w, spacing_px, step)

            if pos is not None:
                if (best_overall is None or
                        by < best_overall[0] or
                        (by == best_overall[0] and bx < best_overall[1])):
                    best_overall = (by, bx, pos[0], pos[1], piece, rpath, angle)

        if best_overall is None:
            print(f"  Warning: could not place region {pi}, using fallback",
                  file=sys.stderr)
            piece = prepare_piece(path, area, grid_res, spacing_px)
            gx, gy = spacing_px, spacing_px
            rpath = path
            angle = 0
        else:
            _, _, gx, gy, piece, rpath, angle = best_overall

        # Mark occupied
        dmask = piece['dmask']
        occupied[gy:gy + dmask.shape[0], gx:gx + dmask.shape[1]] |= dmask

        # Translate to placed position
        bbox = piece['bbox']
        target_x = gx * grid_res
        target_y = gy * grid_res
        dx = target_x - bbox[0]
        dy = target_y - bbox[2]
        translated = rpath.translated(dx + dy * 1j)
        placed_paths.append(translated)

        rot_str = f" rot={angle}°" if angle else ""
        print(f"  Placed region {pi} at ({target_x:.1f}, {target_y:.1f}){rot_str}")

    # Compute actual bounds from occupied grid
    occ_rows = np.any(occupied, axis=1)
    occ_cols = np.any(occupied, axis=0)
    if np.any(occ_rows) and np.any(occ_cols):
        max_row = np.max(np.where(occ_rows)) + spacing_px
        max_col = np.max(np.where(occ_cols)) + spacing_px
        total_w = max_col * grid_res
        total_h = max_row * grid_res
    else:
        total_w = 210
        total_h = 297

    return placed_paths, total_w, total_h


def main():
    if len(sys.argv) < 2:
        print("Usage: svgexplode.py input.svg [output.svg]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else input_file.replace('.svg', '_exploded.svg')

    print(f"Reading {input_file}...")
    paths, attributes = spt.svg2paths(input_file)
    print(f"Found {len(paths)} paths")

    # Step 1: Find intersections
    print("Finding intersections...")
    events = find_all_intersections(paths)
    intersection_pts = set()
    for pt, idx, T in events:
        intersection_pts.add(round(pt.real, 2) + round(pt.imag, 2) * 1j)
    print(f"  {len(intersection_pts)} unique intersection points, {len(events)} split events")

    # Step 2: Split paths
    print("Splitting paths at intersections...")
    segments = split_paths_at_events(paths, events)
    print(f"  {len(segments)} sub-segments")

    # Step 3: Build planar graph
    print("Building planar graph...")
    graph_nodes, graph_edges = build_planar_graph(segments)
    print(f"  {len(graph_nodes)} nodes, {len(graph_edges)} edges")

    # Step 4: Prune dangling edges
    pruned_edges = prune_graph(graph_nodes, graph_edges)
    removed = len(graph_edges) - len(pruned_edges)
    if removed:
        print(f"  Pruned {removed} dangling edges, {len(pruned_edges)} remain")

    # Step 5: Find faces
    print("Finding enclosed regions...")
    faces = find_faces(graph_nodes, pruned_edges, segments)
    print(f"  {len(faces)} candidate faces")

    # Step 6: Convert to paths and filter
    result = []
    for face in faces:
        try:
            p = face_to_path(face, segments)
            area = signed_area(p)
            result.append((p, area))
        except Exception as e:
            print(f"  Warning: face conversion failed: {e}", file=sys.stderr)

    if not result:
        print("No closed regions found!")
        sys.exit(0)

    # Remove the outer face (largest absolute area)
    result.sort(key=lambda x: abs(x[1]))
    outer = result.pop()
    print(f"  Outer face area: {abs(outer[1]):.1f}")
    print(f"  {len(result)} inner regions")

    if not result:
        print("No inner closed regions found!")
        sys.exit(0)

    # Step 7: Pack regions tightly with spacing
    SPACING = 5  # gap between shapes in SVG user units (mm)
    GRID_RES = 0.25  # raster resolution in SVG units per pixel
    placed_paths, total_w, total_h = pack_regions(result, SPACING, GRID_RES)

    out_paths = []
    out_attrs = []
    colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3',
              '#ff7f00', '#a65628', '#f781bf', '#999999']

    for i, p in enumerate(placed_paths):
        color = colors[i % len(colors)]
        out_paths.append(p)
        out_attrs.append({
            'style': f'fill:{color};fill-opacity:0.4;stroke:{color};stroke-width:0.5',
            'id': f'region_{i}',
        })

    spt.wsvg(out_paths, filename=output_file, attributes=out_attrs,
             svg_attributes={
                 'viewBox': f'0 0 {total_w:.1f} {total_h:.1f}',
                 'width': f'{total_w:.1f}mm',
                 'height': f'{total_h:.1f}mm',
             })

    # svgpathtools doesn't write Z for closed paths; add it
    import re
    with open(output_file, 'r') as f:
        svg_text = f.read()
    svg_text = re.sub(
        r'(<path d="[^"]*)("\s+id="region_)',
        r'\1 Z\2',
        svg_text
    )
    with open(output_file, 'w') as f:
        f.write(svg_text)

    print(f"Written to {output_file}")


if __name__ == '__main__':
    main()
