#!/usr/bin/env python3
"""Inkscape extension wrapper for svgexplode."""

import sys
import os
import tempfile

# Locate project dir: the extension files should be siblings of svgexplode.py,
# or svgexplode.py should be in the same directory. If not, fall back to
# the directory containing this extension script.
EXT_DIR = os.path.dirname(os.path.abspath(__file__))

# Add project dir (where svgexplode.py lives) to the path
# Check if svgexplode.py is alongside this extension file
if os.path.isfile(os.path.join(EXT_DIR, "svgexplode.py")):
    PROJECT_DIR = EXT_DIR
else:
    # Fallback: try a config file that points to the project
    PROJECT_DIR = EXT_DIR
if PROJECT_DIR not in sys.path:
    sys.path.insert(0, PROJECT_DIR)

# Add venv site-packages so we can import svgpathtools/numpy/scipy.
# Look for .venv in the project dir first, then alongside svgexplode.py.
VENV_BASE = os.path.join(PROJECT_DIR, ".venv")
if os.path.isdir(VENV_BASE):
    # Windows: .venv/Lib/site-packages
    win_sp = os.path.join(VENV_BASE, "Lib", "site-packages")
    if os.path.isdir(win_sp) and win_sp not in sys.path:
        sys.path.insert(0, win_sp)
    # Linux/macOS: .venv/lib/pythonX.Y/site-packages
    lib_dir = os.path.join(VENV_BASE, "lib")
    if os.path.isdir(lib_dir):
        for d in os.listdir(lib_dir):
            sp = os.path.join(lib_dir, d, "site-packages")
            if os.path.isdir(sp) and sp not in sys.path:
                sys.path.insert(0, sp)

import inkex
import svgpathtools as spt
from svgexplode import (
    find_all_intersections,
    split_paths_at_events,
    build_planar_graph,
    prune_graph,
    find_faces,
    face_to_path,
    signed_area,
    pack_regions,
)


class SvgExplodeEffect(inkex.EffectExtension):

    def add_arguments(self, pars):
        pars.add_argument("--spacing", type=float, default=5.0)
        pars.add_argument("--grid_res", type=float, default=0.25)
        pars.add_argument("--rotations", type=str, default="cardinal")
        pars.add_argument("--pack", type=inkex.Boolean, default=True)

    def effect(self):
        spacing = self.options.spacing
        grid_res = self.options.grid_res
        do_pack = self.options.pack

        rot_map = {
            "none": (0,),
            "cardinal": (0, 90, 180, 270),
            "fine": tuple(range(0, 360, 45)),
            "finest": tuple(range(0, 360, 15)),
        }
        rotations = rot_map.get(self.options.rotations, (0, 90, 180, 270))

        # Collect paths from selection, or all paths if nothing selected
        svg_paths = []
        path_elements = []

        if self.svg.selection:
            for elem in self.svg.selection.filter(inkex.PathElement):
                path_elements.append(elem)
        else:
            for elem in self.svg.xpath("//svg:path"):
                path_elements.append(elem)

        if len(path_elements) < 2:
            inkex.errormsg("SVG Explode needs at least 2 paths to find intersections.")
            return

        # Convert inkex paths to svgpathtools paths
        for elem in path_elements:
            d = elem.get("d")
            if d:
                try:
                    svg_paths.append(spt.parse_path(d))
                except Exception:
                    pass

        if len(svg_paths) < 2:
            inkex.errormsg("Could not parse enough valid paths.")
            return

        # Run the core algorithm
        events = find_all_intersections(svg_paths)
        segments = split_paths_at_events(svg_paths, events)
        graph_nodes, graph_edges = build_planar_graph(segments)
        pruned_edges = prune_graph(graph_nodes, graph_edges)
        faces = find_faces(graph_nodes, pruned_edges, segments)

        result = []
        for face in faces:
            try:
                p = face_to_path(face, segments)
                area = signed_area(p)
                result.append((p, area))
            except Exception:
                pass

        if not result:
            inkex.errormsg("No closed regions found.")
            return

        # Remove outer face
        result.sort(key=lambda x: abs(x[1]))
        result.pop()

        if not result:
            inkex.errormsg("No inner closed regions found.")
            return

        # Pack or keep in place
        if do_pack:
            placed_paths, total_w, total_h = pack_regions(
                result, spacing, grid_res, rotations)
        else:
            placed_paths = [p for p, a in result]

        # Remove original path elements from the SVG
        parent = self.svg.get_current_layer()
        for elem in path_elements:
            elem.getparent().remove(elem)

        # Add the new region paths
        colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3',
                  '#ff7f00', '#a65628', '#f781bf', '#999999']

        for i, p in enumerate(placed_paths):
            color = colors[i % len(colors)]
            elem = inkex.PathElement()
            d_str = p.d()
            if p.isclosed():
                d_str += " Z"
            elem.set("d", d_str)
            elem.set("id", f"region_{i}")
            elem.style = inkex.Style({
                "fill": color,
                "fill-opacity": "0.4",
                "stroke": color,
                "stroke-width": "0.5",
            })
            parent.append(elem)


if __name__ == "__main__":
    SvgExplodeEffect().run()
