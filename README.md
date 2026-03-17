# svgexplode

Break apart intersecting SVG paths into separate closed regions.

Given an SVG with overlapping bezier curves and lines, svgexplode finds every intersection point, splits the paths, and extracts each enclosed face as its own closed SVG path. The output regions are packed tightly for CNC cutting from sheet material.

## How it works

1. **Find intersections** between all path pairs (proper crossings, shared knot points, and T-junctions)
2. **Split paths** at intersection points into sub-segments
3. **Build a planar graph** from the sub-segments
4. **Prune dangling edges** that can't form closed regions
5. **Traverse faces** using the planar face enumeration algorithm
6. **Pack regions** using raster-based collision detection with rotation search

## Installation

Requires Python 3.8+.

```bash
git clone <repo-url>
cd svgexplode
```

**Linux / macOS:**

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install svgpathtools numpy scipy
```

**Windows (PowerShell):**

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install svgpathtools numpy scipy
```

## Usage

### Standalone CLI

```bash
python svgexplode.py input.svg                  # writes input_exploded.svg
python svgexplode.py input.svg output.svg       # explicit output path
```

### Inkscape extension

**Option A: Add the repo as an extensions directory (recommended)**

1. Open Inkscape
2. Go to **Edit > Preferences > System**
3. Under **User extensions**, click the folder icon and add the path to this repo
4. Restart Inkscape

**Option B: Symlink into your existing extensions directory**

Linux:
```bash
ln -s /path/to/svgexplode/svgexplode.inx ~/.config/inkscape/extensions/
ln -s /path/to/svgexplode/svgexplode_ext.py ~/.config/inkscape/extensions/
ln -s /path/to/svgexplode/svgexplode.py ~/.config/inkscape/extensions/
```

macOS:
```bash
ln -s /path/to/svgexplode/svgexplode.inx ~/Library/Application\ Support/org.inkscape.Inkscape/config/inkscape/extensions/
ln -s /path/to/svgexplode/svgexplode_ext.py ~/Library/Application\ Support/org.inkscape.Inkscape/config/inkscape/extensions/
ln -s /path/to/svgexplode/svgexplode.py ~/Library/Application\ Support/org.inkscape.Inkscape/config/inkscape/extensions/
```

Windows (run as admin):
```cmd
mklink %APPDATA%\inkscape\extensions\svgexplode.inx C:\path\to\svgexplode\svgexplode.inx
mklink %APPDATA%\inkscape\extensions\svgexplode_ext.py C:\path\to\svgexplode\svgexplode_ext.py
mklink %APPDATA%\inkscape\extensions\svgexplode.py C:\path\to\svgexplode\svgexplode.py
```

Then restart Inkscape.

The extension appears under **Extensions > Generate from Path > SVG Explode - Extract Closed Regions**.

Options:
- **Spacing** between packed pieces (default 5mm)
- **Grid resolution** for collision detection (smaller = more accurate)
- **Rotation search** (none, 90deg, 45deg, or 15deg increments)
- **Pack pieces** toggle (uncheck to keep regions in their original positions)

Operates on selected paths, or all paths if nothing is selected.

## Examples

The repo includes example SVGs:

- `example.svg` - 3 intersecting paths forming 1 enclosed region
- `example2.svg` - 7 paths (curves + polylines) forming 12 enclosed regions

## Dependencies

- [svgpathtools](https://github.com/mathandy/svgpathtools) - SVG path parsing, bezier math, intersections
- [NumPy](https://numpy.org/) - geometry and raster operations
- [SciPy](https://scipy.org/) - binary dilation for spacing masks
