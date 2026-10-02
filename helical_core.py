# Helical Plotter — helical wheel figures for peptides
# Copyright (C) 2026 Oleg V. Krokhin, Alexandre Préfontaine
# University of Manitoba, Manitoba Centre for Proteomics and Systems Biology
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
"""
Program: helical_core.py
Author: Alexandre Prefontaine
Affiliation: University of Manitoba, Manitoba Centre for Proteomics and Systems Biology

Core module for Helical Plotter: input cleaning, helical wheel drawing and the zip
writer. No Streamlit dependency, so it can be used from scripts:

    from helical_core import prepare_inputs, zip_bytes, OUTPUTS
    peptides, notes = prepare_inputs(["KHPDASVNFSEFSK", "ATCIGNNSAAAVSM[ox]LK"])
    open("wheels.zip", "wb").write(zip_bytes(peptides, set(OUTPUTS)))

Geometry
--------
Residue i sits at i * 100 degrees (alpha-helix periodicity, 3.6 residues per turn),
starting at the top and going clockwise. Positions repeat every 18 residues (5 turns),
so peptides of 19-36 residues are drawn as two concentric rings (1-18 inside,
19-36 outside). Longer peptides are refused.

Residue classes
---------------
Colors follow the PepColors "Residue class" scheme; each class also has a shape so
the figures still read in grayscale. Anything that is not one of the 20 standard
residues (a tagged modification such as M[ox], or a letter like X) is drawn as '*'.
"""

import io
import math
import re
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, RegularPolygon

# --- Limits ---
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
DEGREES_PER_RESIDUE = 100
PER_RING = 18          # 18 * 100 deg = 5 full turns, then positions repeat
MAX_PEPTIDES = 16
MAX_LENGTH = 36        # two rings
DPI = 200

# --- Residue classes: label, residues, color (PepColors), shape ---
GROUPS = [
    {"label": "Aliphatic hydrophobic",   "residues": "ILVMA",  "color": "#D62728", "shape": "triangle"},
    {"label": "Aromatic hydrophobic",    "residues": "WFY",    "color": "#C71585", "shape": "square"},
    {"label": "Basic",                   "residues": "KRH",    "color": "#E6B800", "shape": "hexagon"},
    {"label": "Acidic",                  "residues": "ED",     "color": "#2CA02C", "shape": "hexagon"},
    {"label": "Proline",                 "residues": "P",      "color": "#000000", "shape": None},
    {"label": "Polar / small",           "residues": "QNSTCG", "color": "#1F5FD8", "shape": "circle"},
    {"label": "Modified / non-standard", "residues": "*",      "color": "#7F7F7F", "shape": "dashed"},
]
RESIDUE = {aa: g for g in GROUPS for aa in g["residues"]}

# --- Wheel layout ---
R = 1.5             # inner ring radius (residues 20 degrees apart do not collide)
RING_STEP = 0.95    # gap between concentric rings
RING_GREY = "#B8B8B8"
PATH_GREY = "#E0E0E0"
NUMBER_GREY = "#666666"
TERMINUS = "#111111"

# --- Header names to skip in a CSV ---
KNOWN_HEADERS = {"seq", "sequence", "sequences", "peptide", "peptides", "pep"}

# --- Outputs offered in the app: key -> label ---
OUTPUTS = {
    "individuals":      "Individuals",
    "matrix":           "Matrix",
    "individuals_bare": "Individuals - bare",
    "matrix_bare":      "Matrix - bare",
    "legend":           "Separate Legend",
}

_TAG = re.compile(r"[\[\(][^\]\)]*[\]\)]")


# =============================================================== input
def parse_sequence(seq: str) -> list[str]:
    """
    Split an (upper-case) sequence into residues. A bracketed tag right after a
    residue (M[OX], S(PH), C[+57]) makes that residue one '*'. Any other character
    outside the 20 standard residues is also '*'.
    """
    s = re.sub(r"[A-Z]?" + _TAG.pattern, "*", str(seq).strip().upper().replace(" ", ""))
    return [c if c in VALID_AA else "*" for c in s]


def prepare_inputs(entries) -> tuple[list[str], dict]:
    """
    Clean raw entries (CSV column or typed lines).

    Returns (peptides, notes). Peptides are upper-case strings ready to plot.
    notes holds user-facing problems:
        lowercase - entries with lower-case letters (converted; NOT read as modifications)
        too_long  - entries over MAX_LENGTH residues (refused)
        dropped   - entries beyond MAX_PEPTIDES (ignored)
    """
    notes = {"lowercase": [], "too_long": [], "dropped": []}
    peptides = []
    for e in entries:
        raw = str(e).strip().replace(" ", "")
        if not raw or raw.lower() in KNOWN_HEADERS or raw.lower() == "nan":
            continue
        if any(c.islower() for c in _TAG.sub("", raw)):
            notes["lowercase"].append(raw)
        seq = raw.upper()
        if len(parse_sequence(seq)) > MAX_LENGTH:
            notes["too_long"].append(seq)
            continue
        peptides.append(seq)
    if len(peptides) > MAX_PEPTIDES:
        notes["dropped"] = peptides[MAX_PEPTIDES:]
        peptides = peptides[:MAX_PEPTIDES]
    return peptides, notes


def read_csv_entries(source) -> list[str]:
    """First column of a CSV (path or file-like), read without assuming a header."""
    import pandas as pd
    df = pd.read_csv(source, header=None, dtype=str, skip_blank_lines=True)
    return df.iloc[:, 0].dropna().tolist()


# =============================================================== drawing
def _shape_patch(shape, xy, color):
    kw = dict(facecolor="white", edgecolor=color, linewidth=2.2, zorder=3)
    if shape == "triangle":
        return RegularPolygon(xy, 3, radius=0.25, orientation=0, **kw)
    if shape == "square":
        return RegularPolygon(xy, 4, radius=0.21, orientation=math.pi / 4, **kw)
    if shape == "hexagon":
        return RegularPolygon(xy, 6, radius=0.20, orientation=math.pi / 6, **kw)
    if shape == "dashed":
        return Circle(xy, 0.17, linestyle=(0, (2, 1.5)), **{**kw, "linewidth": 1.8})
    if shape == "circle":
        return Circle(xy, 0.17, **kw)
    # proline: no border, just mask the connecting line behind the letter
    return Circle(xy, 0.16, facecolor="white", edgecolor="none", zorder=3)


def _position(i):
    ang = math.radians(90 - DEGREES_PER_RESIDUE * i)
    r = R + RING_STEP * (i // PER_RING)
    return r * math.cos(ang), r * math.sin(ang), ang, r


def wheel_extent(n: int) -> float:
    """Half-width of the axes needed for a wheel of n residues."""
    rings = max(1, math.ceil(n / PER_RING))
    return R + RING_STEP * (rings - 1) + 0.95


def draw_wheel(ax, residues: list[str], extent: float | None = None):
    """Draw one helical wheel (no title; callers add the sequence when wanted)."""
    n = len(residues)
    rings = max(1, math.ceil(n / PER_RING))
    pts = [_position(i) for i in range(n)]
    lim = extent or wheel_extent(n)

    ax.set_aspect("equal")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.axis("off")

    for k in range(rings):
        ax.add_patch(Circle((0, 0), R + RING_STEP * k, fill=False, edgecolor=RING_GREY,
                            linestyle=(0, (3, 3)), linewidth=1.0, zorder=1))
    ax.plot([p[0] for p in pts], [p[1] for p in pts], color=PATH_GREY, linewidth=0.9, zorder=2)

    for i, (aa, (x, y, ang, r)) in enumerate(zip(residues, pts)):
        g = RESIDUE[aa]
        ax.add_patch(_shape_patch(g["shape"], (x, y), g["color"]))
        star = aa == "*"
        ax.text(x, y - (0.07 if star else 0), aa, ha="center", va="center",
                fontsize=17 if star else 13, fontweight="bold", color=g["color"], zorder=4)

        # position number: outside the outermost ring, inside for an inner ring
        on_outer = (i // PER_RING) == rings - 1
        rn = r + 0.40 if on_outer else r - 0.36
        ax.text(rn * math.cos(ang), rn * math.sin(ang), str(i + 1), ha="center", va="center",
                fontsize=8 if rings == 1 else 7.5, color=NUMBER_GREY, zorder=4)

    # termini: large bold N / C one step beyond the position number
    for i, tag in ((0, "N"), (n - 1, "C")):
        if tag == "C" and n < 2:
            continue
        x, y, ang, r = pts[i]
        on_outer = (i // PER_RING) == rings - 1
        rt = r + 0.70 if on_outer else r - 0.68
        ax.text(rt * math.cos(ang), rt * math.sin(ang), tag, ha="center", va="center",
                fontsize=16, fontweight="bold", color=TERMINUS, zorder=4)


# =============================================================== legend
def _legend_handles(scale=1.0):
    def marker(m, ms, g):
        return Line2D([], [], linestyle="", marker=m, markersize=ms * scale,
                      markerfacecolor="white", markeredgecolor=g["color"],
                      markeredgewidth=1.8 * scale, label=g["label"])
    glyph = {"triangle": ("^", 11), "square": ("s", 9), "hexagon": ("h", 11),
             "circle": ("o", 9)}
    handles = []
    for g in GROUPS:
        if g["shape"] in glyph:
            handles.append(marker(*glyph[g["shape"]], g))
        else:  # proline letter, or the modification star
            m = r"$\mathbf{P}$" if g["shape"] is None else r"$\ast$"
            handles.append(Line2D([], [], linestyle="", marker=m,
                                  markersize=(9 if g["shape"] is None else 11) * scale,
                                  color=g["color"], label=g["label"]))
    return handles


def _add_legend(fig, y, ncol, scale=1.0):
    fig.legend(handles=_legend_handles(scale), loc="lower center", ncol=ncol, frameon=False,
               fontsize=9 * scale, bbox_to_anchor=(0.5, y), handletextpad=0.3,
               columnspacing=1.2)


def save_legend(out, dpi=DPI):
    """The legend on its own, one row."""
    fig = plt.figure(figsize=(11, 0.55))
    _add_legend(fig, 0.05, ncol=len(GROUPS))
    fig.savefig(out, dpi=dpi, facecolor="white")
    plt.close(fig)


# =============================================================== figures
TITLE_H = 0.45      # inches reserved above each wheel for its sequence


def _title(ax, peptide):
    ax.set_title(peptide, fontsize=12, fontfamily="monospace", pad=2)


def save_individual(peptide: str, out, bare=False, dpi=DPI):
    """One wheel. Default: sequence on top and legend underneath. bare: the wheel only."""
    res = parse_sequence(peptide)
    rings = max(1, math.ceil(len(res) / PER_RING))
    side = 5.0 + 1.8 * (rings - 1)
    t_h = 0.0 if bare else TITLE_H
    leg_h = 0.0 if bare else 0.9
    height = side + t_h + leg_h
    fig, ax = plt.subplots(figsize=(side, height))
    draw_wheel(ax, res)
    fig.subplots_adjust(bottom=leg_h / height, top=1 - t_h / height, left=0.0, right=1.0)
    if not bare:
        _title(ax, peptide)
        _add_legend(fig, 0.01, ncol=3)
    fig.savefig(out, dpi=dpi, facecolor="white")
    plt.close(fig)


def grid_shape(n: int) -> tuple[int, int]:
    """Up to 9 peptides: 3 columns. 10 or more: 4 columns."""
    cols = min(n, 3 if n <= 9 else 4)
    return math.ceil(n / cols), cols


def save_grid(peptides: list[str], out, bare=False, dpi=DPI):
    """
    All wheels in a matrix, read left to right in input order.
    Default: each wheel has its sequence on top, legend at the bottom. bare: wheels only.
    """
    n = len(peptides)
    rows, cols = grid_shape(n)
    parsed = [parse_sequence(p) for p in peptides]
    extent = max(wheel_extent(len(p)) for p in parsed)   # same scale in every panel
    cell = 4.8 * extent / wheel_extent(1)

    scale = max(1.0, cols / 3)                 # bigger legend for wider grids
    ncol = len(GROUPS) if cols >= 3 else (4 if cols == 2 else 3)
    t_h = 0.0 if bare else TITLE_H
    leg_h = 0.0 if bare else (0.6 if ncol == len(GROUPS) else 0.9) * scale
    gap = 0.05                                 # inches between panels
    height = rows * (cell + t_h) + (rows - 1) * gap + leg_h

    fig, axes = plt.subplots(rows, cols, figsize=(cell * cols, height), squeeze=False)
    for k, ax in enumerate(axes.flat):
        if k < n:
            draw_wheel(ax, parsed[k], extent=extent)
            if not bare:
                _title(ax, peptides[k])
        else:
            ax.axis("off")
    fig.subplots_adjust(bottom=leg_h / height, top=1 - t_h / height, left=0.005, right=0.995,
                        hspace=(t_h + gap) / cell, wspace=0.04)
    if not bare:
        _add_legend(fig, 0.005, ncol=ncol, scale=scale)
    fig.savefig(out, dpi=dpi, facecolor="white")
    plt.close(fig)


# =============================================================== outputs
def png_bytes(func, *args, **kw) -> bytes:
    buf = io.BytesIO()
    func(*args, buf, **kw)
    return buf.getvalue()


def _safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_+\-]", "_", name)


def zip_bytes(peptides: list[str], selected: set) -> bytes:
    """One zip holding every selected output (keys of OUTPUTS)."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        if "matrix" in selected:
            z.writestr("matrix.png", png_bytes(save_grid, peptides))
        if "matrix_bare" in selected:
            z.writestr("matrix_bare.png", png_bytes(save_grid, peptides, bare=True))
        for key, folder, bare in (("individuals", "individuals", False),
                                  ("individuals_bare", "individuals_bare", True)):
            if key in selected:
                for k, p in enumerate(peptides, 1):
                    z.writestr(f"{folder}/{k:02d}_{_safe(p)}.png",
                               png_bytes(save_individual, p, bare=bare))
        if "legend" in selected:
            z.writestr("legend.png", png_bytes(save_legend))
        # order of the panels (bare figures carry no sequence labels)
        z.writestr("peptides.txt", "".join(f"{k}\t{p}\n" for k, p in enumerate(peptides, 1)))
    return buf.getvalue()


if __name__ == "__main__":
    import sys
    from pathlib import Path

    if len(sys.argv) < 2:
        sys.exit("usage: python helical_core.py peptides.csv")
    src = Path(sys.argv[1])
    peptides, notes = prepare_inputs(read_csv_entries(src))
    for key, items in notes.items():
        if items:
            print(f"{key}: {', '.join(items)}")
    if not peptides:
        sys.exit("No peptides to plot.")
    out = src.with_name(src.stem + "_helical.zip")
    out.write_bytes(zip_bytes(peptides, set(OUTPUTS)))
    print(f"Wrote {out} ({len(peptides)} peptides)")
