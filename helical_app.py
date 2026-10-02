"""
Program: helical_app.py
Author: Alexandre Prefontaine
Affiliation: University of Manitoba, Manitoba Centre for Proteomics and Systems Biology

Streamlit web application for Helical Plotter.

Run with:
    streamlit run helical_app.py
"""

import base64
from pathlib import Path

import streamlit as st

from helical_core import (
    GROUPS, MAX_LENGTH, MAX_PEPTIDES, OUTPUTS,
    png_bytes, prepare_inputs, read_csv_entries, save_grid, zip_bytes,
)

st.set_page_config(page_title="Helical Plotter", layout="wide")


@st.cache_data(show_spinner="Drawing figures…")
def cached_zip(peptides: tuple, selected: tuple) -> bytes:
    return zip_bytes(list(peptides), set(selected))


@st.cache_data(show_spinner=False)
def cached_preview(peptides: tuple) -> bytes:
    # lower resolution for the page; downloads are 200 dpi
    return png_bytes(save_grid, list(peptides), legend=True, dpi=100)


def shape_svg(shape: str | None, color: str) -> str:
    """Small outline symbol matching the figure shapes, for the color key."""
    if shape is None:
        return _as_img('<svg xmlns="http://www.w3.org/2000/svg" width="18" height="18"><text x="9" y="14" text-anchor="middle" '
                f'font-weight="bold" font-family="Arial" font-size="14" fill="{color}">P</text></svg>')
    attr = f'fill="white" stroke="{color}" stroke-width="2"'
    body = {
        "triangle": f'<polygon points="9,2 16.5,15.5 1.5,15.5" {attr}/>',
        "square":   f'<rect x="3" y="3" width="12" height="12" {attr}/>',
        "hexagon":  f'<polygon points="9,1.5 15.5,5.25 15.5,12.75 9,16.5 2.5,12.75 2.5,5.25" {attr}/>',
        "circle":   f'<circle cx="9" cy="9" r="6.5" {attr}/>',
        "dashed":   f'<circle cx="9" cy="9" r="6.5" {attr} stroke-dasharray="2.5,2"/>',
    }[shape]
    return _as_img(f'<svg xmlns="http://www.w3.org/2000/svg" width="18" height="18">{body}</svg>')


def _as_img(svg: str) -> str:
    # Streamlit strips inline <svg>, so embed it as an image
    return f'<img width="18" height="18" src="data:image/svg+xml;base64,{base64.b64encode(svg.encode()).decode()}">'


# White panel so colors read the same in light and dark mode (as in PepColors)
BOX = ("background:#FFFFFF;color:#222222;border:1px solid #DDDDDD;border-radius:8px;"
       "padding:12px 16px;display:inline-block;min-width:300px;font-family:sans-serif;"
       "font-size:15px")

st.title("Helical Plotter — helical wheel figures for peptides")
st.write(
    "Enter up to 16 peptides and download helical wheel projections (100° per residue) "
    "with residues colored and shaped by class, ready for publication."
)

# --- Color key ---
st.subheader("Color key")
rows = "".join(
    "<div style='margin:4px 0;display:flex;align-items:center;gap:10px'>"
    f"{shape_svg(g['shape'], g['color'])}"
    f"<span style='color:{g['color']};font-weight:bold;font-family:Courier New,monospace;"
    f"min-width:64px'>{g['residues']}</span>"
    f"<span>{g['label']}</span></div>"
    for g in GROUPS
)
st.html(f"<div style='{BOX}'>{rows}</div>")

st.divider()

# --- Input ---
st.subheader("1. Peptides")
source = st.radio("Input", ["Upload a CSV file", "Type sequences"], horizontal=True,
                  label_visibility="collapsed")

entries, stem = [], "helical_plotter"
if source == "Upload a CSV file":
    up = st.file_uploader("CSV with sequences in the first column", type=["csv"])
    if up is not None:
        try:
            entries = read_csv_entries(up)
            stem = Path(up.name).stem
        except Exception as e:
            st.error(f"Could not read the file: {e}")
else:
    text = st.text_area("One sequence per line", height=180,
                        placeholder="KHPDASVNFSEFSK\nYATLSLFNTYK\nATCIGNNSAAAVSM[ox]LK")
    entries = text.splitlines()

st.caption(
    f"Up to {MAX_PEPTIDES} peptides, each up to {MAX_LENGTH} residues "
    "(19–36 residues are drawn as a second ring). Sequences are converted to upper case. "
    "Mark a modification with a tag after the residue, e.g. `M[ox]` or `S(ph)`; it is "
    "drawn as `*`. Lower-case letters are **not** read as modifications."
)

peptides, notes = prepare_inputs(entries)

if notes["lowercase"]:
    st.warning(
        "Converted to upper case: " + ", ".join(notes["lowercase"]) +
        ". Lower-case letters are not treated as modifications; use a tag such as "
        "`M[ox]` to mark a modified residue."
    )
if notes["too_long"]:
    st.error(f"Not plotted (longer than {MAX_LENGTH} residues): " + ", ".join(notes["too_long"]))
if notes["dropped"]:
    st.warning(f"Only the first {MAX_PEPTIDES} peptides are used. Ignored: "
               + ", ".join(notes["dropped"]))

st.divider()

# --- Outputs ---
st.subheader("2. Files to include in the download")
selected = set()
for col, (key, label) in zip(st.columns(len(OUTPUTS)), OUTPUTS.items()):
    with col:
        if st.checkbox(label, value=(key == "matrix"), key=key):
            selected.add(key)
st.caption("“Bare” figures have no legend. Figures carry no sequence labels; the zip "
           "includes peptides.txt listing the order (matrix reads left to right).")

if peptides:
    if selected:
        st.download_button(
            "Download figures (.zip)",
            data=cached_zip(tuple(peptides), tuple(sorted(selected))),
            file_name=f"{stem}_helical.zip",
            mime="application/zip",
            type="primary",
        )
    else:
        st.info("Tick at least one box to build the download.")

    st.divider()

    # --- Preview: always the matrix with legend ---
    st.subheader("Preview")
    st.caption("Peptides in order, left to right: " +
               " · ".join(f"{k}. {p}" for k, p in enumerate(peptides, 1)))
    st.image(cached_preview(tuple(peptides)), width="stretch")

# --- Footer ---
st.divider()
col1, col2 = st.columns([1, 3])
with col1:
    here = Path(__file__).parent
    for logo in (here / "assets" / "UM-logo-horizontal-CMYK.jpg",
                 here / "UM-logo-horizontal-CMYK.jpg"):
        if logo.exists():
            st.image(str(logo), width=120)
            break
with col2:
    st.markdown(
        "**Alexandre Préfontaine** — Krokhin Laboratory<br>"
        "Manitoba Centre for Proteomics and Systems Biology<br>"
        "[University of Manitoba](https://umanitoba.ca) · Department of Internal Medicine",
        unsafe_allow_html=True,
    )
