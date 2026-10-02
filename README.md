# Helical Plotter

Helical wheel figures for peptides, ready for publication. Residues are colored with the [PepColors](https://github.com/Krokhin-Lab/PepColors) residue-class scheme and drawn inside a shape for each class, so the amphipathic faces of a peptide stand out (and still read in grayscale).

![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)
![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)

**Launch the interactive app: [Helical Plotter on Streamlit](https://LINK-TO-ADD.streamlit.app/)**

**Authors:** Oleg V. Krokhin, Alexandre Préfontaine  
**Affiliation:** Manitoba Centre for Proteomics and Systems Biology, University of Manitoba  
**License:** GNU General Public License v3.0 — see [LICENSE](LICENSE).

## What it does

Upload a CSV with peptide sequences in the first column, or type them in (one per line). Up to **16 peptides**, each up to **36 residues**. Tick the outputs you want and download them as a single `.zip`:

| Output | File(s) in the zip |
|---|---|
| **Individuals** | `individuals/01_<sequence>.png`, … one wheel per peptide, with its sequence and the legend |
| **Matrix** | `matrix.png`, all wheels in one figure, each with its sequence, plus the legend |
| **Individuals - bare** | `individuals_bare/…`, one wheel per peptide, wheel only |
| **Matrix - bare** | `matrix_bare.png`, all wheels, wheels only |
| **Separate Legend** | `legend.png` |

Every zip also contains `peptides.txt`, listing the peptides in order. Individuals and Matrix show each sequence above its wheel with the legend underneath; the "bare" versions are the wheels only, for publication. The matrix reads left to right, top to bottom. The matrix uses 3 columns for up to 9 peptides and 4 columns for 10–16. Figures are 200 dpi PNG on a white background.

The app preview always shows the Matrix (sequences and legend).

## How the wheel is drawn

- Residue *i* sits at *i* × 100° (3.6 residues per turn), starting at the top and going clockwise. Every residue is numbered; a large **N** and **C** mark the termini.
- Positions repeat every 18 residues (5 turns). Peptides of 19–36 residues get a second, outer ring for residues 19–36. Longer peptides are refused.
- A faint line joins residues in sequence order.

## Residue classes

| Residues | Class | Color | Shape |
|---|---|---|---|
| I L V M A | Aliphatic hydrophobic | Red | Triangle |
| W F Y | Aromatic hydrophobic | Violet red | Square |
| K R H | Basic | Yellow | Hexagon |
| E D | Acidic | Green | Hexagon |
| P | Proline | Black | none |
| Q N S T C G | Polar / small | Blue | Circle |
| * | Modified / non-standard | Gray | Dashed circle |

## Modifications and input rules

- Sequences are converted to upper case. **Lower-case letters are not read as modifications**; the app warns when it converts them.
- Mark a modified residue with a tag right after it, e.g. `M[ox]`, `S(ph)`, `C[+57]`. That residue is drawn as `*`.
- Any other character that is not one of the 20 standard residues (e.g. `X`, `B`, `Z`, `U`) is also drawn as `*`.
- A header row (`sequence`, `peptide`, …) and blank rows are ignored.

## Running locally

```bash
pip install -r requirements.txt
streamlit run helical_app.py
```

Command line, without the app (writes every output to `peptides_helical.zip`):

```bash
python helical_core.py peptides.csv
```
