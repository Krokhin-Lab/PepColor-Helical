# Helical Plotter

Helical wheel figures for peptides, ready for publication. Residues are colored with the [PepColors](https://github.com/Krokhin-Lab/PepColors) residue-class scheme, and each class is drawn inside its own shape. This makes the amphipathic faces of a peptide easy to see, and the figures still read in grayscale.

![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)
![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)

**Launch the interactive app: [Helical Plotter on Streamlit](https://pepcolor-helical.streamlit.app/)**

**Authors:** Oleg V. Krokhin, Alexandre Préfontaine  
**Affiliation:** Manitoba Centre for Proteomics and Systems Biology, University of Manitoba  
**License:** GNU General Public License v3.0. See [LICENSE](LICENSE).

## What it does

1. Upload a CSV with peptide sequences in the first column, or type them in, one per line. You can enter up to **16 peptides**, each up to **36 residues** long.
2. Tick the outputs you want.
3. Download everything as a single `.zip`.

| Output | File(s) in the zip | Contents |
|---|---|---|
| **Individuals** | `individuals/01_<sequence>.png`, … | One wheel per peptide, with its sequence on top and the legend underneath |
| **Matrix** | `matrix.png` | All wheels in one figure, each with its sequence, and the legend at the bottom |
| **Individuals - bare** | `individuals_bare/01_<sequence>.png`, … | One wheel per peptide, wheel only |
| **Matrix - bare** | `matrix_bare.png` | All wheels in one figure, wheels only |
| **Separate Legend** | `legend.png` | The legend on its own |

Figures are 200 dpi PNGs on a white background.

The matrix reads left to right, then top to bottom. It uses 3 columns for up to 9 peptides and 4 columns for 10 to 16. Every zip also contains `peptides.txt`, which lists the peptides in order. Use it to identify the panels in the bare figures, since those have no sequence labels.

The app preview always shows the Matrix, with sequences and legend.

## How the wheel is drawn

- Residue *i* sits at *i* × 100° (3.6 residues per turn). Residue 1 is at the top, and the wheel runs clockwise.
- Every residue is numbered. A large **N** and **C** mark the termini.
- Positions repeat every 18 residues (5 turns). For peptides of 19 to 36 residues, residues 19–36 go on a second, outer ring. Longer peptides are refused.
- A faint line joins the residues in sequence order.

## Residue classes

| Residues | Class | Color | Shape |
|---|---|---|---|
| I L V M A | Aliphatic hydrophobic | Red | Triangle |
| W F Y | Aromatic hydrophobic | Violet red | Square |
| K R H | Basic | Yellow | Hexagon |
| E D | Acidic | Green | Hexagon |
| P | Proline | Black | None (letter only) |
| Q N S T C G | Polar / small | Blue | Circle |
| `*` | Modified / non-standard | Gray | Dashed circle |

Basic and acidic residues share the hexagon and are told apart by color. As in PepColors, the yellow is slightly darkened (`#E6B800`) so it stays readable on white.

## Modifications and input rules

- Sequences are converted to upper case. **Lower-case letters are not read as modifications.** The app warns you when it converts them.
- To mark a modified residue, put a tag right after it, e.g. `M[ox]`, `S(ph)` or `C[+57]`. That residue is drawn as `*`.
- Any other character outside the 20 standard residues (e.g. `X`, `B`, `Z`, `U`) is also drawn as `*`.
- A header row (`sequence`, `peptide`, …) and blank rows are ignored.
- If you enter more than 16 peptides, only the first 16 are used. Peptides longer than 36 residues are listed and not plotted.

`example_peptides.csv` contains nine test peptides you can use to try the app.

## Running locally

```bash
pip install -r requirements.txt
streamlit run helical_app.py
```

Command line, without the app:

```bash
python helical_core.py peptides.csv    # writes every output to peptides_helical.zip
```

## Files

| File | Purpose |
|---|---|
| `helical_core.py` | Input cleaning, wheel drawing, zip writer. No Streamlit dependency |
| `helical_app.py` | Streamlit front end |
| `example_peptides.csv` | Example input |
| `assets/` | University of Manitoba logo for the app footer |
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
