# Latent Chemical Recognition in an Evolved Periplasmic Binding Protein Family for Facile Diversification of Biosensors

Data, analysis code, and figures for the manuscript *Latent Chemical Recognition in an Evolved Periplasmic Binding Protein Family for Facile Diversification of Biosensors.*

Eighteen OpuBC/cpGFP fluorescent biosensors — seventeen from a single evolved
lineage plus one sequence-outgroup negative control — were screened against 63
chemically diverse ligands. Selected non-native hits were followed with
dose-response measurements, and response profiles were related to sequence
differences across the scaffold family.

## Quick start

```bash
pip install -r requirements.txt
python analysis/make_figures.py        # Figures 2-4 and analysis/tables/
```

`make_figures.py` regenerates the committed main-text figures from the processed
tables. Figure 1 is a composite of structural renders and is assembled separately
(see below).

## Repository layout

```
├── analysis/
│   ├── make_figures.py                Figures 2-4 from the processed tables
│   ├── compose_figure1.py             assembles Figure 1 from the rendered panels
│   ├── render_panel_a.py              Figure 1a  scaffold architecture      (PyMOL)
│   ├── render_panel_b.py              Figure 1b  sensing-mechanism schematic
│   ├── build_panel_c.py               Figure 1c  docked ciprofloxacin       (PyMOL)
│   ├── render_panel_d.py              Figure 1d  prediction vs crystal      (PyMOL)
│   ├── build_supplement_tables.py     rebuilds the supplementary workbook
│   └── tables/                        per-figure source tables (written by make_figures.py)
├── data/
│   ├── processed/
│   │   ├── response_matrix.csv        63 x 18 single-concentration dF/F0 matrix
│   │   ├── response_sem_matrix.csv    matched propagated SEM matrix
│   │   ├── dose_response_fits.csv     Hill-fit parameters for 49 curves
│   │   ├── dose_response_exclusions.csv
│   │   └── dose_responses/            49 per-curve CSVs
│   ├── metadata/
│   │   ├── ligand_smiles.csv          ligand names and SMILES
│   │   ├── ligand_categories.csv      ligand scope classes
│   │   └── sensor_sequences.csv       sensor sequences and lineage labels
│   └── structures/                    7S7U crystal, ESMFold prediction, docked poses
├── figures/                           main-text figures (PNG at 600 dpi, vector PDF)
│   └── panels/                        individual Figure 1 panels
├── supplement/                        Supplementary Tables S1-S7 (CSV + XLSX)
└── manuscript/                        manuscript text
```

## Dataset

The primary dataset is a 63-ligand x 18-sensor matrix of single-concentration
ΔF/F₀ responses (1,134 ligand-sensor pairs) with a matched propagated SEM matrix.

| Quantity | Value |
|----------|-------|
| Ligands | 63 |
| Sensors | 18 (17 scaffold + 1 outgroup control) |
| Ligand-sensor pairs | 1,134 |
| Pairs with ΔF/F₀ > 0.3 | 124 |
| Pairs with ΔF/F₀ > 1.0 | 50 |
| Ligands with ≥1 response > 0.3 | 24 |
| Ligands with ≥1 response > 1.0 | 8 |
| Response range | −0.908 to 12.144 |

The 17 scaffold sequences share a common analyzed length of 521 residues and vary
at 31 positions; pairwise distances span 1-26 substitutions (95.0-99.8% identity).

## Conventions

- **Response thresholds.** ΔF/F₀ > 0.3 is a permissive engineering hit; ΔF/F₀ > 1.0
  is a strong response. Ligand hit breadth counts responding sensors; sensor hit
  breadth counts detected ligands.
- **Outgroup control.** Fent2 436L shares 8.7% identity with the scaffold family.
  It is retained in the chemical screen as a negative control and excluded from all
  scaffold sequence statistics.
- **Excluded measurements.** Ciprofloxacin well H2 at 200 µM was not dispensed by
  the liquid handler; that concentration is recomputed from the remaining two
  replicates. Four weak or negative curves are excluded from the processed set
  (MEHP/cc93, aspartame/iLevaSnFR, bilirubin/V4.8.1.2, theobromine/Tap1.0). All
  exclusions are listed in `data/processed/dose_response_exclusions.csv`.
- **Unsaturated dose responses.** Curves without a plateau in the measured range
  are treated as lower bounds. Their fitted EC50 values are numerical
  extrapolations and are not interpreted as potencies; figures mark them with
  dashed lines and arrowheads.
- **Apparent sensitivity.** EC50 and ΔFmax/EC50 derived from these fits combine
  ligand binding, PBP conformational change, coupling to cpGFP, and fluorescence
  output. They are not binding affinities; no orthogonal binding assay was
  performed.

## Reproducing the figures and supplement

Figures 2-4 regenerate from the processed tables:

```bash
python analysis/make_figures.py
```

Figure 1 is composed from four rendered panels. The panel PNGs are committed under
`figures/panels/`, so the composite can be rebuilt without PyMOL:

```bash
python analysis/compose_figure1.py
```

Panel (b) is a matplotlib schematic and needs no extra dependencies:

```bash
python analysis/render_panel_b.py
```

Regenerating the three structural panels requires PyMOL:

```bash
cd analysis
python render_panel_a.py     # scaffold architecture
python build_panel_c.py      # docked ciprofloxacin
python render_panel_d.py     # ESMFold prediction vs crystal
```

Panel-specific details, including the choice of crystal structure and the
provenance of the docked pose, are documented in `figures/FIGURE1_PANELS.md`.

The supplementary tables and workbook are built separately:

```bash
python analysis/build_supplement_tables.py
```

Running all seven scripts on a clean checkout reproduces every generated file
byte for byte: 8 main-text figure PNGs and PDFs, 5 Figure 1 panels, 4 analysis
tables, and the 9 supplement files, for 26 in total. Embedded timestamps are
pinned in the PDF and workbook writers, and PyMOL ray tracing is pinned to one
thread, so rebuilt outputs can be compared directly against the committed
copies. `figures/FIGURE1_PANELS.md` is written documentation, not a generated
file, and is not part of that count.

## Requirements

Python 3.11 or newer; see `requirements.txt`. The PyMOL panel scripts additionally
require `pymol-open-source`, which is not needed for the main figure pipeline.

## Data provenance

Raw plate-reader workbooks are not included in this deposition. The processed
tables under `data/processed/` are the analysis inputs and are the source of every
number in the manuscript figures.

## Licence

Code is released under the MIT License. Data are released under CC BY 4.0.
