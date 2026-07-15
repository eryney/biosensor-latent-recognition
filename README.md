# Latent Chemical Recognition Diversity in an Evolved Periplasmic Binding Protein Biosensor Family

Data, analysis code, and figures for the manuscript "Latent Chemical Recognition
Diversity in an Evolved Periplasmic Binding Protein Biosensor Family."

This repository screens 18 related OpuBC/cpGFP fluorescent biosensors against 63
structurally diverse ligands, follows selected non-native hits with dose-response
measurements, and relates response profiles to sequence differences across the
scaffold family.

## Repository layout

```
biosensor-latent-recognition/
├── analysis/                     Analysis pipeline (two scripts)
│   ├── make_figures.py               Regenerate all main-text figures
│   └── build_supplement_tables.py    Rebuild the supplementary tables workbook
├── data/
│   ├── processed/               Processed response matrices and dose-response curves
│   │   ├── response_matrix.csv       63 x 18 single-concentration dF/F0 matrix
│   │   ├── response_sem_matrix.csv   Matched propagated SEM matrix
│   │   ├── dose_response_fits.csv    Hill-fit parameters for 49 curves
│   │   ├── dose_response_exclusions.csv
│   │   └── dose_responses/           49 per-curve CSVs
│   └── metadata/
│       ├── ligand_smiles.csv         Ligand names and SMILES
│       ├── ligand_categories.csv     Ligand scope classes
│       └── sensor_sequences.csv      Sensor sequences and lineage labels
├── figures/                      Main-text figures (PNG and PDF)
├── structures/                   Structure models behind Figure 1
│   ├── cipro_docked_7S7T.sdf         Docked ciprofloxacin pose (PDB 7S7T)
│   ├── cipro_docked_7S7U.sdf         Docked ciprofloxacin pose (PDB 7S7U)
│   └── esmfold_cc93_fulllength.pdb   ESMFold model of the cc93 scaffold
└── supplement/                   Supplementary tables (workbook and per-table CSVs)
```

## Reproducing the analysis

The processed tables, figures, and supplement are all included. To regenerate
the figures and supplement from the processed data, run the two scripts in order
from the repository root.

```bash
python analysis/make_figures.py               # regenerates all main-text figures
python analysis/build_supplement_tables.py     # rebuilds the supplement workbook
```

Run `make_figures.py` first. It reads the processed tables under `data/processed/`
and the metadata under `data/metadata/`, writes the main-text figures to
`figures/`, and writes intermediate summary tables (for example
`figure2_ligand_scope_summary.csv`) to `analysis/`. `build_supplement_tables.py`
then reads those summary tables together with the processed data to assemble
`supplement/Supplementary_Tables.xlsx`.

The processed tables that ship in this repository were derived from raw
plate-reader workbooks. Those raw workbooks are available from the authors on
request and are not redistributed here.

Requirements: Python 3.10 or later with numpy, pandas, scipy, scikit-learn,
matplotlib, openpyxl, and rdkit.

## Data and code license

Code is released under the MIT License (see LICENSE). Data are released under
CC BY 4.0.
