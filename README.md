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
├── analysis/                     Analysis pipeline (three scripts)
│   ├── rebuild_processed_data.py     Rebuild processed tables from raw workbooks
│   ├── make_figures.py               Regenerate all main-text figures
│   └── build_supplement_tables.py    Rebuild the supplementary tables workbook
├── data/
│   ├── raw/                      Raw plate-reader workbooks (available on request)
│   │   ├── raw_screening_data/
│   │   └── raw_dose_response_data/
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
├── supplement/                   Supplementary tables (workbook and per-table CSVs)
└── manuscript/                   Manuscript document
```

## Reproducing the analysis

The processed tables, figures, and supplement are all included. To regenerate
them from source, run the three scripts in order from the repository root.

```bash
python analysis/rebuild_processed_data.py     # requires the raw workbooks
python analysis/make_figures.py
python analysis/build_supplement_tables.py
```

`rebuild_processed_data.py` reads the raw plate-reader workbooks from
`data/raw/`. Those workbooks are available from the authors on request and are
not redistributed here. `make_figures.py` and `build_supplement_tables.py` run
against the processed tables that ship in this repository, so figures and the
supplement can be regenerated without the raw data.

Requirements: Python 3.10 or later with numpy, pandas, scipy, scikit-learn,
matplotlib, openpyxl, and rdkit.

## Data and code license

Code is released under the MIT License (see LICENSE). Data are released under
CC BY 4.0.
