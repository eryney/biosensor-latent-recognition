# Latent Chemical Recognition Diversity in an Evolved PBP Biosensor Family

Data and analysis code for the manuscript:

> **Latent Chemical Recognition Diversity in an Evolved Periplasmic Binding
> Protein Biosensor Family.**
> Eryney Marrogi, Anand Muthusamy, Aaron Nichols, Henry Lester.

This study screens a family of 18 closely related OpuBC/cpGFP fluorescent
biosensors — evolved for nicotine, nicotinic agonists, SSRIs, opioids, and
related drug-like ligands — against 63 chemically diverse compounds, then
follows selected non-native hits with dose-response measurements and
sequence-function analysis. The central finding is that a small family of
highly similar proteins contains surprisingly diverse ligand-recognition
profiles, and that broad empirical screening can expose practical starting
points for new sensors. Ciprofloxacin is the strongest near-term application
lead.

## Repository layout

```
data/
  processed/
    response_matrix.csv            63 ligand x 18 sensor single-concentration ΔF/F0 matrix
    response_sem_matrix.csv        matched propagated SEM matrix
    ligand_smiles.csv              ligand names + SMILES
    ligand_categories.csv          paper-facing scope classes
    dose_response_fits.csv         Hill-fit parameters, fit quality, exclusions
    dose_response_exclusions.csv   excluded curves/wells + reasons
    dose_responses/                49 processed analyte-sensor dose-response CSVs
  supplement/
    Supplementary_Tables_v3.xlsx   full supplementary workbook
analysis/
  rebuild_v3_processed_from_raw.py builds processed tables from raw workbooks
  make_v3_figures.py               generates all main figures
  build_v3_supplement_tables.py    builds the supplementary workbook
  build_biorxiv_docx.py            assembles the submission .docx from the manuscript
figures/
  Figure1.{png,pdf}                    scaffold structure, mechanism, pocket, prediction overlay
  fig1_screen_scope_heatmap.{png,pdf}  Figure 2: compound screen grouped by ligand scope
  fig2_sequence_activity_cliffs.{png,pdf}  Figure 3: sequence-function cliffs / lineage
  fig3_dose_response_leads.{png,pdf}   Figure 4: dose-response leads
  fig4_translational_relevance.{png,pdf}  Figure 5: application gaps and 1-mutation gains
structures/
  esmfold_cc93_fulllength.pdb      full-length ESMFold prediction of cc93/iNicSnFR3a
  cipro_docked_7S7U.sdf            illustrative docked ciprofloxacin pose (7S7U frame)
  cipro_docked_7S7T.sdf            illustrative docked ciprofloxacin pose (7S7T frame)
manuscript/
  manuscript_biorxiv.md            manuscript source (Markdown)
  manuscript_biorxiv.docx          bioRxiv submission document (figures embedded)
```

> **Note on figure file names.** The `figN_` prefixes are historical. In the
> current manuscript numbering the structure figure is **Figure 1**
> (`Figure1.*`), and the four data figures are **Figures 2–5**
> (`fig1_`→Fig 2, `fig2_`→Fig 3, `fig3_`→Fig 4, `fig4_`→Fig 5).

## Reproducing the analysis

```bash
python -m pip install -r requirements.txt
# Rebuild processed tables from raw workbooks (raw data not included here; see below)
python analysis/rebuild_v3_processed_from_raw.py
# Regenerate all figures from processed data
python analysis/make_v3_figures.py
# Rebuild the supplementary workbook
python analysis/build_v3_supplement_tables.py
```

The figure and supplement scripts run from the **processed** CSVs included in
this repository, so they reproduce the paper figures without the raw data.

## Data scope

This repository contains **processed** data and analysis code. The raw
plate-reader workbooks (single-concentration screen and dose-response plates)
are large binary Excel files and are not included here; they are available from
the authors on request or via a data archive (to be deposited). The
`rebuild_v3_processed_from_raw.py` script documents exactly how the processed
tables are derived from those workbooks.

## Key definitions

- **Meaningful / permissive engineering hit:** ΔF/F0 > 0.3
- **Strong response:** ΔF/F0 > 1.0
- Fluorescence EC50 and ΔFmax/EC50 are reported as **apparent sensitivity**;
  they combine ligand binding, PBP conformational switching, cpGFP coupling, and
  fluorescence output. Direct binding affinity is not measured here.
- The ciprofloxacin H2 well at 200 µM is a known bad well; it is excluded and
  the 200 µM point recomputed from the remaining 2 replicates (documented in
  `dose_response_exclusions.csv`).

## License

Code is released under the MIT License (see `LICENSE`). Data files
(`data/`, `figures/`, `structures/`) are released under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Citation

See `CITATION.cff`.
