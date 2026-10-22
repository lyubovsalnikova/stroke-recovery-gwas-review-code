# From Genetic Discovery to Clinical Utility in Stroke Recovery: code and reporting checklist

Code, small input tables and the reporting checklist for the review *From Genetic Discovery to Clinical Utility in Stroke
Recovery: Collider Bias, Phenotype Coupling, and Data Limitations* (Frontiers in Neurology).

Every number in the review that was computed (rather than taken from the literature) is produced by a script here.
Reference outputs are in `results/`, so each number can be checked without re-running anything.

## What reproduces what

| Review item | Script | Inputs | Output (`results/`) | Runs as shipped? |
|---|---|---|---|---|
| Table 1, Figure 2, Supplementary Note S3 (change-score and collider simulation; closed forms checked by Monte Carlo, 2 x 10^6 draws per cell) | `simulation/changescore_rtm_simulation.py` | `data/T7_recovery_index_feasibility.csv` (identity check only) | `tables/A_changescore_bias.csv`, `B_identity_check.csv`, `C_collider_induced.csv`; `figures/changescore_rtm_simulation.png/.pdf` | yes |
| Figure 1 (three constraints; panel B simulated) | `figures/constraints_figure.py` | none | `figures/constraints.png/.pdf` | yes |
| Figure 3 (phase, compartment, construct) | `figures/multiomic_mismatch_figure.py` | `data/phase_coverage.csv` | `figures/multiomic_mismatch.png/.pdf` | yes |
| Table 3 (GWAS Catalog associations) | `catalog/table3_catalog_summary.py` | `data/gwas_catalog_stroke_outcome_severity.csv` | `tables/table3_catalogued_outcome_associations.csv` | yes |
| Table 6 (reporting standard) | - | - | `reporting_checklist/` (Word and CSV, to fill in) | - |
| Supplementary Note S1 (GSE296792 expression indices); sampling-time audit | `derivations/01_recovery_index_tables.py` | companion outputs (PROJECT_ROOT) | `tables/T7_recovery_index_feasibility.csv`, `T9_longitudinal_design_audit.csv` | needs PROJECT_ROOT |
| Supplementary Note S2 (GSE163256 worked example) | `derivations/02_gse163256_worked_example.py` | companion outputs; `data/glycolysis_HIF_panel.csv` | `tables/gse163256_worked_example.csv` | needs PROJECT_ROOT |
| Figure 3A, Supplementary Table S3 (sampling-time coverage) | `derivations/03_phase_coverage.py` (after 01) | companion outputs; `data/comparator_sampling_lag.csv` | `tables/phase_coverage.csv` | needs PROJECT_ROOT |
| Section 7 (plasma proteome vs blood transcriptome measured space) | `derivations/04_proteome_overlap.py` | companion outputs; PRIDE PXD015006 `proteinGroups.txt` | `tables/T11_proteome_overlap_feasibility.csv` | needs PROJECT_ROOT and PXD015006_DIR (proteinGroups.txt, 46.7 MB, from the PRIDE FTP archive) |

## Quick start

```
pip install -r requirements.txt
python simulation/changescore_rtm_simulation.py
python figures/constraints_figure.py
python figures/multiomic_mismatch_figure.py
python catalog/table3_catalog_summary.py
```

Each script writes into `results/` and overwrites the reference copy; `git diff` then shows any difference.

## Derivations

The four scripts in `derivations/` compute the inputs that come from the authors' integrative analysis of public
stroke blood transcriptomes. They read that analysis's outputs from the directory given by `PROJECT_ROOT`, laid out
as the companion repository `stroke-transcriptomics-code` writes them (`harmonization/metadata/`, `diff_expr/`,
`cellcomp/tables/`, `sorted_cell/tables/`, `trajectories/temporal_programs/`, `GSE163256/`).
`04_proteome_overlap.py` also needs `proteinGroups.txt` of the PRIDE project PXD015006 in the directory given by
`PXD015006_DIR`.

```
export PROJECT_ROOT=/path/to/stroke-transcriptomics-code/outputs
python derivations/01_recovery_index_tables.py
python derivations/02_gse163256_worked_example.py
python derivations/03_phase_coverage.py
export PXD015006_DIR=/path/to/PXD015006
python derivations/04_proteome_overlap.py
```

## Data sources

All primary data are public. Blood transcriptomes: GEO GSE16561, GSE22255, GSE37587, GSE58294, GSE73378, GSE125512,
GSE158312, GSE296792 and GSE163256 (processing code in the companion repository). Plasma proteome: PRIDE PXD015006.
GWAS associations: NHGRI-EBI GWAS Catalog, trait "stroke outcome severity measurement" (78 associations from three
studies). `data/glycolysis_HIF_panel.csv` and `data/comparator_sampling_lag.csv` are taken from the supplementary
tables of the companion study (Kashatnikova et al., 2026).

## Software

Python 3.14 with numpy, pandas, scipy and matplotlib (`requirements.txt`). Figures use Arial when it is installed.
All source files are ASCII.

## License and citation

Code: MIT (`LICENSE`). Please cite the article (`CITATION.cff`).
