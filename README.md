# Predicting Language Lateralization from Resting-State fMRI Connectivity

Master's thesis (Data Science, UPC Barcelona)

This repository documents the analysis code used to produce the results in
the thesis: testing whether resting-state fMRI asymmetric functional
connectivity (AFC) predicts hemispheric language dominance in a
non-right-handed healthy adult cohort (n = 258, cohorts EZM2 and EZM3), with
ground-truth labels derived from task-based fMRI (Villar-Rodríguez et al.,
2024).

## Scope of this repository

This repository is provided for transparency and code review. It is **not**
packaged to run end-to-end from a fresh clone: subject-level fMRI data and
intermediate feature files are not included due to size and data-sharing
restrictions, and several scripts use file paths configured for the author's
local machine. Each folder below documents one stage of the analysis as it
was actually run.

## Repository structure

```
master-thesis/
├── matlab_scripts/
│   ├── 01_preprocessing/
│   ├── 02_functional_connectivity/
│   ├── 03_irreversability/          
│   └── 04_calculate_features/
├── R_scripts/
│   ├── 01_pca/
│   ├── 02_modeling/
│   ├── 03_backprojection/
│   ├── 04_visualization/
│   ├── irreversability/            
│   └── lookup/
├── python_scripts/
│   ├── analyze_demographics/
│   └── create_plots/
├── atlas/
│   ├── FINAL_left_SENS_HCP_Tian.nii
│   ├── FINAL_flipped_SENS_HCP_Tian.nii
│   ├── atlas_SENS_HCP_Tian_gm.nii
│   └── region_lookup.xlsx
├── requirements.txt
└── .gitignore
```

## Code organization

**`matlab_scripts/`** — fMRI preprocessing and connectivity computation
(CONN toolbox 25.b with SPM12, run separately per cohort).
- `01_preprocessing` — SLIC supervoxel clustering
- `02_functional_connectivity` — ROI signal
  extraction and FC/AFC matrix computation.
- `03_irreversability` — temporal irreversibility computation (INSIDEOUT
  framework).
- `04_calculate_features` — graph metric and feature export.

**`R_scripts/`** — statistical modeling and visualization.
- `01_pca` — PCA compression of the 9,900-feature AFC matrix.
- `02_modeling` — Strategy 1 (PCA + ridge logistic regression) and Strategy 2
  (interpretable features + ridge) classification models.
- `03_backprojection` — back-projection of ridge coefficients to AFC
  connections and ROIs.
- `04_visualization` — NIfTI map generation from back-projected weights and
  NodeAI coefficients.
- `irreversability` — IrrLI computation and group comparisons (independent
  analysis, does not feed the classification models above).
- `lookup` — ROI label lookup table used across scripts.

**`python_scripts/`** 
- `analyze_demographics` descriptive statistics or summary tables
  over the demographics of the dataset.
- `create_plots` — demographic plots.

## Atlas

Composite symmetric 200-ROI atlas (100 ROIs per hemisphere) combining:
- SENSAAS (Labache et al., 2019)
- HCP-MMP1.0 (Glasser et al., 2016)
- Tian subcortical parcellation (Tian et al., 2020)

## Requirements

- **R**: `glmnet`, `pROC`, `caret`, `emmeans`, `dplyr`, `RNifti`, `ggplot2`,
  `gridExtra`, `tidyverse`
- **MATLAB**: CONN toolbox 25.b, SPM12
- **Python**: see `requirements.txt`

## Data availability

Subject-level fMRI data and derived feature matrices are not included.
Demographic and outcome variables are generated locally from the EZM2/EZM3
cohorts.

## Notes
- Preprocessing was done using the CONN toolbox
- Atlas filenames in this folder (`FINAL_left_SENS_HCP_Tian.nii`, etc.) do
  not currently match the filenames read by the R visualization scripts
  (`left_composite_atlas.nii`, `region_lookup_100.csv`). One side needs to be
  updated to match the other.

## Citation

(Predicting Language Lateralization from Resting-State fMRI Connectivity, Victoria Margrethe Christensen, UPC, 2026)
