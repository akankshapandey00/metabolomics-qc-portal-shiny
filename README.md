# Metabolomics QC Portal (Shiny)

A lightweight QC + PCA portal for metabolomics **feature tables** (peak-area intensity matrices) with optional metadata.  
Built for core-facility style workflows: fast sanity checks, consistent exports, and a one-click HTML report that’s easy to share.

## What this app does
- Upload a **feature table** (CSV/TSV) and optional **metadata** (CSV/TSV)
- Detect or manually select **sample intensity columns**
- Validate sample IDs against metadata (**join check**)
- Generate:
  - QC by sample
  - QC by feature
  - PCA (select PCs + color by metadata column) + top loadings
- Export:
  - QC tables (CSV)
  - PCA scores (CSV)
  - Top loadings (CSV)
  - Filtered feature table (CSV)
  - Self-contained HTML QC report

---

## Project layout
```
metabolomics-qc-portal-shiny/
  app.R
  report_template.Rmd
  README.md
  data/
    tima_example_features.csv
    tima_example_metadata.tsv
  docs/
    screenshots/
  R/
    eda.R
    io.R
    qc.R
    report.R
```

> Note: The app currently runs from `app.R`. The `R/` folder is reserved for splitting logic into modules later.

---

## Input formats

### 1) Feature table (required)
A table where rows are **features** and columns include:
- descriptive columns (optional): `row ID`, `row m/z`, `row retention time`, etc.
- sample intensity columns (required): numeric peak areas

Supported:
- `.csv`, `.tsv`, `.txt`

#### Sample column detection
Two modes:
- **Auto (TIMA-style):** columns ending in ` Peak area`
- **Manual select:** you choose the intensity columns yourself

### 2) Metadata (optional but recommended)
A table with one row per sample. Must include a join column that matches sample IDs.

Common pattern:
- join column: `filename`
- other columns: `ATTRIBUTE_species`, `batch`, `run`, `group`, etc.

---

## Run locally

### Option A: Quick install (recommended for MVP)
Open R in the repo folder and run:

```r
install.packages(c(
  "shiny", "DT", "readr", "dplyr", "stringr", "tibble", "ggplot2",
  "rmarkdown", "knitr"
))

shiny::runApp()
```

### Option B: Reproducible setup with renv (recommended for sharing)
If you want everyone to get the same package versions:

```r
install.packages("renv")
renv::init()     # run once
renv::snapshot() # records package versions into renv.lock
```

On a fresh machine:

```r
install.packages("renv")
renv::restore()
shiny::runApp()
```

---

## Using the app (typical workflow)
1. **Load demo dataset** (button in sidebar) to confirm everything works
2. Upload your feature table + metadata
3. Choose sample intensity columns (Auto or Manual)
4. Check **Join check** tab to confirm sample IDs match metadata
5. Review **QC** and **PCA**
6. Go to **Exports** and download:
   - QC tables
   - PCA scores / loadings
   - Filtered feature table
   - HTML QC report

---

## Outputs
All downloads are named consistently using:
- dataset name (from uploaded features filename)
- date stamp (YYYY-MM-DD)

Exports include:
- `*_qc_by_sample_YYYY-MM-DD.csv`
- `*_qc_by_feature_YYYY-MM-DD.csv`
- `*_pca_scores_YYYY-MM-DD.csv`
- `*_pca_top_loadings_PC1_YYYY-MM-DD.csv`
- `*_filtered_features_miss20_nz20_YYYY-MM-DD.csv`
- `*_qc_report_YYYY-MM-DD.html`

---

## Notes on methods (high level)
- QC by sample: missing%, zero%, median intensity, median non-zero intensity, simple flags
- QC by feature: missing%, non-zero%, median intensity (preview shown in app; full table downloadable)
- PCA:
  - intensities -> median-impute missing values (per sample)
  - remove all-zero features and zero-variance features
  - log1p transform
  - PCA with centering + scaling
- HTML report:
  - self-contained HTML generated from `report_template.Rmd`
  - includes QC tables, PCA plot, top loadings, and session info

---

## Troubleshooting
- **No sample columns detected:** switch to “Manual select” and choose the intensity columns.
- **Join check shows many FALSE:** confirm the metadata join column (often `filename`) and whether sample columns include suffixes like ` Peak area`.
- **Report export fails:** install `rmarkdown` and `knitr`:
  ```r
  install.packages(c("rmarkdown", "knitr"))
  ```

---
