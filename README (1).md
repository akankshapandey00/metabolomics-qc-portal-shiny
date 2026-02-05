# Metabolomics QC Portal (Shiny)

A lightweight, facility-style QC + PCA web app for metabolomics **feature tables**.  
Designed for fast triage: validate uploads, spot batch/run effects, and export clean outputs for downstream analysis.

## What you can do
- **Upload** a feature table (CSV/TSV) and optional metadata (CSV/TSV)
- **Auto-detect** TIMA-style intensity columns ending with ` Peak area`, or **manually select** columns
- **Join check** to confirm sample IDs match metadata (default join key: `filename`)
- **QC by sample** and **QC by feature**
- **PCA** with variance explained + top loadings
- **Batch QC**: group summary table + signal distribution plot by any metadata column
- **Exports**: CSVs + an **HTML QC report**

> This is meant to complement tools like MetaboAnalyst (not replace them). The goal here is fast, repeatable QC and shareable outputs inside a facility/lab workflow.

---

## Quick start (local)
From the repo root:

```r
source("run_local.R")
```

Or:

```r
shiny::runApp()
```

## Input expectations
### Feature table
- Rows = features
- Columns include:
  - descriptor fields (optional, e.g. `row m/z`, `row retention time`, etc.)
  - **sample intensity columns** (peak areas)

If your sample columns end with ` Peak area`, the app can auto-detect them. Otherwise, use **Manual select**.

### Metadata (optional)
- Must contain a join column with sample identifiers (default recommended: `filename`)
- Other columns can be anything you want to group/color by (e.g. `batch`, `run`, `plate`, `operator`, `ATTRIBUTE_species`, etc.)

---

## Repo structure
```
metabolomics-qc-portal-shiny/
├─ app.R
├─ run_local.R
├─ report_template.Rmd
├─ DEPLOYMENT.md
├─ R/
│  ├─ io.R        # file reading, sample column detection, validation helpers
│  ├─ qc.R        # QC computations + group-level (batch/run) QC summary
│  ├─ eda.R       # PCA + loadings logic
│  └─ report.R    # exports + HTML report render helper
├─ data/
│  └─ (demo files)
└─ docs/
   └─ screenshots/
```

---

## Screenshots
Add two quick screenshots so a reviewer understands the app in ~10 seconds.

Place these files in `docs/screenshots/`:
- `docs/screenshots/pca.png` — PCA tab (colored by a metadata column)
- `docs/screenshots/exports.png` — Exports tab (showing the report + CSV download buttons)

![PCA view](docs/screenshots/pca.png)
![Exports view](docs/screenshots/exports.png)

---

## Notes for hosting
See `DEPLOYMENT.md` for internal hosting options (Posit Connect / Shiny Server / container).

## License
For a portfolio repo, an MIT license is usually fine. (Add one if you plan to share publicly.)
