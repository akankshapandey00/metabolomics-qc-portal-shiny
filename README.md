# Metabolomics QC Portal (Shiny)

A lightweight, facility-style web app for **quick QC + PCA** on metabolomics **feature tables**.  
Built to help a metabolomics core (or any lab) sanity-check data fast, spot technical effects (batch/run), and export clean tables + a shareable HTML report.

This project complements tools like **MetaboAnalyst** (it is *not* meant to replace full downstream analysis pipelines). The value here is: **fast upload → validation → QC triage → shareable outputs**.

---

## Table of contents
- [What this app does](#what-this-app-does)
- [Typical workflow](#typical-workflow)
- [Inputs](#inputs)
  - [Feature table](#feature-table)
  - [Metadata table](#metadata-table)
  - [Common format issues](#common-format-issues)
- [Tabs and outputs](#tabs-and-outputs)
  - [Features / Metadata](#features--metadata)
  - [Join check](#join-check)
  - [QC](#qc)
  - [PCA](#pca)
  - [Batch QC](#batch-qc)
  - [Exports](#exports)
- [QC metrics: definitions](#qc-metrics-definitions)
- [PCA: what it’s doing](#pca-what-its-doing)
- [HTML report](#html-report)
- [Install and run](#install-and-run)
- [Repo structure](#repo-structure)
- [Development notes](#development-notes)


---

## What this app does
**Core features**
- Upload a **feature table** (CSV/TSV) and optional **metadata** (CSV/TSV)
- **Auto-detect** TIMA-style intensity columns ending with ` Peak area` (or manually select them)
- Validate that chosen sample columns are numeric intensities (helps prevent “wrong column selected” issues)
- **Join check** to confirm sample IDs in the feature table match metadata IDs
- Compute:
  - **QC by sample**
  - **QC by feature**
- Run **PCA** with:
  - variance explained (%)
  - top loadings table (features driving a selected PC)
  - optional coloring by any metadata column
  - optional custom axes (e.g., PC1 vs PC3)
- **Batch QC**:
  - group summary table (median missing%, median zero%, median signal per group)
  - signal distribution plot by group
- Export:
  - QC tables as CSV
  - PCA scores/loadings as CSV
  - filtered feature table as CSV (simple facility-friendly filter)
  - **HTML QC report** for sharing

---

## Typical workflow
1. **Load demo dataset** (optional) to confirm everything works end-to-end.
2. Upload your **feature table** and (optional) **metadata**.
3. Confirm sample intensity columns:
   - Auto mode (if columns end with ` Peak area`)
   - or Manual select (choose your intensity columns)
4. Use **Join check** to verify sample IDs match metadata (fix join key if needed).
5. Review **QC** + **PCA**:
   - look for outliers, high missingness, unusual signal distributions
   - check PCA separation by run/batch/plate/operator if available in metadata
6. Use **Batch QC** to quickly confirm whether technical groups differ.
7. Export CSV outputs and/or the **HTML report**.

---

## Inputs

### Feature table
A “feature table” is expected in **wide format**:

- Rows = features (m/z features, annotated peaks, etc.)
- Columns = feature descriptors (optional) + sample intensity columns

**Common descriptor columns (optional):**
- `row ID`
- `row m/z`
- `row retention time`
- anything else that describes the feature

**Intensity columns:**
- numeric intensities for each sample (often peak areas)

**Auto-detection**
If your intensity columns end with ` Peak area`, the app will detect them in Auto mode.
Example:
- `191105_AR_Panax_Pos.mzML Peak area`
- `191105_AR_Ginkgo_Pos.mzML Peak area`

If your file uses different naming, switch to **Manual select** and pick the correct columns yourself.

### Metadata table
Metadata is optional but strongly recommended if you want to:
- color PCA by group/batch/run
- do Batch QC summaries

Requirements:
- A join column that matches your sample IDs (default recommended: `filename`)
- Any other columns are optional and can be used for grouping/coloring (examples):
  - `batch`, `run`, `plate`, `operator`, `instrument`
  - biological groups (`case_control`, `treatment`, `timepoint`)
  - any label you use in your workflow

### Common format issues
- **Mismatch between feature table sample names and metadata IDs**  
  Fix by changing the join key, or normalize naming upstream.
- **Non-numeric intensity columns**  
  If you accidentally select a text column as an “intensity” column, the app will warn you.
- **Comma vs tab** hookup issues  
  Use the correct file extension (`.csv`, `.tsv`, `.txt`) so parsing is consistent.

---

## Tabs and outputs

### Features / Metadata
Quick previews of the uploaded data. Useful for:
- verifying the file loaded correctly
- checking column names
- confirming your join column exists

### Join check
Shows whether detected/selected sample IDs are present in metadata.

- `in_metadata = TRUE` means the sample ID was found in your metadata join column.
- If most rows are `FALSE`, your join key is likely wrong (or sample naming differs).

### QC
Two tables:

**QC by sample**
- one row per sample
- highlights samples with high missingness or extreme zeros

**QC by feature**
- one row per feature
- helps spot features that are mostly missing or mostly zero

### PCA
- PCA scatter plot
- variance explained table (PC1…PC10)
- top loadings table (features driving a selected PC)

You can:
- color by any metadata column
- toggle sample labels
- choose x/y axes (PC1 vs PC2 is default)

### Batch QC
Designed for facility-style “quick triage”.

- Select a metadata column to group by (e.g., `batch`, `run`, `plate`, `operator`)
- View:
  - group summary QC table
  - signal distribution plot (sample median intensity) per group

### Exports
Download:
- QC by sample (CSV)
- QC by feature (CSV)
- PCA scores (CSV)
- top PCA loadings (CSV)
- filtered feature table (CSV)
- HTML QC report

---

## QC metrics: definitions
These are meant to be intuitive “triage” metrics (not a full statistical pipeline).

**Missing %**
- fraction of values that are NA
- high missing% can indicate parsing, alignment, or detection issues

**Zero %**
- fraction of **zeros** among **non-missing** values
- zeros are common in metabolomics (not detected), but extremely high zero% can indicate low signal or poor detection

**Median intensity**
- per sample: median across features
- a quick signal-level sanity check

**Median non-zero intensity**
- per sample: median of non-zero values
- useful when many features are zero for a sample

Flags (simple heuristics)
- `flag_high_missing`: sample missing% above a basic threshold
- `flag_many_zeros`: sample zero% above a basic threshold

(Thresholds are intentionally simple. Facilities often adjust them to their own assays.)

---

## PCA: what it’s doing
PCA is run on the intensity matrix after a lightweight preprocessing step:

- Select intensity columns
- Convert to numeric (with validation)
- Transform: `log1p(intensity)` (stabilizes scale and reduces heavy tails)
- Drop features that are:
  - all zeros
  - zero variance
- Run PCA with centering + scaling
- Plot requested PCs (default PC1 vs PC2)
- Variance explained is computed from PCA standard deviations

Interpretation tip:
- PCA separates samples by the **largest sources of variation**.
- If separation aligns with **batch/run/plate/operator**, it often suggests technical effects.
- If separation aligns with **biology/treatment**, that’s often expected (but still confirm QC).

---

## HTML report
The app renders an HTML report from `report_template.Rmd`.  
The report includes:
- QC by sample
- QC by feature
- PCA plot and variance
- top loadings table
- session info (package versions)

This is meant to be something you can attach in email/slack or drop into a project folder.

---

## Install and run

### 1) Quick run (recommended)
From the repo root:

```r
source("run_local.R")
```

This script checks for packages and installs missing ones once.

### 2) Standard run
```r
shiny::runApp()
```

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

```

---

## Development notes
Why modularize?
- Keeping `app.R` mostly UI + wiring makes it easier to:
  - test functions
  - reuse computation logic
  - maintain the app when multiple people touch it

Where to add new features:
- parsing & validation → `R/io.R`
- QC metrics → `R/qc.R`
- PCA / exploration → `R/eda.R`
- exports / reports → `R/report.R`

---

