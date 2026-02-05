# Deployment notes (internal / facility use)

This repo is a Shiny app meant to be run locally for development and demo, and hosted internally when ready.

## Quick local run (recommended for reviewers)
From the project root:

```r
source("run_local.R")
```

Or:

```r
shiny::runApp()
```

## What the app expects
- A **feature table** (CSV/TSV) where sample intensity columns are either:
  - auto-detected as TIMA-style columns ending with ` Peak area`, or
  - manually selected by the user.
- Optional **metadata** (CSV/TSV) with a join column (default: `filename`) to match sample IDs.

## Hosting options (common in academic medical centers)

### 1) Posit Connect (best “managed” option)
- Upload the repo as a Shiny application bundle.
- Configure a system R + package library.
- Set an app access group (e.g., metabolomics facility users).

### 2) Shiny Server (open source)
- Place the repo folder under the server’s app directory (often `/srv/shiny-server/`).
- Ensure the server user can read the folder and installed packages.

### 3) Container (future-proof)
If your environment prefers containers, package this as a Docker image:
- Base on `rocker/r-ver` (or a standard internal base image)
- Install system deps + R packages
- Expose 3838

(We can add a `Dockerfile` later if you want.)

## Basic security / hygiene
- This app does **not** store uploaded files; it keeps them in memory for the current session.
- Avoid uploading PHI / identifiers unless your hosting environment is approved for it.
- If you need audit logging (who uploaded what), host behind your institution’s auth (SSO/VPN) and add server-side logs.

## Performance notes
- The app is designed for “feature tables” that fit in memory.
- If you expect very large uploads (hundreds of MB+), consider:
  - tighter file size limits
  - streaming reads / chunking
  - moving PCA to a background job queue

## Suggested next “production” upgrades (optional)
- Support MS vendor exports beyond TIMA naming conventions (more robust sample column detection).
- Add “injection order” trend plots if metadata contains `injection_order`.
- Add a small “admin” config for facility defaults (preferred join key, preferred group columns).
