# Export + report helpers for Metabolomics QC Portal
# Imported by app.R via: source("R/report.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_stub <- function(x) {
  x <- tolower(x %||% "dataset")
  x <- stringr::str_replace_all(x, "\\.[A-Za-z0-9]+$", "")
  x <- stringr::str_replace_all(x, "[^a-z0-9]+", "_")
  x <- stringr::str_replace_all(x, "_+", "_")
  x <- stringr::str_replace_all(x, "^_|_$", "")
  if (nchar(x) == 0) x <- "dataset"
  x
}

timestamp_tag <- function() format(Sys.Date(), "%Y-%m-%d")

filter_features <- function(features_df, sample_cols, missing_max = 20, nonzero_min = 20) {
  # Apply a simple, facility-friendly filter:
  # - keep features with <= missing_max missing %
  # - keep features with >= nonzero_min non-zero % among non-missing values
  int_mat <- get_intensity_matrix(features_df, sample_cols)

  miss_pct <- apply(int_mat, 1, function(x) mean(is.na(x)) * 100)
  nonzero_pct <- apply(int_mat, 1, function(x) {
    denom <- sum(!is.na(x))
    if (denom == 0) return(0)
    (sum(x[!is.na(x)] != 0) / denom) * 100
  })

  keep <- (miss_pct <= missing_max) & (nonzero_pct >= nonzero_min)

  descriptor_cols <- setdiff(names(features_df), sample_cols)
  features_df[keep, c(descriptor_cols, sample_cols), drop = FALSE]
}

render_qc_report <- function(template_path,
                             dataset_name,
                             qc_sample_df,
                             qc_feature_df,
                             pca_scores_df,
                             pca_top_loadings_df,
                             pca_color_col = "",
                             pca_x_pc = "PC1",
                             pca_y_pc = "PC2",
                             top_n_loadings = 10) {

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required for report export. Install it with install.packages('rmarkdown').")
  }
  if (!file.exists(template_path)) {
    stop("Missing report_template.Rmd. Put it in the same folder as app.R.")
  }

  tmp <- tempfile("qc_report_")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)

  qc_sample_path <- file.path(tmp, "qc_by_sample.csv")
  qc_feature_path <- file.path(tmp, "qc_by_feature.csv")
  pca_scores_path <- file.path(tmp, "pca_scores.csv")
  pca_loadings_path <- file.path(tmp, "pca_top_loadings.csv")

  readr::write_csv(qc_sample_df, qc_sample_path)
  readr::write_csv(qc_feature_df, qc_feature_path)
  readr::write_csv(pca_scores_df, pca_scores_path)
  readr::write_csv(pca_top_loadings_df, pca_loadings_path)

  out_html <- file.path(tmp, "report.html")

  rmarkdown::render(
    input = template_path,
    output_file = out_html,
    params = list(
      report_date = Sys.Date(),
      dataset_name = dataset_name,
      qc_sample_csv = qc_sample_path,
      qc_feature_csv = qc_feature_path,
      pca_scores_csv = pca_scores_path,
      pca_loadings_csv = pca_loadings_path,
      pca_color_col = pca_color_col,
      pca_x_pc = pca_x_pc,
      pca_y_pc = pca_y_pc,
      top_n_loadings = top_n_loadings
    ),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )

  list(html_path = out_html, workdir = tmp)
}
