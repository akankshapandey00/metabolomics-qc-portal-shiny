# I/O + validation helpers for Metabolomics QC Portal
# Imported by app.R via: source("R/io.R")

read_table_any <- function(path, ext) {
  # Keep read behavior predictable for end users:
  # - TSV/TXT -> read_tsv
  # - CSV or anything else -> read_csv
  if (ext %in% c("tsv", "txt")) readr::read_tsv(path, show_col_types = FALSE)
  else readr::read_csv(path, show_col_types = FALSE)
}

detect_sample_cols <- function(df) {
  # TIMA-style: sample intensity columns often end with " Peak area"
  names(df)[stringr::str_detect(names(df), " Peak area$")]
}

strip_peak_area_suffix <- function(x) stringr::str_remove(x, " Peak area$")

normalize_sample_ids <- function(sample_cols) {
  # For join checks, we want sample IDs without the " Peak area" suffix.
  out <- sample_cols
  if (any(stringr::str_detect(out, " Peak area$"))) out <- strip_peak_area_suffix(out)
  out
}

get_intensity_matrix <- function(features_df, sample_cols) {
  # Pull sample columns and coerce to numeric.
  # We don't throw here — we validate separately so we can show a clean message in the UI.
  features_df %>%
    dplyr::select(dplyr::all_of(sample_cols)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ suppressWarnings(as.numeric(.x))))
}

build_feature_labels <- function(features_df) {
  # Human-readable feature labels for the loadings table.
  has_id <- "row ID" %in% names(features_df)
  has_mz <- "row m/z" %in% names(features_df)
  has_rt <- "row retention time" %in% names(features_df)

  if (has_id || has_mz || has_rt) {
    id_part <- if (has_id) paste0("ID=", features_df[["row ID"]]) else paste0("idx=", seq_len(nrow(features_df)))
    mz_part <- if (has_mz) paste0(" mz=", round(features_df[["row m/z"]], 4)) else ""
    rt_part <- if (has_rt) paste0(" rt=", round(features_df[["row retention time"]], 3)) else ""
    return(paste0(id_part, mz_part, rt_part))
  }
  paste0("feature_idx=", seq_len(nrow(features_df)))
}

validate_intensity_matrix <- function(int_mat, sample_cols) {
  # Validate that selected sample columns are usable numeric intensity columns.
  # Returns a list with:
  # - ok: TRUE/FALSE
  # - message: user-facing string
  # - bad_cols: character vector of problematic columns
  bad <- character(0)

  for (i in seq_along(sample_cols)) {
    col_name <- sample_cols[i]
    x <- int_mat[[i]]

    # If everything became NA after coercion, it's almost certainly not an intensity column.
    if (all(is.na(x))) {
      bad <- c(bad, col_name)
      next
    }

    # If most values are NA, it's probably the wrong column (e.g. text)
    na_frac <- mean(is.na(x))
    if (is.finite(na_frac) && na_frac > 0.90) {
      bad <- c(bad, col_name)
    }
  }

  if (length(bad) > 0) {
    msg <- paste0(
      "Some selected sample columns don't look numeric (or became mostly NA after conversion). ",
      "Please re-check your selection. Problem columns: ",
      paste(utils::head(bad, 10), collapse = ", "),
      if (length(bad) > 10) " …" else ""
    )
    return(list(ok = FALSE, message = msg, bad_cols = bad))
  }

  list(ok = TRUE, message = NULL, bad_cols = character(0))
}
