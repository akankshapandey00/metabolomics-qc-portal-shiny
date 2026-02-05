# QC utilities for Metabolomics QC Portal
# Imported by app.R via: source("R/qc.R")

qc_by_sample <- function(int_mat, sample_cols_display,
                         missing_flag_pct = 20,
                         zeros_flag_pct = 80) {

  # Zeros are common in metabolomics (not detected).
  # Missing values (NA) usually indicate parsing/alignment issues.

  missing_n <- sapply(int_mat, function(x) sum(is.na(x)))
  missing_pct <- sapply(int_mat, function(x) mean(is.na(x)) * 100)
  non_missing_n <- sapply(int_mat, function(x) sum(!is.na(x)))

  zero_n <- sapply(int_mat, function(x) sum(!is.na(x) & x == 0))
  zero_pct <- sapply(int_mat, function(x) {
    denom <- sum(!is.na(x))
    if (denom == 0) return(NA_real_)
    (sum(x[!is.na(x)] == 0) / denom) * 100
  })

  median_all <- sapply(int_mat, function(x) median(x, na.rm = TRUE))
  median_nonzero <- sapply(int_mat, function(x) {
    y <- x[!is.na(x) & x != 0]
    if (length(y) == 0) return(NA_real_)
    median(y)
  })

  tibble::tibble(
    sample_id = sample_cols_display,
    non_missing_n = non_missing_n,
    missing_n = missing_n,
    missing_pct = round(missing_pct, 2),
    zero_n = zero_n,
    zero_pct = round(zero_pct, 2),
    median_intensity = round(median_all, 3),
    median_nonzero_intensity = round(median_nonzero, 3),
    flag_high_missing = missing_pct > missing_flag_pct,
    flag_many_zeros = zero_pct > zeros_flag_pct
  ) %>%
    dplyr::arrange(dplyr::desc(flag_high_missing), dplyr::desc(flag_many_zeros),
                  dplyr::desc(missing_pct), dplyr::desc(zero_pct))
}

qc_by_feature <- function(int_mat, feature_labels = NULL) {

  tbl <- tibble::tibble(
    feature_idx = seq_len(nrow(int_mat)),
    missing_n = apply(int_mat, 1, function(x) sum(is.na(x))),
    missing_pct = apply(int_mat, 1, function(x) mean(is.na(x)) * 100),
    median_intensity = apply(int_mat, 1, function(x) median(x, na.rm = TRUE)),
    nonzero_n = apply(int_mat, 1, function(x) sum(!is.na(x) & x != 0)),
    nonzero_pct = apply(int_mat, 1, function(x) {
      denom <- sum(!is.na(x))
      if (denom == 0) return(NA_real_)
      (sum(x[!is.na(x)] != 0) / denom) * 100
    })
  ) %>%
    dplyr::mutate(
      missing_pct = round(missing_pct, 2),
      median_intensity = round(median_intensity, 3),
      nonzero_pct = round(nonzero_pct, 2)
    )

  if (!is.null(feature_labels) && length(feature_labels) == nrow(tbl)) {
    tbl <- tbl %>%
      dplyr::mutate(feature = feature_labels) %>%
      dplyr::select(feature, dplyr::everything())
  }

  tbl %>% dplyr::arrange(dplyr::desc(missing_pct), dplyr::desc(nonzero_pct))
}

qc_by_group <- function(qc_sample_df, group_col) {
  # Quick facility-style summary: "are batches/runs behaving differently?"
  # Expects qc_sample_df to already have the metadata column present.

  if (!group_col %in% names(qc_sample_df)) {
    stop("Grouping column not found in qc_sample_df: ", group_col)
  }

  qc_sample_df %>%
    dplyr::mutate(.group = as.character(.data[[group_col]])) %>%
    dplyr::filter(!is.na(.group) & .group != "") %>%
    dplyr::group_by(.group) %>%
    dplyr::summarise(
      n_samples = dplyr::n(),
      median_missing_pct = round(stats::median(missing_pct, na.rm = TRUE), 2),
      median_zero_pct = round(stats::median(zero_pct, na.rm = TRUE), 2),
      median_signal = round(stats::median(median_intensity, na.rm = TRUE), 3),
      median_nonzero_signal = round(stats::median(median_nonzero_intensity, na.rm = TRUE), 3),
      flagged_high_missing = sum(flag_high_missing %in% TRUE, na.rm = TRUE),
      flagged_many_zeros = sum(flag_many_zeros %in% TRUE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(median_missing_pct), dplyr::desc(median_zero_pct))
}

