# EDA / PCA utilities for Metabolomics QC Portal
# Imported by app.R via: source("R/eda.R")

median_impute <- function(x) {
  if (all(is.na(x))) return(rep(0, length(x)))
  x[is.na(x)] <- stats::median(x, na.rm = TRUE)
  x
}

prep_for_pca <- function(int_mat_df, feature_labels) {

  m <- as.matrix(int_mat_df)

  # Impute missing values per sample (column-wise) using the sample median.
  for (j in seq_len(ncol(m))) {
    m[, j] <- median_impute(m[, j])
  }

  # Drop features with no signal anywhere.
  keep_feature <- rowSums(m != 0) > 0
  m <- m[keep_feature, , drop = FALSE]
  feature_labels <- feature_labels[keep_feature]

  # PCA expects samples x features.
  x <- t(m)

  # Reduce skew from very large peaks.
  x <- log1p(x)

  # Remove zero-variance features.
  sds <- apply(x, 2, stats::sd)
  keep_var <- sds > 0
  x <- x[, keep_var, drop = FALSE]
  feature_labels <- feature_labels[keep_var]

  list(x = x, feature_labels = feature_labels)
}

run_pca <- function(sample_by_feature_matrix) {
  pca <- stats::prcomp(sample_by_feature_matrix, center = TRUE, scale. = TRUE)
  scores <- as.data.frame(pca$x)
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
  list(model = pca, scores = scores, variance = var_expl)
}

pc_index <- function(pc_name) {
  as.integer(stringr::str_remove(pc_name, "^PC"))
}
