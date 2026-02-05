# Metabolomics QC Portal (Facility-ready MVP)
# Modularized: QC + PCA logic lives in R/qc.R and R/eda.R

options(shiny.maxRequestSize = 300 * 1024^2)

library(shiny)
library(DT)
library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(ggplot2)

# ---- Source modular code ----
# Keep these files in: metabolomics-qc-portal-shiny/R/
source("R/qc.R")
source("R/eda.R")

# ---------- small helpers ----------

`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_stub <- function(x) {
  x <- tolower(x %||% "dataset")
  x <- str_replace_all(x, "\\.[A-Za-z0-9]+$", "")
  x <- str_replace_all(x, "[^a-z0-9]+", "_")
  x <- str_replace_all(x, "_+", "_")
  x <- str_replace_all(x, "^_|_$", "")
  if (nchar(x) == 0) x <- "dataset"
  x
}

# ---------- I/O + column helpers ----------

read_table_any <- function(path, ext) {
  if (ext %in% c("tsv", "txt")) readr::read_tsv(path, show_col_types = FALSE)
  else readr::read_csv(path, show_col_types = FALSE)
}

detect_sample_cols <- function(df) {
  names(df)[str_detect(names(df), " Peak area$")]
}

strip_peak_area_suffix <- function(x) str_remove(x, " Peak area$")

normalize_sample_ids <- function(sample_cols) {
  out <- sample_cols
  if (any(str_detect(out, " Peak area$"))) out <- strip_peak_area_suffix(out)
  out
}

get_intensity_matrix <- function(features_df, sample_cols) {
  features_df %>%
    select(all_of(sample_cols)) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))
}

build_feature_labels <- function(features_df) {
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

# ---------- UI ----------

ui <- fluidPage(
  titlePanel("Metabolomics QC Portal (Facility-ready MVP)"),

  sidebarLayout(
    sidebarPanel(
      h4("Upload"),
      fileInput("features_file", "Feature table (CSV/TSV)", accept = c(".csv", ".tsv", ".txt")),
      fileInput("metadata_file", "Metadata (CSV/TSV)", accept = c(".csv", ".tsv", ".txt")),

      hr(),
      h4("Column mapping"),
      uiOutput("join_key_ui"),

      radioButtons(
        "sample_mode",
        "Sample columns",
        choices = c(
          "Auto (TIMA-style: ends with ' Peak area')" = "tima",
          "Manual select" = "manual"
        ),
        selected = "tima"
      ),
      uiOutput("manual_sample_cols_ui"),

      hr(),
      h4("PCA settings"),
      uiOutput("pca_color_ui"),
      checkboxInput("pca_show_labels", "Show sample labels", value = FALSE),
      uiOutput("pca_axis_ui"),

      hr(),
      h4("Loadings"),
      selectInput("loading_pc", "Show loadings for", choices = c("PC1", "PC2"), selected = "PC1"),
      numericInput("loading_top_n", "Top N features", value = 10, min = 5, max = 200, step = 5),

      hr(),
      h4("Exports"),
      numericInput("export_feature_missing_max", "Max missing % per feature", value = 20, min = 0, max = 100, step = 5),
      numericInput("export_feature_nonzero_min", "Min non-zero % per feature", value = 20, min = 0, max = 100, step = 5),

      hr(),
      h4("Demo"),
      actionButton("load_demo", "Load demo dataset (for testing)"),

      hr(),
      h4("Status"),
      verbatimTextOutput("status")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Home",
          h3("Metabolomics QC Portal"),
          tags$p("A lightweight QC + PCA web tool for metabolomics feature tables. Built for fast facility triage: check data quality, spot batch/run effects, and export a shareable report."),

          wellPanel(
            strong("Quick start"),
            tags$ol(
              tags$li("Click “Load demo dataset” to test the app end-to-end."),
              tags$li("Upload your feature table and (optional) metadata."),
              tags$li("Confirm sample columns (Auto or Manual select)."),
              tags$li("Use “Join check” to verify metadata matches sample IDs."),
              tags$li("Review QC, PCA, and Batch QC, then export CSVs / HTML report.")
            )
          ),

          hr(),
          h4("For reviewers / collaborators"),
          tags$ul(
            tags$li(tags$b("One-click run:"), " in the project root, run ", tags$code('source("run_local.R")')),
            tags$li(tags$b("Hosting notes:"), " see ", tags$code("DEPLOYMENT.md"), " for internal deployment options")
          ),

          hr(),
          h4("What this tool is (and isn’t)"),
          tags$ul(
            tags$li(tags$b("Is:"), " quick QC + PCA + batch/run checks and clean exports."),
            tags$li(tags$b("Isn’t:"), " a full metabolomics pipeline replacement (MetaboAnalyst and other tools still have broader downstream stats/biology workflows).")
          )
        ),

        tabPanel("Features", DTOutput("features_preview")),
        tabPanel("Metadata", DTOutput("metadata_preview")),

        tabPanel("Join check",
          uiOutput("join_help_ui"),
          DTOutput("join_check")
        ),

        tabPanel("QC",
          uiOutput("qc_help_ui"),
          h4("QC by sample"),
          DTOutput("qc_sample_table"),
          hr(),
          h4("QC by feature (first 200)"),
          DTOutput("qc_feature_table")
        ),

        tabPanel("PCA",
          uiOutput("pca_help_ui"),
          plotOutput("pca_plot", height = 420),
          hr(),
          h4("Variance explained"),
          DTOutput("pca_variance_table"),
          hr(),
          h4("Top features driving PCs"),
          DTOutput("pca_loadings_table")
        ),

        tabPanel("Exports",
          h4("Export summary"),
          wellPanel(verbatimTextOutput("export_summary")),
          hr(),
          h4("Download outputs"),
          downloadButton("download_qc_sample", "Download QC by sample (CSV)"),
          tags$br(), tags$br(),
          downloadButton("download_qc_feature", "Download QC by feature (CSV)"),
          tags$br(), tags$br(),
          downloadButton("download_pca_scores", "Download PCA scores (CSV)"),
          tags$br(), tags$br(),
          downloadButton("download_pca_loadings", "Download top PCA loadings (CSV)"),
          tags$br(), tags$br(),
          downloadButton("download_filtered_features", "Download filtered feature table (CSV)"),
          tags$br(), tags$br(),
          downloadButton("download_html_report", "Download HTML QC report")
        ),

        tabPanel("About / Methods",
          h3("What this portal does"),
          tags$p("This is a lightweight QC + PCA portal for metabolomics feature tables. It’s designed to help a facility (or a lab) quickly sanity-check data, spot technical issues, and generate shareable outputs."),
          hr(),
          h4("QC metrics"),
          tags$ul(
            tags$li(tags$b("Missing %:"), " fraction of NA values in a sample/feature. This usually indicates data parsing or alignment issues."),
            tags$li(tags$b("Zero %:"), " fraction of zeros among non-missing values. Zeros are common in metabolomics (not detected), but extremely high zero % can indicate low signal."),
            tags$li(tags$b("Median intensity:"), " quick indicator of global signal level across samples.")
          ),
          hr(),
          h4("PCA"),
          tags$ul(
            tags$li("PCA reduces high-dimensional data to a small number of components (PCs) that explain the largest variance."),
            tags$li("We use log1p intensities, drop all-zero and zero-variance features, then run PCA with centering + scaling."),
            tags$li("The loadings table highlights which features contribute most to a selected PC.")
          ),
          hr(),
          h4("Reproducibility"),
          tags$p("The HTML report includes session information so results can be traced back to package versions.")
        )
      )
    )
  )
)

# ---------- Server ----------

server <- function(input, output, session) {

  rv <- reactiveValues(features = NULL, metadata = NULL, sample_cols = NULL)

  dataset_stub <- reactive({
    nm <- if (!is.null(input$features_file$name)) input$features_file$name else "demo_dataset"
    safe_stub(nm)
  })

  timestamp_tag <- reactive({
    format(Sys.Date(), "%Y-%m-%d")
  })

  observeEvent(input$load_demo, {
    feat <- read_table_any("data/tima_example_features.csv", "csv")
    meta <- read_table_any("data/tima_example_metadata.tsv", "tsv")

    rv$features <- feat
    rv$metadata <- meta
    rv$sample_cols <- detect_sample_cols(feat)
  })

  observeEvent(input$features_file, {
    ext <- tools::file_ext(input$features_file$name)
    feat <- read_table_any(input$features_file$datapath, ext)

    rv$features <- feat
    rv$sample_cols <- detect_sample_cols(feat)
  })

  observeEvent(input$metadata_file, {
    ext <- tools::file_ext(input$metadata_file$name)
    rv$metadata <- read_table_any(input$metadata_file$datapath, ext)
  })

  output$join_key_ui <- renderUI({
    req(rv$metadata)
    cols <- names(rv$metadata)
    default <- if ("filename" %in% cols) "filename" else cols[1]
    selectInput("join_key", "Metadata join column", choices = cols, selected = default)
  })

  output$manual_sample_cols_ui <- renderUI({
    req(rv$features)
    if (input$sample_mode != "manual") return(NULL)

    selectInput(
      "manual_sample_cols",
      "Select sample columns (intensity columns)",
      choices = names(rv$features),
      selected = rv$sample_cols,
      multiple = TRUE
    )
  })

  sample_cols_selected <- reactive({
    req(rv$features)
    if (input$sample_mode == "manual") {
      req(input$manual_sample_cols)
      return(input$manual_sample_cols)
    }
    detect_sample_cols(rv$features)
  })

  output$pca_color_ui <- renderUI({
    if (is.null(rv$metadata)) return(helpText("Upload metadata to color PCA by group/batch/etc."))

    join_key <- input$join_key
    cols <- names(rv$metadata)
    if (!is.null(join_key) && join_key %in% cols) cols <- setdiff(cols, join_key)
    if (length(cols) == 0) return(helpText("No metadata columns available for coloring PCA."))

    selectInput("pca_color_col", "Color PCA by", choices = cols, selected = cols[1])
  })

  output$status <- renderText({
    if (is.null(rv$features) && is.null(rv$metadata)) return("No data loaded yet.")

    features_msg <- if (!is.null(rv$features)) {
      paste0("Features: ", nrow(rv$features), " rows × ", ncol(rv$features), " columns")
    } else "Features: not loaded"

    sample_msg <- if (!is.null(rv$features)) {
      paste0("Sample columns (current): ", length(sample_cols_selected()))
    } else "Sample columns (current): n/a"

    meta_msg <- if (!is.null(rv$metadata)) {
      paste0("Metadata: ", nrow(rv$metadata), " rows; columns: ", paste(names(rv$metadata), collapse = ", "))
    } else "Metadata: not loaded"

    paste(features_msg, sample_msg, meta_msg, sep = "\n")
  })

  output$join_help_ui <- renderUI({
    if (is.null(rv$features)) return(helpText("Upload a feature table to start."))
    if (is.null(rv$metadata)) return(helpText("Upload metadata to run the join check (or use demo data)."))
    NULL
  })

  output$qc_help_ui <- renderUI({
    if (is.null(rv$features)) return(helpText("Upload a feature table (or load the demo) to see QC tables."))
    if (length(sample_cols_selected()) == 0) return(helpText("No sample columns detected/selected yet."))
    NULL
  })

  output$pca_help_ui <- renderUI({
    if (is.null(rv$features)) return(helpText("Upload a feature table (or load the demo) to run PCA."))
    if (length(sample_cols_selected()) < 2) return(helpText("PCA needs at least 2 sample columns."))
    NULL
  })

  output$features_preview <- renderDT({
    req(rv$features)
    datatable(head(rv$features, 20), options = list(scrollX = TRUE, pageLength = 10))
  })

  output$metadata_preview <- renderDT({
    req(rv$metadata)
    datatable(head(rv$metadata, 200), options = list(scrollX = TRUE, pageLength = 10))
  })

  output$join_check <- renderDT({
    req(rv$features, rv$metadata)
    req(input$join_key)

    sample_cols <- sample_cols_selected()
    if (length(sample_cols) == 0) {
      return(datatable(data.frame(note = "No sample columns selected/detected yet."), options = list(dom = "t")))
    }

    meta_ids <- rv$metadata[[input$join_key]]
    if (is.null(meta_ids)) {
      return(datatable(data.frame(note = "Selected join key not found in metadata."), options = list(dom = "t")))
    }

    sample_ids <- normalize_sample_ids(sample_cols)
    out <- tibble(sample_id = sample_ids, in_metadata = sample_ids %in% meta_ids) %>%
      arrange(in_metadata)

    datatable(out, options = list(pageLength = 15, scrollX = TRUE))
  })

  qc_sample_tbl <- reactive({
    req(rv$features)
    sample_cols <- sample_cols_selected()
    req(length(sample_cols) > 0)

    int_mat <- get_intensity_matrix(rv$features, sample_cols)
    sample_ids <- normalize_sample_ids(sample_cols)
    qc_by_sample(int_mat, sample_ids)
  })

  qc_feature_tbl <- reactive({
    req(rv$features)
    sample_cols <- sample_cols_selected()
    req(length(sample_cols) > 0)

    int_mat <- get_intensity_matrix(rv$features, sample_cols)
    feat_labels <- build_feature_labels(rv$features)
    qc_by_feature(int_mat, feat_labels)
  })

  output$qc_sample_table <- renderDT({
    if (is.null(rv$features) || length(sample_cols_selected()) == 0) {
      return(datatable(data.frame(note = "Upload/select sample columns to view QC."), options = list(dom = "t")))
    }
    datatable(qc_sample_tbl(), options = list(pageLength = 10, scrollX = TRUE))
  })

  output$qc_feature_table <- renderDT({
    if (is.null(rv$features) || length(sample_cols_selected()) == 0) {
      return(datatable(data.frame(note = "Upload/select sample columns to view QC."), options = list(dom = "t")))
    }
    datatable(head(qc_feature_tbl(), 200), options = list(pageLength = 10, scrollX = TRUE))
  })

  pca_result <- reactive({
    req(rv$features)
    sample_cols <- sample_cols_selected()
    req(length(sample_cols) >= 2)

    int_mat <- get_intensity_matrix(rv$features, sample_cols)
    feature_labels <- build_feature_labels(rv$features)

    prep <- prep_for_pca(int_mat, feature_labels)
    res <- run_pca(prep$x)

    sample_ids <- normalize_sample_ids(sample_cols)
    scores <- res$scores
    scores$sample_id <- sample_ids

    if (!is.null(rv$metadata) &&
        !is.null(input$join_key) &&
        input$join_key %in% names(rv$metadata)) {
      scores <- scores %>%
        left_join(rv$metadata, by = setNames(input$join_key, "sample_id"))
    }

    rotation <- as.data.frame(res$model$rotation)
    rotation$feature <- prep$feature_labels

    list(scores = scores, variance = res$variance, loadings = rotation)
  })

  output$pca_axis_ui <- renderUI({
    res <- try(pca_result(), silent = TRUE)
    if (inherits(res, "try-error")) return(helpText("Axes: PC1 vs PC2"))

    pcs <- names(res$scores)
    pcs <- pcs[str_detect(pcs, "^PC\\d+$")]
    if (length(pcs) == 0) return(helpText("Axes: PC1 vs PC2"))

    tagList(
      selectInput("pca_x_pc", "X-axis", choices = pcs, selected = "PC1"),
      selectInput("pca_y_pc", "Y-axis", choices = pcs, selected = "PC2")
    )
  })

  output$pca_plot <- renderPlot({
    req(rv$features)
    req(length(sample_cols_selected()) >= 2)

    res <- pca_result()
    scores <- res$scores
    v <- res$variance

    pcs <- names(scores)
    pcs <- pcs[str_detect(pcs, "^PC\\d+$")]

    x_pc <- if (!is.null(input$pca_x_pc) && input$pca_x_pc %in% pcs) input$pca_x_pc else "PC1"
    y_pc <- if (!is.null(input$pca_y_pc) && input$pca_y_pc %in% pcs) input$pca_y_pc else "PC2"
    if (x_pc == y_pc) y_pc <- if ("PC2" %in% pcs) "PC2" else pcs[min(2, length(pcs))]

    x_i <- pc_index(x_pc)
    y_i <- pc_index(y_pc)

    x_lab <- paste0(x_pc, " (", round(v[x_i] * 100, 1), "% variance)")
    y_lab <- paste0(y_pc, " (", round(v[y_i] * 100, 1), "% variance)")

    color_col <- NULL
    if (!is.null(rv$metadata) &&
        !is.null(input$pca_color_col) &&
        input$pca_color_col %in% names(scores)) {
      color_col <- input$pca_color_col
    }

    p <- if (is.null(color_col)) {
      ggplot(scores, aes(x = .data[[x_pc]], y = .data[[y_pc]])) +
        geom_point(size = 3, alpha = 0.9) +
        labs(x = x_lab, y = y_lab)
    } else {
      ggplot(scores, aes(x = .data[[x_pc]], y = .data[[y_pc]], color = .data[[color_col]])) +
        geom_point(size = 3, alpha = 0.9) +
        labs(x = x_lab, y = y_lab, color = color_col)
    }

    if (isTRUE(input$pca_show_labels)) {
      p <- p + geom_text(aes(label = sample_id), vjust = -0.8, size = 3, show.legend = FALSE)
    }

    p + theme_bw() + theme(panel.grid.minor = element_blank())
  })

  output$pca_variance_table <- renderDT({
    if (is.null(rv$features) || length(sample_cols_selected()) < 2) {
      return(datatable(data.frame(note = "Upload/select at least 2 sample columns to run PCA."), options = list(dom = "t")))
    }
    v <- pca_result()$variance
    tbl <- tibble(PC = paste0("PC", seq_along(v)), variance_explained_pct = round(v * 100, 2))
    datatable(head(tbl, 10), options = list(pageLength = 10, scrollX = TRUE))
  })

  output$pca_loadings_table <- renderDT({
    if (is.null(rv$features) || length(sample_cols_selected()) < 2) {
      return(datatable(data.frame(note = "Upload/select at least 2 sample columns to run PCA."), options = list(dom = "t")))
    }

    res <- pca_result()
    loadings <- res$loadings
    pc <- input$loading_pc
    top_n <- input$loading_top_n

    if (!pc %in% names(loadings)) {
      return(datatable(data.frame(note = "Selected PC not found in loadings."), options = list(dom = "t")))
    }

    out <- loadings %>%
      transmute(feature = feature, loading = .data[[pc]], abs_loading = abs(.data[[pc]])) %>%
      arrange(desc(abs_loading)) %>%
      head(top_n)

    datatable(out, options = list(pageLength = min(10, nrow(out)), scrollX = TRUE))
  })

  output$export_summary <- renderText({
    feat_loaded <- !is.null(rv$features)
    meta_loaded <- !is.null(rv$metadata)

    n_feat <- if (feat_loaded) nrow(rv$features) else NA
    n_cols <- if (feat_loaded) ncol(rv$features) else NA

    sample_cols <- if (feat_loaded) sample_cols_selected() else character(0)
    n_samples <- length(sample_cols)

    join_key <- if (meta_loaded && !is.null(input$join_key)) input$join_key else "n/a"
    color_by <- input$pca_color_col %||% "n/a"
    x_pc <- input$pca_x_pc %||% "PC1"
    y_pc <- input$pca_y_pc %||% "PC2"

    miss_max <- input$export_feature_missing_max %||% NA
    nonzero_min <- input$export_feature_nonzero_min %||% NA

    paste(
      paste0("Dataset: ", dataset_stub()),
      paste0("Features table: ", if (feat_loaded) paste0(n_feat, " rows × ", n_cols, " columns") else "not loaded"),
      paste0("Samples detected/selected: ", n_samples),
      paste0("Metadata loaded: ", if (meta_loaded) "yes" else "no"),
      paste0("Join key: ", join_key),
      paste0("PCA: ", x_pc, " vs ", y_pc, " | color by: ", color_by),
      paste0("Filter thresholds: max missing% = ", miss_max, " | min non-zero% = ", nonzero_min),
      sep = "\n"
    )
  })

  output$download_qc_sample <- downloadHandler(
    filename = function() paste0(dataset_stub(), "_qc_by_sample_", timestamp_tag(), ".csv"),
    content  = function(file) readr::write_csv(qc_sample_tbl(), file)
  )

  output$download_qc_feature <- downloadHandler(
    filename = function() paste0(dataset_stub(), "_qc_by_feature_", timestamp_tag(), ".csv"),
    content  = function(file) readr::write_csv(qc_feature_tbl(), file)
  )

  output$download_pca_scores <- downloadHandler(
    filename = function() paste0(dataset_stub(), "_pca_scores_", timestamp_tag(), ".csv"),
    content  = function(file) readr::write_csv(pca_result()$scores, file)
  )

  output$download_pca_loadings <- downloadHandler(
    filename = function() paste0(dataset_stub(), "_pca_top_loadings_", input$loading_pc, "_", timestamp_tag(), ".csv"),
    content  = function(file) {
      res <- pca_result()
      loadings <- res$loadings
      pc <- input$loading_pc
      top_n <- input$loading_top_n

      out <- loadings %>%
        transmute(feature = feature, loading = .data[[pc]], abs_loading = abs(.data[[pc]])) %>%
        arrange(desc(abs_loading)) %>%
        head(top_n)

      readr::write_csv(out, file)
    }
  )

  output$download_filtered_features <- downloadHandler(
    filename = function() paste0(
      dataset_stub(), "_filtered_features_miss", input$export_feature_missing_max,
      "_nz", input$export_feature_nonzero_min, "_", timestamp_tag(), ".csv"
    ),
    content = function(file) {
      req(rv$features)
      sample_cols <- sample_cols_selected()
      req(length(sample_cols) > 0)

      int_mat <- get_intensity_matrix(rv$features, sample_cols)

      miss_pct <- apply(int_mat, 1, function(x) mean(is.na(x)) * 100)
      nonzero_pct <- apply(int_mat, 1, function(x) {
        denom <- sum(!is.na(x))
        if (denom == 0) return(0)
        (sum(x[!is.na(x)] != 0) / denom) * 100
      })

      keep <- (miss_pct <= input$export_feature_missing_max) &
        (nonzero_pct >= input$export_feature_nonzero_min)

      descriptor_cols <- setdiff(names(rv$features), sample_cols)
      out <- rv$features[keep, c(descriptor_cols, sample_cols), drop = FALSE]

      readr::write_csv(out, file)
    }
  )

  output$download_html_report <- downloadHandler(
    filename = function() paste0(dataset_stub(), "_qc_report_", timestamp_tag(), ".html"),
    content = function(file) {
      if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        stop("Package 'rmarkdown' is required for report export. Install it with install.packages('rmarkdown').")
      }

      tmp <- tempfile("qc_report_")
      dir.create(tmp, recursive = TRUE, showWarnings = FALSE)

      qc_sample_path <- file.path(tmp, "qc_by_sample.csv")
      qc_feature_path <- file.path(tmp, "qc_by_feature.csv")
      pca_scores_path <- file.path(tmp, "pca_scores.csv")
      pca_loadings_path <- file.path(tmp, "pca_top_loadings.csv")

      readr::write_csv(qc_sample_tbl(), qc_sample_path)
      readr::write_csv(qc_feature_tbl(), qc_feature_path)
      readr::write_csv(pca_result()$scores, pca_scores_path)

      res <- pca_result()
      loadings <- res$loadings
      pc <- input$loading_pc
      top_n <- input$loading_top_n

      top_loadings <- loadings %>%
        transmute(feature = feature, loading = .data[[pc]], abs_loading = abs(.data[[pc]])) %>%
        arrange(desc(abs_loading)) %>%
        head(top_n)

      readr::write_csv(top_loadings, pca_loadings_path)

      template_path <- "report_template.Rmd"
      if (!file.exists(template_path)) {
        stop("Missing report_template.Rmd. Put it in the same folder as app.R.")
      }

      out_html <- file.path(tmp, "report.html")

      rmarkdown::render(
        input = template_path,
        output_file = out_html,
        params = list(
          report_date = Sys.Date(),
          dataset_name = (input$features_file$name %||% "Demo dataset"),
          qc_sample_csv = qc_sample_path,
          qc_feature_csv = qc_feature_path,
          pca_scores_csv = pca_scores_path,
          pca_loadings_csv = pca_loadings_path,
          pca_color_col = (input$pca_color_col %||% ""),
          pca_x_pc = (input$pca_x_pc %||% "PC1"),
          pca_y_pc = (input$pca_y_pc %||% "PC2"),
          top_n_loadings = top_n
        ),
        envir = new.env(parent = globalenv()),
        quiet = TRUE
      )

      file.copy(out_html, file, overwrite = TRUE)
    }
  )
}

shinyApp(ui, server)
