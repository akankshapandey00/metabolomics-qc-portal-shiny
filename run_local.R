# One-click launcher for the Metabolomics QC Portal (Shiny)
# Run this from the project root:
#   source("run_local.R")
#
# Tip: If you double-click this file in RStudio, it’ll open and you can click “Source”.

options(shiny.maxRequestSize = 300 * 1024^2)

needed <- c("shiny","DT","readr","dplyr","stringr","tibble","ggplot2","rmarkdown","knitr")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) > 0) {
  message("Missing packages: ", paste(missing, collapse = ", "))
  message("Installing… (this is one-time)")
  install.packages(missing)
}

shiny::runApp(launch.browser = TRUE)
