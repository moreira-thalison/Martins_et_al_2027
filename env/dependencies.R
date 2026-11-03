#!/usr/bin/env Rscript
# ============================================================================
# env/dependencies.R  -  Idempotent dependency installer for the AF293 DIA
# proteomics pipeline. Installs any MISSING CRAN + Bioconductor packages via
# BiocManager (which resolves both repositories). Safe to re-run.
#
#   Rscript env/dependencies.R          # install what is missing
#   Rscript env/dependencies.R --check  # only report status, install nothing
#
# Developed & last run on R 4.5.2 (see env/sessionInfo.txt for the exact
# package versions of a known-good run).
# ============================================================================
check_only <- "--check" %in% commandArgs(TRUE)

PKGS <- c(
  # ---- CRAN ----
  "data.table", "ggplot2", "matrixStats", "pheatmap", "ggrepel", "UpSetR",
  "circlize", "VennDiagram", "futile.logger", "tidygraph", "dplyr", "ggraph",
  "igraph", "patchwork", "scales", "stringr", "RCurl", "ggforce", "svglite", "ragg",
  # ---- Bioconductor ----
  "limma", "clusterProfiler", "enrichplot", "KEGGREST", "ComplexHeatmap",
  "EnhancedVolcano", "ggkegg", "STRINGdb", "fgsea"
)

have <- vapply(PKGS, requireNamespace, logical(1), quietly = TRUE)
cat(sprintf("R %s | %d/%d dependencies present\n",
            getRversion(), sum(have), length(PKGS)))
missing <- PKGS[!have]
if (length(missing)) cat("Missing:", paste(missing, collapse = ", "), "\n")

if (check_only) quit(save = "no", status = if (length(missing)) 1L else 0L)

if (length(missing)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(missing, update = FALSE, ask = FALSE)
  still <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still)) stop("Failed to install: ", paste(still, collapse = ", "))
  cat("All dependencies installed.\n")
} else {
  cat("Nothing to install.\n")
}
