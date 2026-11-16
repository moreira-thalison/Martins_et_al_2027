#!/usr/bin/env Rscript
# =============================================================================
# 02b_kegg_bubbleplots.R
# Bubble plots for the KEGG enrichment of 02_kegg_enrichment.R, read from its
# *_all.tsv tables (no re-enrichment). Per contrast:
#   ORA  : rows = direction (up / down), top terms by nominal p,
#          x = fold enrichment, size = overlap Count
#   GSEA : top sets by nominal p, x = NES (sign = direction), size = setSize
# Both: fill = -log10(nominal p), black outline = BH FDR <= 0.05 (same encoding
# as the GO 05b and COG 05c bubbles). Up and down are stacked as rows because
# KEGG terms differ by direction.
#
# Input  : results/prot__enrich__<contrast>*__kegg_{ora,gsea}_all.tsv
# Outputs: figures/prot__enrich__<contrast>__kegg_bubble_{ora,gsea}.{png,pdf}
# Run    : Rscript scripts/02b_kegg_bubbleplots.R   (after 02)
#
# Steps
#   0. Setup
#   1. Bubble plots per contrast
#      1.1 ORA bubble
#      1.2 GSEA bubble
# =============================================================================


# 0. Setup ---------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
SDIR <- if (length(script_file)) dirname(normalizePath(script_file)) else "scripts"
source(file.path(SDIR, "config.R"))
source(file.path(SDIR, "plot_utils.R"))
source(file.path(SDIR, "enrich_utils.R"))
PADJ     <- CFG$padj
TOP_ORA  <- 10L    # terms per direction
TOP_GSEA <- 15L    # sets


# 1. Bubble plots per contrast -------------------------------------------------
for (ct in CFG$contrasts) {
  cn <- ct$name
  tg <- ct$test
  dl_low <- dir_labels(tg, lower = TRUE)
  message("=== KEGG bubble plots: ", cn, " ===")

  ## 1.1 ORA bubble ----
  dir_lab <- setNames(sub(" in ", "\nin ", dir_labels(tg)), c("up", "down"))   # 2-line strips
  ora <- rbindlist(lapply(names(dir_lab), function(dir) {
    x <- read_result_table(sprintf("prot__enrich__%s_%s__kegg_ora_all.tsv", cn, dir), CFG)
    if (is.null(x) || !nrow(x)) return(NULL)
    x[, direction := dir_lab[[dir]]]
    x[order(pvalue)][seq_len(min(TOP_ORA, .N))]
  }), fill = TRUE)

  if (nrow(ora)) {
    ora[, direction := factor(direction, levels = dir_lab)]
    ora[, sig := factor(p.adjust <= PADJ, levels = c(TRUE, FALSE))]
    ora[, Description := tidy_label(Description)]
    ora <- order_ykey(ora, c("direction", "FoldEnrichment"))
    p <- ggplot(ora, aes(FoldEnrichment, ykey, fill = -log10(pvalue), size = Count, colour = sig)) +
      geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
      bubble_layers("overlap\n(Count)") +
      facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
      labs(title = sprintf("KEGG ORA - %s", cn),
           subtitle = sprintf(paste0("top %d pathways / direction by nominal p; universe = detected proteome\n",
                                     "fill = -log10 nominal p; black outline = BH FDR<=0.05"), TOP_ORA),
           x = "fold enrichment", y = NULL) +
      enrich_theme()
    save_fig(p, "enrich", cn, "kegg_bubble_ora", CFG,
             width_mm = CFG$width_double_mm * 0.85, height_mm = 35 + 6.2 * nrow(ora))
    message(sprintf("  ORA bubble: %d terms (%d pass FDR)", nrow(ora), sum(ora$sig == "TRUE")))
  } else {
    message("  [empty] no KEGG ORA tables for ", cn)
  }

  ## 1.2 GSEA bubble ----
  gse <- read_result_table(sprintf("prot__enrich__%s__kegg_gsea_all.tsv", cn), CFG)
  if (!is.null(gse) && nrow(gse)) {
    gse <- gse[!is.na(pvalue)][order(pvalue)][seq_len(min(TOP_GSEA, .N))]
    gse[, sig := factor(p.adjust <= PADJ, levels = c(TRUE, FALSE))]
    gse[, Description := tidy_label(Description)]
    gse <- order_ykey(gse, "NES")
    p <- ggplot(gse, aes(NES, ykey, fill = -log10(pvalue), size = setSize, colour = sig)) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
      bubble_layers("set size") +
      labs(title = sprintf("KEGG GSEA - %s", cn),
           subtitle = sprintf(paste0("top %d pathways by nominal p; ranked by limma log2FC\n",
                                     "NES>0 = %s, NES<0 = %s; black outline = BH FDR<=0.05"),
                              TOP_GSEA, dl_low[["Up"]], dl_low[["Down"]]),
           x = "normalised enrichment score (NES)", y = NULL) +
      enrich_theme()
    save_fig(p, "enrich", cn, "kegg_bubble_gsea", CFG,
             width_mm = CFG$width_double_mm * 0.85, height_mm = 40 + 6.2 * nrow(gse))
    message(sprintf("  GSEA bubble: %d sets (%d pass FDR)", nrow(gse), sum(gse$sig == "TRUE")))
  } else {
    message("  [empty] no KEGG GSEA table for ", cn)
  }
}

message("\nDONE: KEGG bubble plots")
