#!/usr/bin/env Rscript
# =============================================================================
# 05b_go_bubbleplots.R
# Bubble plots for the GO enrichment of 05_go_enrichment.R, read from its
# *_all.tsv tables (no re-enrichment). Per contrast:
#   ORA  : facets aspect (BP / MF / CC rows) x direction (up / down columns),
#          x = fold enrichment, size = overlap Count
#   GSEA : facets aspect (rows), x = NES (sign = direction), size = setSize
# Both: fill = -log10(nominal p) shows the nominal trend; black outline marks
# what survives BH FDR <= 0.05 (same padj as the rest of the pipeline).
#
# Input  : results/prot__enrich__<contrast>*__go_<aspect>_{ora,gsea}_all.tsv
# Outputs: figures/prot__enrich__<contrast>__go_bubble_{ora,gsea}.{png,pdf}
# Run    : Rscript scripts/05b_go_bubbleplots.R   (after 05)
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
PADJ <- CFG$padj
TOP_ORA  <- 6L     # terms per aspect x direction
TOP_GSEA <- 8L     # sets per aspect

# GO keeps vertical aspect strips (unlike enrich_theme(), used by KEGG / COG).
go_theme <- list(theme_mo(8),
                 theme(panel.grid.major.y = element_line(colour = "grey92"),
                       strip.text = element_text(face = "bold")))


# 1. Bubble plots per contrast -------------------------------------------------
for (ct in CFG$contrasts) {
  cn <- ct$name
  tg <- ct$test
  dl <- dir_labels(tg)
  dl_low <- dir_labels(tg, lower = TRUE)
  message("=== bubble plots: ", cn, " ===")

  ## 1.1 ORA bubble ----
  dir_lab <- setNames(unname(dl), c("up", "down"))
  ora <- rbindlist(lapply(names(GO_ASPECT_LABELS), function(asp) {
    rbindlist(lapply(names(dir_lab), function(dir) {
      x <- read_result_table(sprintf("prot__enrich__%s_%s__go_%s_ora_all.tsv", cn, dir, tolower(asp)), CFG)
      if (is.null(x) || !nrow(x)) return(NULL)
      x[, `:=`(aspect = GO_ASPECT_LABELS[[asp]], direction = dir_lab[[dir]])]
      x[order(pvalue)][seq_len(min(TOP_ORA, .N))]
    }), fill = TRUE)
  }), fill = TRUE)

  if (nrow(ora)) {
    ora[, aspect := factor(aspect, levels = GO_ASPECT_LABELS)]
    ora[, direction := factor(direction, levels = unname(dl))]
    ora[, sig := p.adjust <= PADJ]
    ora[, neglogp := -log10(pvalue)]
    ora[, Description := tidy_label(Description)]
    ora <- order_ykey(ora, c("aspect", "direction", "FoldEnrichment"))
    p <- ggplot(ora, aes(FoldEnrichment, ykey, fill = neglogp, size = Count, colour = sig)) +
      geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
      bubble_layers("overlap\n(Count)") +
      facet_grid(aspect ~ direction, scales = "free_y", space = "free_y") +
      labs(title = sprintf("GO ORA - %s", cn),
           subtitle = sprintf("top %d terms / panel by nominal p; fill = -log10 p; outline = BH FDR<=0.05",
                              TOP_ORA),
           x = "fold enrichment", y = NULL) +
      go_theme
    save_fig(p, "enrich", cn, "go_bubble_ora", CFG,
             width_mm = CFG$width_double_mm, height_mm = CFG$height_mm * 2.4)
    message(sprintf("  ORA bubble: %d terms (%d pass FDR)", nrow(ora), sum(ora$sig)))
  } else {
    message("  [empty] no ORA tables for ", cn)
  }

  ## 1.2 GSEA bubble ----
  gse <- go_gsea_top(cn, TOP_GSEA, CFG)
  if (nrow(gse)) {
    p <- ggplot(gse, aes(NES, ykey, fill = neglogp, size = setSize, colour = sig)) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
      bubble_layers("set size") +
      facet_grid(aspect ~ ., scales = "free_y", space = "free_y") +
      labs(title = sprintf("GO GSEA - %s", cn),
           subtitle = sprintf("top %d sets / aspect by nominal p; NES>0 = %s, NES<0 = %s; outline = BH FDR<=0.05",
                              TOP_GSEA, dl_low[["Up"]], dl_low[["Down"]]),
           x = "normalised enrichment score (NES)", y = NULL) +
      go_theme
    save_fig(p, "enrich", cn, "go_bubble_gsea", CFG,
             width_mm = CFG$width_double_mm, height_mm = CFG$height_mm * 2.4)
    message(sprintf("  GSEA bubble: %d sets (%d pass FDR)", nrow(gse), sum(gse$sig)))
  } else {
    message("  [empty] no GSEA tables for ", cn)
  }
}

message("\nDONE: GO bubble plots")
