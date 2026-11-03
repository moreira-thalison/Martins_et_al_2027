#!/usr/bin/env Rscript
# ============================================================================
# 02b_kegg_bubbleplots.R  -  Bubble (dot) plots for the KEGG enrichment produced
# by 02_kegg_enrichment.R (replaces the removed nominal-trend barplots). Reads the
# written *_all.tsv tables (no re-enrichment) and renders, per contrast:
#   (1) ORA bubble  : rows = direction (up / down), top terms by nominal p,
#                     x = fold enrichment, size = overlap Count,
#                     fill = -log10(nominal p), black outline if FDR<=0.05.
#   (2) GSEA bubble : top sets by nominal p, x = NES (sign = direction),
#                     size = setSize, fill = -log10(nominal p), outline if FDR.
# Same encoding as the GO (05b) and COG (05c) bubbles. Up and down are stacked
# as rows (not side-by-side columns) because KEGG terms differ by direction.
#
# Run:  Rscript scripts/02b_kegg_bubbleplots.R      (after 02_kegg_enrichment.R)
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
args <- commandArgs(FALSE); self <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
SDIR <- dirname(self); source(file.path(SDIR, "config.R")); source(file.path(SDIR, "plot_utils.R"))
PADJ  <- CFG$padj
RES   <- CFG$dir_results
TOP_ORA  <- 10L    # terms per direction
TOP_GSEA <- 15L    # sets

# read a written enrichment table if it exists (else NULL)
rd <- function(f) { p <- file.path(RES, f); if (file.exists(p)) fread(p) else NULL }

# cap label length, then wrap (base strwrap; ASCII "..." renders on every device)
tidy_label <- function(x, max = 80, width = 42) {
  x <- ifelse(nchar(x) > max, paste0(substr(x, 1, max - 3), "..."), x)
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
}

# common bubble geom/scale block (fill = -log10 nominal p, outline = FDR pass)
bubble_layers <- function(size_name) {
  list(
    # show.legend = TRUE: draw the FDR key even when no term passes (ggplot2 >= 3.5)
    geom_point(shape = 21, stroke = 0.5, alpha = 0.9, show.legend = TRUE),
    scale_fill_gradient(low = "grey85", high = "#d73027",
                        name = expression(-log[10]~italic(p))),
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey70"),
                        labels = c(`TRUE` = "FDR<=0.05", `FALSE` = "n.s."),
                        name = "BH", drop = FALSE),
    scale_size_continuous(range = c(2, 8), name = size_name),
    guides(colour = guide_legend(override.aes = list(size = 3.5, fill = "white", alpha = 1)),
           fill = guide_colourbar(raster = FALSE)),
    theme_mo(8),
    theme(plot.title.position = "plot",
          panel.grid.major.y = element_line(colour = "grey92"),
          strip.text = element_text(face = "bold"),
          strip.text.y = element_text(angle = 0, hjust = 0))
  )
}
# unique y key per row so ordering holds even if a term recurs across panels
mk_ykey <- function(dt, ord_cols) {
  setorderv(dt, ord_cols)
  dt[, ykey := factor(.I, levels = .I, labels = Description)]
  dt
}

for (ct in CFG$contrasts) {
  cn <- ct$name; tg <- ct$test; cg <- ct$control
  message("=== KEGG bubble plots: ", cn, " ===")

  # ---------------- (1) ORA: up / down ---------------------------------------
  dir_lab <- setNames(sub(" in ", "\nin ", dir_labels(tg)), c("up", "down"))  # 2-line strips
  ora <- rbindlist(lapply(names(dir_lab), function(dir) {
    x <- rd(sprintf("prot__enrich__%s_%s__kegg_ora_all.tsv", cn, dir))
    if (is.null(x) || !nrow(x)) return(NULL)
    x[, direction := dir_lab[[dir]]]
    x[order(pvalue)][seq_len(min(TOP_ORA, .N))]
  }), fill = TRUE)

  if (nrow(ora)) {
    ora[, direction := factor(direction, levels = dir_lab)]
    ora[, sig := factor(p.adjust <= PADJ, levels = c(TRUE, FALSE))]
    ora[, Description := tidy_label(Description)]
    ora <- mk_ykey(ora, c("direction", "FoldEnrichment"))
    p <- ggplot(ora, aes(FoldEnrichment, ykey, fill = -log10(pvalue), size = Count, colour = sig)) +
      geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
      bubble_layers("overlap\n(Count)") +
      facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
      labs(title = sprintf("KEGG ORA - %s", cn),
           subtitle = sprintf("top %d pathways / direction by nominal p; universe = detected proteome\nfill = -log10 nominal p; black outline = BH FDR<=0.05",
                              TOP_ORA),
           x = "fold enrichment", y = NULL)
    save_fig(p, "enrich", cn, "kegg_bubble_ora", CFG,
             width_mm = CFG$width_double_mm * 0.85, height_mm = 35 + 6.2 * nrow(ora))
    message(sprintf("  ORA bubble: %d terms (%d pass FDR)", nrow(ora), sum(ora$sig == "TRUE")))
  } else message("  [empty] no KEGG ORA tables for ", cn)

  # ---------------- (2) GSEA (NES sign = direction) --------------------------
  gse <- rd(sprintf("prot__enrich__%s__kegg_gsea_all.tsv", cn))
  if (!is.null(gse) && nrow(gse)) {
    gse <- gse[!is.na(pvalue)][order(pvalue)][seq_len(min(TOP_GSEA, .N))]
    gse[, sig := factor(p.adjust <= PADJ, levels = c(TRUE, FALSE))]
    gse[, Description := tidy_label(Description)]
    gse <- mk_ykey(gse, "NES")
    p <- ggplot(gse, aes(NES, ykey, fill = -log10(pvalue), size = setSize, colour = sig)) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
      bubble_layers("set size") +
      labs(title = sprintf("KEGG GSEA - %s", cn),
           subtitle = sprintf("top %d pathways by nominal p; ranked by limma log2FC\nNES>0 = %s, NES<0 = %s; black outline = BH FDR<=0.05",
                              TOP_GSEA, dir_labels(tg, lower = TRUE)[["Up"]], dir_labels(tg, lower = TRUE)[["Down"]]),
           x = "normalised enrichment score (NES)", y = NULL)
    save_fig(p, "enrich", cn, "kegg_bubble_gsea", CFG,
             width_mm = CFG$width_double_mm * 0.85, height_mm = 40 + 6.2 * nrow(gse))
    message(sprintf("  GSEA bubble: %d sets (%d pass FDR)", nrow(gse), sum(gse$sig == "TRUE")))
  } else message("  [empty] no KEGG GSEA table for ", cn)
}

message("\nDONE: KEGG bubble plots")
