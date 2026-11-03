#!/usr/bin/env Rscript
# ============================================================================
# 05b_go_bubbleplots.R  -  Bubble (dot) plots for the GO enrichment produced by
# 05_go_enrichment.R. Reads the written *_all.tsv tables (no re-enrichment) and
# renders, per contrast:
#   (1) ORA bubble  : facet aspect (BP/MF/CC rows) x direction (up/down cols),
#                     x = fold enrichment, size = overlap Count,
#                     fill = -log10(nominal p), black outline if FDR<=0.05.
#   (2) GSEA bubble : facet aspect (rows), x = NES (sign = direction),
#                     size = setSize, fill = -log10(nominal p), outline if FDR.
# Exploratory-honest: colour is the nominal trend; the FDR outline marks what
# actually survives BH (same padj<=0.05 as the rest of the pipeline).
#
# Run:  Rscript scripts/05b_go_bubbleplots.R      (after 05_go_enrichment.R)
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
args <- commandArgs(FALSE); self <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
SDIR <- dirname(self); source(file.path(SDIR, "config.R")); source(file.path(SDIR, "plot_utils.R"))
PADJ <- CFG$padj
RES  <- CFG$dir_results
ASP_LEVELS <- c(BP = "Biological process", MF = "Molecular function", CC = "Cellular component")

# read a written enrichment table if it exists (else NULL)
rd <- function(f) { p <- file.path(RES, f); if (file.exists(p)) fread(p) else NULL }

# keep y-axis labels from blowing out the panel: cap length, then wrap to lines
# <= `width` chars (base strwrap, no stringr dependency). Long GO names (e.g. the
# rRNA-processing terms) otherwise consume the whole figure width.
tidy_label <- function(x, max = 80, width = 42) {
  x <- ifelse(nchar(x) > max, paste0(substr(x, 1, max - 3), "..."), x)
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
}

# common bubble geom/scale block (fill = -log10 nominal p, outline = FDR pass)
bubble_layers <- function(size_name) {
  list(
    geom_point(shape = 21, stroke = 0.5, alpha = 0.9),
    scale_fill_gradient(low = "grey85", high = "#d73027",
                        name = expression(-log[10]~italic(p))),
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey70"),
                        labels = c(`TRUE` = "FDR<=0.05", `FALSE` = "n.s."),
                        name = "BH"),
    scale_size_continuous(range = c(2, 8), name = size_name),
    guides(colour = guide_legend(override.aes = list(size = 3.5))),
    theme_mo(8),
    theme(panel.grid.major.y = element_line(colour = "grey92"),
          strip.text = element_text(face = "bold"))
  )
}
# unique per-facet y key so ordering is correct even if a term recurs across panels
mk_ykey <- function(dt, ord_cols) {
  setorderv(dt, ord_cols)
  dt[, ykey := factor(.I, levels = .I, labels = Description)]
  dt
}

for (ct in CFG$contrasts) {
  cn <- ct$name; tg <- ct$test; cg <- ct$control
  message("=== bubble plots: ", cn, " ===")

  # ---------------- (1) ORA: up / down, per aspect -------------------------
  ora <- rbindlist(lapply(names(ASP_LEVELS), function(asp) {
    ak <- tolower(asp)
    rbindlist(lapply(c("up", "down"), function(dir) {
      x <- rd(sprintf("prot__enrich__%s_%s__go_%s_ora_all.tsv", cn, dir, ak))
      if (is.null(x) || !nrow(x)) return(NULL)
      x[, `:=`(aspect = ASP_LEVELS[[asp]],
               direction = if (dir == "up") dir_labels(tg)[["Up"]] else dir_labels(tg)[["Down"]])]
      x[order(pvalue)][seq_len(min(6L, .N))]           # top 6 per aspect x direction
    }), fill = TRUE)
  }), fill = TRUE)

  if (nrow(ora)) {
    ora[, aspect := factor(aspect, levels = ASP_LEVELS)]
    ora[, direction := factor(direction,
        levels = unname(dir_labels(tg)))]
    ora[, sig := p.adjust <= PADJ]
    ora[, neglogp := -log10(pvalue)]
    ora[, Description := tidy_label(Description)]
    ora <- mk_ykey(ora, c("aspect", "direction", "FoldEnrichment"))
    p <- ggplot(ora, aes(FoldEnrichment, ykey, fill = neglogp, size = Count, colour = sig)) +
      geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
      bubble_layers("overlap\n(Count)") +
      facet_grid(aspect ~ direction, scales = "free_y", space = "free_y") +
      labs(title = sprintf("GO ORA - %s", cn),
           subtitle = "top 6 terms / panel by nominal p; fill = -log10 p; outline = BH FDR<=0.05",
           x = "fold enrichment", y = NULL)
    save_fig(p, "enrich", cn, "go_bubble_ora", CFG,
             width_mm = CFG$width_double_mm, height_mm = CFG$height_mm * 2.4)
    message(sprintf("  ORA bubble: %d terms (%d pass FDR)", nrow(ora), sum(ora$sig)))
  } else message("  [empty] no ORA tables for ", cn)

  # ---------------- (2) GSEA: per aspect (NES sign = direction) ------------
  gse <- rbindlist(lapply(names(ASP_LEVELS), function(asp) {
    ak <- tolower(asp)
    x <- rd(sprintf("prot__enrich__%s__go_%s_gsea_all.tsv", cn, ak))
    if (is.null(x) || !nrow(x)) return(NULL)
    x[, aspect := ASP_LEVELS[[asp]]]
    x[order(pvalue)][seq_len(min(8L, .N))]              # top 8 per aspect
  }), fill = TRUE)

  if (nrow(gse)) {
    gse[, aspect := factor(aspect, levels = ASP_LEVELS)]
    gse[, sig := p.adjust <= PADJ]
    gse[, neglogp := -log10(pvalue)]
    gse[, Description := tidy_label(Description)]
    gse <- mk_ykey(gse, c("aspect", "NES"))
    p <- ggplot(gse, aes(NES, ykey, fill = neglogp, size = setSize, colour = sig)) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
      bubble_layers("set size") +
      facet_grid(aspect ~ ., scales = "free_y", space = "free_y") +
      labs(title = sprintf("GO GSEA - %s", cn),
           subtitle = sprintf("top 8 sets / aspect by nominal p; NES>0 = %s, NES<0 = %s; outline = BH FDR<=0.05",
                              dir_labels(tg, lower = TRUE)[["Up"]], dir_labels(tg, lower = TRUE)[["Down"]]),
           x = "normalised enrichment score (NES)", y = NULL)
    save_fig(p, "enrich", cn, "go_bubble_gsea", CFG,
             width_mm = CFG$width_double_mm, height_mm = CFG$height_mm * 2.4)
    message(sprintf("  GSEA bubble: %d sets (%d pass FDR)", nrow(gse), sum(gse$sig)))
  } else message("  [empty] no GSEA tables for ", cn)
}

message("\nDONE: GO bubble plots")
