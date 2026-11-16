#!/usr/bin/env Rscript
# =============================================================================
# plot_utils.R
# Shared figure helpers: theme, saving, QC / DE plot templates and heatmap
# helpers. Every figure is saved as PNG (figures/) plus PDF (figures/Figures
# PDF/) at CFG$dpi. Sourced after config.R.
#
# Sections
#   1. Theme and units       theme_mo(), mm2in(), svg_inline_styles()
#   2. Saving figures        fig_path(), save_fig(), save_heatmap(),
#                            save_complex_heatmap()
#   3. QC plots              pca_scores(), plot_pca(), plot_corr_heatmap()
#   4. DE plots              plot_volcano(), plot_ma(), plot_pval_hist(),
#                            plot_count_bar(), plot_overlap()
#   5. Heatmap helpers       row_zscore(), status_ordered_proteins()
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })
.have <- function(pkg) requireNamespace(pkg, quietly = TRUE)


# 1. Theme and units -----------------------------------------------------------
theme_mo <- function(base = 9) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = base + 1),
          legend.key.size = unit(3.5, "mm"))
}

mm2in <- function(mm) mm / 25.4

# svglite keeps the default shape style (fill: none; stroke: black; round caps
# and joins) in a CSS <style> block. Affinity and some other editors ignore that
# block, so unstyled lines and outlines fall back to the SVG default (solid
# black fill). This writes the defaults into every shape's own style attribute.
svg_inline_styles <- function(path) {
  x <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  if (grepl("<g class='glyphgroup'", x, fixed = TRUE))    # real groups, not the CSS rule
    warning("SVG has glyph groups; their paths get fill: none too: ", path)
  defaults <- c("fill" = "none", "stroke" = "#000000", "stroke-linecap" = "round",
                "stroke-linejoin" = "round", "stroke-miterlimit" = "10.00")
  fix_tag <- function(tag) {
    has_style <- grepl("style='", tag, fixed = TRUE)
    style <- if (has_style) sub("^.*style='([^']*)'.*$", "\\1", tag) else ""
    add <- defaults[!vapply(names(defaults), function(k)
      grepl(sprintf("(^|;)\\s*%s\\s*:", k), style, perl = TRUE), logical(1))]
    if (!length(add)) return(tag)
    extra <- paste0(names(add), ": ", add, ";", collapse = " ")
    if (has_style) sub("style='([^']*)'", paste0("style='\\1 ", extra, "'"), tag)
    else sub("^(<[a-z]+)", paste0("\\1 style='", extra, "'"), tag)
  }
  m <- gregexpr("<(line|polyline|polygon|path|rect|circle)\\b[^>]*>", x, perl = TRUE)
  regmatches(x, m) <- list(vapply(regmatches(x, m)[[1]], fix_tag, character(1), USE.NAMES = FALSE))
  writeLines(x, path, useBytes = TRUE)
  invisible(path)
}


# 2. Saving figures ------------------------------------------------------------
# PNG goes to figures/, PDF to figures/Figures PDF/.
fig_path <- function(file_base, ext, CFG)
  file.path(if (ext == "pdf") CFG$dir_figures_pdf else CFG$dir_figures,
            paste0(file_base, ".", ext))

# ggplot -> PNG + PDF named prot__<stage>__<key>__<name>.
save_fig <- function(p, stage, key, name, CFG,
                     width_mm = CFG$width_single_mm, height_mm = CFG$height_mm) {
  fb <- sprintf("prot__%s__%s__%s", stage, key, name)
  for (ext in c("png", "pdf"))
    ggsave(fig_path(fb, ext, CFG), p, width = mm2in(width_mm),
           height = mm2in(height_mm), dpi = CFG$dpi, limitsize = FALSE)
  invisible(fb)
}

# pheatmap object -> PNG + PDF named prot__<stage>__<key>__<name>.
save_heatmap <- function(ph, stage, key, name, CFG,
                         width_mm = CFG$width_single_mm, height_mm = CFG$height_mm) {
  fb <- sprintf("prot__%s__%s__%s", stage, key, name)
  png(fig_path(fb, "png", CFG), width = width_mm, height = height_mm, units = "mm", res = CFG$dpi)
  grid::grid.newpage(); grid::grid.draw(ph$gtable); invisible(dev.off())
  pdf(fig_path(fb, "pdf", CFG), width = mm2in(width_mm), height = mm2in(height_mm))
  grid::grid.newpage(); grid::grid.draw(ph$gtable); invisible(dev.off())
  invisible(fb)
}

# ComplexHeatmap -> PNG + PDF named <name>, legends merged on the right.
save_complex_heatmap <- function(ht, name, width_mm, height_mm, CFG) {
  draw_ht <- function()
    ComplexHeatmap::draw(ht, merge_legends = TRUE, heatmap_legend_side = "right",
                         annotation_legend_side = "right")
  png(fig_path(name, "png", CFG), width = width_mm, height = height_mm, units = "mm", res = CFG$dpi)
  draw_ht(); invisible(dev.off())
  pdf(fig_path(name, "pdf", CFG), width = mm2in(width_mm), height = mm2in(height_mm))
  draw_ht(); invisible(dev.off())
  invisible(name)
}


# 3. QC plots ------------------------------------------------------------------
# PCA of samples. Column means fill missing values for the projection only
# (the DE input is untouched); zero-variance proteins are dropped.
# Returns list(scores = data.frame(PC1, PC2, sample), var_pct).
pca_scores <- function(mat) {
  m <- mat
  for (j in seq_len(ncol(m))) m[is.na(m[, j]), j] <- mean(m[, j], na.rm = TRUE)
  m  <- m[matrixStats::rowVars(m) > 0, , drop = FALSE]
  pc <- prcomp(t(m), scale. = TRUE)
  list(scores  = data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], sample = colnames(mat)),
       var_pct = 100 * pc$sdev^2 / sum(pc$sdev^2))
}

plot_pca <- function(mat, group, cols, title = "Proteome PCA") {
  pca <- pca_scores(mat)
  d <- pca$scores
  d$group <- group
  ggplot(d, aes(PC1, PC2, colour = group, label = sample)) +
    geom_point(size = 2.4) +
    ggrepel::geom_text_repel(size = 2.3, show.legend = FALSE, max.overlaps = 20) +
    scale_colour_manual(values = cols) +
    labs(title = title, x = sprintf("PC1 (%.1f%%)", pca$var_pct[1]),
         y = sprintf("PC2 (%.1f%%)", pca$var_pct[2])) +
    theme_mo()
}

plot_corr_heatmap <- function(mat, group, cols, title = "Sample correlation") {
  cm  <- cor(mat, use = "pairwise.complete.obs")
  ann <- data.frame(group = group)
  rownames(ann) <- colnames(mat)
  pheatmap::pheatmap(cm, main = title, annotation_col = ann,
                     annotation_colors = list(group = cols[names(cols) %in% group]),
                     show_colnames = FALSE, silent = TRUE, fontsize = 7,
                     color = grDevices::colorRampPalette(c("#4575b4", "white", "#d73027"))(50))
}


# 4. DE plots ------------------------------------------------------------------
# Volcano: faded NS background, significant points on top, y-axis trimmed to
# the data, repelled labels on the top n_label hits.
plot_volcano <- function(d, contrast, padj_cut, cols = REG_COLORS,
                         n_label = 10L, label_size = 2.6) {
  if (!.have("ggrepel")) stop("ggrepel is required for volcano labels")

  ## 4.1 Plot data ----
  d <- as.data.table(copy(d))[is.finite(log2FC) & is.finite(padj)]
  d[, neglog10 := -log10(pmax(padj, .Machine$double.xmin))]
  d[, label_plot := label]
  d[is.na(label_plot) | !nzchar(label_plot), label_plot := protein_id]
  d[, status := direction_status(regulation, "NS")]
  setorder(d, -status)                                   # NS drawn first (background)

  status_cols <- setNames(unname(cols[c("Up", "Down", "NS")]), levels(d$status))
  x_lim  <- ceiling(max(abs(d$log2FC)) * 2) / 2          # symmetric, snapped to 0.5
  y_top  <- max(d$neglog10) * 1.12                       # trim empty upper space
  n_up   <- sum(d$regulation == "Up")
  n_down <- sum(d$regulation == "Down")
  top    <- d[regulation != "NS"][order(padj, -abs(log2FC))][seq_len(min(n_label, .N))]
  dl_low <- dir_labels(lower = TRUE)

  ## 4.2 Plot ----
  ggplot(d, aes(log2FC, neglog10)) +
    geom_vline(xintercept = 0, linetype = "22", linewidth = 0.3, colour = "grey60") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "22", linewidth = 0.3, colour = "grey60") +
    geom_point(aes(colour = status, size = status, alpha = status), stroke = 0) +
    ggrepel::geom_text_repel(
      data = top, aes(label = label_plot), size = label_size, colour = "grey15",
      segment.size = 0.25, segment.colour = "grey60", min.segment.length = 0,
      box.padding = 0.5, point.padding = 0.25, max.overlaps = Inf, seed = 42L) +
    scale_colour_manual(values = status_cols, name = NULL, breaks = names(status_cols)) +
    scale_size_manual(values = setNames(c(2.0, 2.0, 0.8), names(status_cols)), guide = "none") +
    scale_alpha_manual(values = setNames(c(0.9, 0.9, 0.3), names(status_cols)), guide = "none") +
    coord_cartesian(xlim = c(-x_lim, x_lim), ylim = c(0, y_top), clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0))) +
    guides(colour = guide_legend(override.aes = list(size = 2.8, alpha = 1))) +
    labs(title = sprintf("Volcano - %s", gsub("_", " ", contrast, fixed = TRUE)),
         subtitle = sprintf("BH padj <= %.2g; no log2FC cutoff", padj_cut),
         caption = sprintf("%d %s | %d %s", n_up, dl_low[["Up"]], n_down, dl_low[["Down"]]),
         x = expression(Log[2]~"fold change"),
         y = expression(-Log[10]~"adjusted p-value")) +
    theme_mo(9) +
    theme(legend.position = "right",
          plot.subtitle = element_text(colour = "grey30"),
          plot.caption  = element_text(colour = "grey35"),
          plot.margin   = margin(6, 10, 6, 6))
}

plot_ma <- function(d, contrast, padj_cut, cols = REG_COLORS) {
  d  <- as.data.table(d)[is.finite(log2FC) & is.finite(AveExpr)]
  dl <- dir_labels()
  ggplot(d, aes(AveExpr, log2FC, colour = regulation)) +
    geom_point(size = 1, alpha = 0.7) +
    scale_colour_manual(values = cols, name = NULL,
                        labels = c(Up = dl[["Up"]], Down = dl[["Down"]], NS = "NS")) +
    geom_hline(yintercept = 0, linetype = 3, colour = "grey60") +
    labs(title = sprintf("MA - %s", contrast),
         x = "mean log2 abundance", y = "log2 fold change") +
    theme_mo()
}

plot_pval_hist <- function(d, contrast) {
  d <- as.data.table(d)[is.finite(pvalue)]
  ggplot(d, aes(pvalue)) +
    geom_histogram(breaks = seq(0, 1, 0.05), fill = "#4575b4", colour = "white") +
    labs(title = sprintf("p-value histogram - %s", contrast),
         x = "raw p-value", y = "proteins") +
    theme_mo()
}

plot_count_bar <- function(summ_dt, title) {
  m <- melt(summ_dt[, .(contrast, Up = up, Down = down)],
            id.vars = "contrast", variable.name = "regulation", value.name = "n")
  ggplot(m, aes(contrast, n, fill = regulation)) +
    geom_col(position = position_dodge(0.7), width = 0.65) +
    geom_text(aes(label = n), position = position_dodge(0.7), vjust = -0.3, size = 2.4) +
    scale_fill_manual(values = REG_COLORS[c("Up", "Down")], name = NULL,
                      labels = c(Up = dir_labels()[["Up"]], Down = dir_labels()[["Down"]])) +
    labs(title = title, x = NULL, y = "significant proteins") +
    theme_mo() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

# UpSet plot of DE proteins shared between contrasts (needs >= 2 non-empty sets).
plot_overlap <- function(sets, title) {
  sets <- sets[vapply(sets, length, 0L) > 0]
  if (length(sets) < 2) return(invisible(NULL))
  UpSetR::upset(UpSetR::fromList(sets), order.by = "freq",
                mainbar.y.label = "shared DE proteins", sets.x.label = "DE per contrast")
}


# 5. Heatmap helpers -----------------------------------------------------------
# Per-protein (row) z-score, clipped to [-clip, clip]. Constant rows become 0
# instead of NaN.
row_zscore <- function(m, clip = Inf) {
  row_mean <- rowMeans(m, na.rm = TRUE)
  row_sd   <- apply(m, 1, sd, na.rm = TRUE)
  z  <- sweep(m, 1, row_mean, "-")
  ok <- is.finite(row_sd) & row_sd > 0
  z[ok, ]  <- sweep(z[ok, , drop = FALSE], 1, row_sd[ok], "/")
  z[!ok, ] <- 0
  z[z < -clip] <- -clip
  z[z >  clip] <-  clip
  z
}

# Every quantified protein of the DE table, with `status` = Higher / Lower in
# <test> / <ns_label>, ordered by status then padj (row order of the
# all-protein heatmap in 08 and 09b).
status_ordered_proteins <- function(de, protein_ids, ns_label = "Not significant") {
  d <- de[protein_id %chin% protein_ids]
  d[, status := direction_status(regulation, ns_label)]
  setorder(d, status, padj)
  d[]
}
