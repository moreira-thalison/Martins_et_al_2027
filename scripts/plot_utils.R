#!/usr/bin/env Rscript
# ============================================================================
# plot_utils.R  -  Self-contained QC/DE plot templates (ggplot2, pheatmap,
# EnhancedVolcano),
# reproducing the reference multi-omics figure set for this A. fumigatus run.
# All figures saved as PNG (figures/) + PDF (figures/Figures PDF/) at the configured DPI.
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(data.table)
})
.have <- function(p) requireNamespace(p, quietly = TRUE)

theme_mo <- function(base = 9) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = base + 1),
          legend.key.size = unit(3.5, "mm"))
}

mm2in <- function(mm) mm / 25.4

# output path for a figure: PNG in figures/, PDF in figures/Figures PDF/
fig_path <- function(file_base, ext, CFG)
  file.path(if (ext == "pdf") CFG$dir_figures_pdf else CFG$dir_figures,
            paste0(file_base, ".", ext))

# save a ggplot as PNG + PDF using the standard naming prot__<stage>__<key>__<name>
save_fig <- function(p, stage, key, name, CFG,
                     width_mm = CFG$width_single_mm, height_mm = CFG$height_mm) {
  fb <- sprintf("prot__%s__%s__%s", stage, key, name)
  for (ext in c("png", "pdf"))
    ggsave(fig_path(fb, ext, CFG), p, width = mm2in(width_mm),
           height = mm2in(height_mm), dpi = CFG$dpi, limitsize = FALSE)
  invisible(fb)
}

save_heatmap <- function(ph, stage, key, name, CFG, width_mm = CFG$width_single_mm,
                         height_mm = CFG$height_mm) {
  fb <- sprintf("prot__%s__%s__%s", stage, key, name)
  png(fig_path(fb, "png", CFG), width = width_mm, height = height_mm, units = "mm", res = CFG$dpi)
  grid::grid.newpage(); grid::grid.draw(ph$gtable); dev.off()
  pdf(fig_path(fb, "pdf", CFG), width = mm2in(width_mm), height = mm2in(height_mm))
  grid::grid.newpage(); grid::grid.draw(ph$gtable); dev.off()
  invisible(fb)
}

# ---- QC --------------------------------------------------------------------
plot_pca <- function(mat, group, cols, title = "Proteome PCA") {
  # impute column means only for the PCA projection (does not touch DE input)
  m <- mat; for (j in seq_len(ncol(m))) m[is.na(m[, j]), j] <- mean(m[, j], na.rm = TRUE)
  m <- m[matrixStats::rowVars(m) > 0, , drop = FALSE]
  pc <- prcomp(t(m), scale. = TRUE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  d <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                  sample = colnames(mat), group = group)
  ggplot(d, aes(PC1, PC2, colour = group, label = sample)) +
    geom_point(size = 2.4) +
    ggrepel::geom_text_repel(size = 2.3, show.legend = FALSE, max.overlaps = 20) +
    scale_colour_manual(values = cols) +
    labs(title = title, x = sprintf("PC1 (%.1f%%)", ve[1]),
         y = sprintf("PC2 (%.1f%%)", ve[2])) + theme_mo()
}

plot_corr_heatmap <- function(mat, group, cols, title = "Sample correlation") {
  cm <- cor(mat, use = "pairwise.complete.obs")
  ann <- data.frame(group = group); rownames(ann) <- colnames(mat)
  pheatmap::pheatmap(cm, main = title, annotation_col = ann,
                     annotation_colors = list(group = cols[names(cols) %in% group]),
                     show_colnames = FALSE, silent = TRUE, fontsize = 7,
                     color = grDevices::colorRampPalette(c("#4575b4", "white", "#d73027"))(50))
}

# ---- DE --------------------------------------------------------------------
# Native ggplot2 volcano: faded NS background, popped significant points,
# y-axis trimmed to the data, tidy repelled labels (top hits only).
plot_volcano <- function(d, contrast, padj_cut, cols = c(Up="#d73027",NS="#bdbdbd",Down="#4575b4"),
                         n_label = 10L, label_size = 2.6) {
  if (!.have("ggrepel")) stop("ggrepel is required for volcano labels")
  d <- as.data.table(copy(d))[is.finite(log2FC) & is.finite(padj)]
  d[, padj_plot := pmax(padj, .Machine$double.xmin)]
  d[, neglog10 := -log10(padj_plot)]
  d[, label_plot := label]
  d[is.na(label_plot) | !nzchar(label_plot), label_plot := protein_id]
  dl <- dir_labels()                                    # Higher / Lower in <test>
  d[, status := factor(
      fifelse(regulation == "Up", dl[["Up"]],
      fifelse(regulation == "Down", dl[["Down"]], "NS")),
      levels = c(dl[["Up"]], dl[["Down"]], "NS"))]
  # draw NS first (background), significant points on top
  setorder(d, -status)

  status_cols <- setNames(unname(cols[c("Up", "Down", "NS")]), c(dl[["Up"]], dl[["Down"]], "NS"))

  x_lim <- ceiling(max(abs(d$log2FC)) * 2) / 2          # symmetric, snapped to 0.5
  y_top <- max(d$neglog10) * 1.12                       # trim empty upper space
  n_up  <- sum(d$regulation == "Up")
  n_down<- sum(d$regulation == "Down")
  top   <- d[regulation != "NS"][order(padj, -abs(log2FC))][seq_len(min(n_label, .N))]

  ggplot(d, aes(log2FC, neglog10)) +
    geom_vline(xintercept = 0, linetype = "22", linewidth = 0.3, colour = "grey60") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "22",
               linewidth = 0.3, colour = "grey60") +
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
         caption = sprintf("%d %s | %d %s", n_up, dir_labels(lower = TRUE)[["Up"]], n_down, dir_labels(lower = TRUE)[["Down"]]),
         x = expression(Log[2]~"fold change"),
         y = expression(-Log[10]~"adjusted p-value")) +
    theme_mo(9) +
    theme(legend.position = "right",
          plot.subtitle = element_text(colour = "grey30"),
          plot.caption  = element_text(colour = "grey35"),
          plot.margin   = margin(6, 10, 6, 6))
}

plot_ma <- function(d, contrast, padj_cut, cols = c(Up="#d73027",NS="#bdbdbd",Down="#4575b4")) {
  d <- as.data.table(d); d <- d[is.finite(log2FC) & is.finite(AveExpr)]
  dl <- dir_labels()
  ggplot(d, aes(AveExpr, log2FC, colour = regulation)) +
    geom_point(size = 1, alpha = 0.7) +
    scale_colour_manual(values = cols, name = NULL,
                        labels = c(Up = dl[["Up"]], Down = dl[["Down"]], NS = "NS")) +
    geom_hline(yintercept = 0, linetype = 3, colour = "grey60") +
    labs(title = sprintf("MA - %s", contrast),
         x = "mean log2 abundance", y = "log2 fold change") + theme_mo()
}

plot_pval_hist <- function(d, contrast) {
  d <- as.data.table(d)[is.finite(pvalue)]
  ggplot(d, aes(pvalue)) +
    geom_histogram(breaks = seq(0, 1, 0.05), fill = "#4575b4", colour = "white") +
    labs(title = sprintf("p-value histogram - %s", contrast),
         x = "raw p-value", y = "proteins") + theme_mo()
}

plot_count_bar <- function(summ_dt, title) {
  m <- melt(summ_dt[, .(contrast, Up = up, Down = down)],
            id.vars = "contrast", variable.name = "regulation", value.name = "n")
  ggplot(m, aes(contrast, n, fill = regulation)) +
    geom_col(position = position_dodge(0.7), width = 0.65) +
    geom_text(aes(label = n), position = position_dodge(0.7), vjust = -0.3, size = 2.4) +
    scale_fill_manual(values = c(Up = "#d73027", Down = "#4575b4"), name = NULL,
                      labels = c(Up = dir_labels()[["Up"]], Down = dir_labels()[["Down"]])) +
    labs(title = title, x = NULL, y = "significant proteins") +
    theme_mo() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

plot_de_heatmap <- function(zmat, row_labels, group, cols, title, show_rownames = TRUE) {
  ann <- data.frame(group = group); rownames(ann) <- colnames(zmat)
  rn <- if (show_rownames) row_labels else rep("", nrow(zmat))
  rownames(zmat) <- make.unique(ifelse(is.na(rn) | rn == "", "", rn))
  pheatmap::pheatmap(zmat, main = title, annotation_col = ann,
                     annotation_colors = list(group = cols[names(cols) %in% group]),
                     show_rownames = show_rownames, cluster_cols = TRUE,
                     silent = TRUE, fontsize = 7, fontsize_row = 5,
                     na_col = "grey85",
                     color = grDevices::colorRampPalette(c("#4575b4","white","#d73027"))(50))
}

# z-score rows NA-robustly (for heatmaps on un-imputed matrices)
row_zscore <- function(m) t(apply(m, 1, function(x) {
  mu <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) x - mu else (x - mu) / s
}))

plot_overlap <- function(sets, title) {
  sets <- sets[vapply(sets, length, 0L) > 0]
  if (length(sets) < 2) return(invisible(NULL))
  UpSetR::upset(UpSetR::fromList(sets), order.by = "freq",
                mainbar.y.label = "shared DE proteins", sets.x.label = "DE per contrast")
}
