#!/usr/bin/env Rscript
# =============================================================================
# 09b_panel_variants.R
# Publication composite panels for Lactate vs Saline, built natively (patchwork
# + grid grobs) from the cached DE / matrix / enrichment tables; nothing is
# re-run.
#
# Sub-panels
#   a  Sample PCA                  (imputed matrix)               [cf. prot__qc__all__pca_limma]
#   b  Volcano                     (limma DE)                     [cf. prot__de__..__volcano_limma]
#   c  Detected-proteome Venn      (VennDiagram, scaled = FALSE)  [cf. prot__detected__..__venn]
#   d  Identified-proteome heatmap (all quantified proteins)      [cf. prot__de__..__heatmap_all_proteins]
#   e  GO GSEA bubble              (BP / MF / CC, NES)            [cf. prot__enrich__..__go_bubble_gsea]
#   f  COG sub-panel, two versions:
#        v1  log2FC strip          (beeswarm per category)        [cf. prot__enrich__..__cog_lfc_strip]
#        v2  over-representation   (bubble, category x direction) [cf. prot__enrich__..__cog_bubble_ora]
#   v3 = same layout without f.
#
# Inputs : results/ tables from 01, 05, 05c and 08
# Outputs: figures/panel/prot__panel__<contrast>__v1_cog_lfc_strip.{svg,png}
#          figures/panel/prot__panel__<contrast>__v2_cog_bubble_ora.{svg,png}
#          figures/panel/prot__panel__<contrast>__v3_no_cog.{svg,png}
#          PNG at 600 dpi; SVG fully vector (vector colourbars, heatmap with
#          use_raster = FALSE), so nothing turns into a black <image> box in
#          Affinity / Illustrator.
# Run    : Rscript scripts/09b_panel_variants.R   (after 01, 05, 05c, 08)
#
# Steps
#   0. Setup
#   1. Shared inputs
#   2. Panel a: PCA
#   3. Panel b: volcano
#   4. Panel c: detected-proteome Venn
#   5. Panel d: identified-proteome heatmap
#   6. Panel e: GO GSEA bubble
#   7. Panel f: COG strip (v1) and COG ORA bubble (v2)
#   8. Assemble and save
# =============================================================================


# 0. Setup ---------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(ComplexHeatmap); library(circlize); library(VennDiagram); library(grid)
  library(ggbeeswarm)
})
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
SDIR <- if (length(script_file)) dirname(normalizePath(script_file)) else "scripts"
source(file.path(SDIR, "config.R"))
source(file.path(SDIR, "io.R"))
source(file.path(SDIR, "plot_utils.R"))
source(file.path(SDIR, "enrich_utils.R"))
invisible(futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))

cn     <- CFG$contrasts[[1]]$name
tg     <- CFG$contrasts[[1]]$test
cg     <- CFG$contrasts[[1]]$control
cn_lab <- gsub("_", " ", cn, fixed = TRUE)
PADJ   <- CFG$padj
dl     <- dir_labels(tg)                                  # Higher / Lower in <test>
dl_low <- dir_labels(tg, lower = TRUE)
panel_dir <- file.path(CFG$dir_figures, "panel")
if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)

# Small-text title / subtitle used by the ggplot sub-panels.
panel_text_theme <- theme(plot.title = element_text(face = "bold"),
                          plot.subtitle = element_text(colour = "grey35", size = 7))


# 1. Shared inputs -------------------------------------------------------------
imp       <- read_imputed_matrix(CFG)                      # columns: Saline then Lactate
abundance <- imp$abundance
condition <- imp$samples$group
de <- read_de_results(CFG, cn, na.strings = c("", "NA"))


# 2. Panel a: PCA --------------------------------------------------------------
# Same PCA as plot_pca(); legend below and labels allowed into the margin
# (clip = off, extra x/y expansion) so long sample names fit the narrow cell.
pca <- pca_scores(abundance)
pca_d <- pca$scores
pca_d$group <- as.character(condition)
p_pca <- ggplot(pca_d, aes(PC1, PC2, colour = group, label = sample)) +
  geom_point(size = 2.8) +
  ggrepel::geom_text_repel(size = 2.3, show.legend = FALSE, max.overlaps = 20,
                           box.padding = 0.5, min.segment.length = 0,
                           segment.size = 0.25, segment.colour = "grey65", seed = 42L) +
  scale_colour_manual(values = GROUP_COLORS) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.08)) +
  coord_cartesian(clip = "off") +
  labs(title = "Sample PCA", x = sprintf("PC1 (%.1f%%)", pca$var_pct[1]),
       y = sprintf("PC2 (%.1f%%)", pca$var_pct[2])) +
  theme_mo(10) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(face = "bold"), axis.text = element_text(size = 8),
        plot.margin = margin(6, 12, 4, 12))


# 3. Panel b: volcano ----------------------------------------------------------
# Fewer and smaller labels than the stand-alone volcano.
p_volc <- plot_volcano(de, cn, PADJ, REG_COLORS, n_label = 6L, label_size = 2.3) +
  labs(title = "Differential abundance", subtitle = NULL) +
  theme(plot.title = element_text(face = "bold"))


# 4. Panel c: detected-proteome Venn -------------------------------------------
# Same geometry as 08's Venn, smaller text, drawn into a grob.
venn_panel_grob <- function(sets, labels, colors, title) {
  colors <- unname(colors)
  vg <- venn.diagram(
    x = sets, filename = NULL, category.names = c("", ""),
    scaled = FALSE, ext.text = FALSE,
    fill = colors, alpha = 0.45, col = colors, lwd = 2.2, lty = "solid",
    cex = 1.25, fontface = "bold", fontfamily = "sans", fontcolour = "grey15",
    cat.cex = 0.01, cat.col = "white", margin = 0.08,
    print.mode = "raw", disable.logging = TRUE)
  # VennDiagram draws the larger set on the left regardless of list order.
  left_idx  <- which.max(lengths(sets))
  right_idx <- setdiff(seq_along(sets), left_idx)
  set_label <- function(idx, x)
    grid.text(sprintf("%s (n=%d)", labels[idx], length(sets[[idx]])), x = x, y = 0.83,
              gp = gpar(fontface = "bold", fontsize = 9.5, col = colors[idx], fontfamily = "sans"))
  grid.grabExpr({
    grid.rect(gp = gpar(fill = "white", col = NA))
    grid.text(title, y = 0.955, gp = gpar(fontface = "bold", fontfamily = "sans", fontsize = 11))
    set_label(left_idx, 0.27)
    set_label(right_idx, 0.73)
    pushViewport(viewport(y = 0.40, height = 0.74, width = 0.96))
    grid.draw(vg)
    upViewport()
  })
}

memb <- fread(file.path(CFG$dir_results, "prot__detected__Lactate_vs_Saline__membership.tsv"), sep = "\t")
det_sets <- list(Lactate = memb[detected_Lactate == TRUE, protein_id],
                 Saline  = memb[detected_Saline  == TRUE, protein_id])
g_venn <- venn_panel_grob(det_sets, labels = c("Lactate", "Saline"),
                          colors = c(GROUP_COLORS[["Lactate"]], GROUP_COLORS[["Saline"]]),
                          title = "Detected proteome overlap")


# 5. Panel d: identified-proteome heatmap --------------------------------------
# Same content as 08's all-protein heatmap, smaller text and bars.
NS_LAB <- "Not significant"
status_cols_hm <- c(setNames(c(REG_COLORS[["Up"]], REG_COLORS[["Down"]]), unname(dl)),
                    setNames(REG_COLORS[["NS"]], NS_LAB))
small_legend <- list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7))

all_de   <- status_ordered_proteins(de, rownames(abundance), NS_LAB)
zall     <- row_zscore(abundance[all_de$protein_id, , drop = FALSE], clip = 2)
n_status <- table(all_de$status)

top_ha <- HeatmapAnnotation(
  Condition = condition,
  col = list(Condition = GROUP_COLORS[CFG$group_levels]),
  annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(Condition = small_legend))
ht_all <- Heatmap(
  zall, name = "row z-score",
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  top_annotation = top_ha,
  left_annotation = rowAnnotation(
    Direction = all_de$status, col = list(Direction = status_cols_hm),
    show_annotation_name = FALSE, annotation_legend_param = list(Direction = small_legend)),
  right_annotation = rowAnnotation(
    log2FC = anno_barplot(all_de$log2FC, baseline = 0,
                          gp = gpar(fill = unname(status_cols_hm[as.character(all_de$status)]), col = NA),
                          width = unit(16, "mm")),
    annotation_name_rot = 0, annotation_name_gp = gpar(fontsize = 8)),
  row_split = all_de$status, cluster_rows = TRUE, cluster_row_slices = FALSE,
  row_title = sprintf("%s\n(n=%d)", names(n_status), as.integer(n_status)),
  row_title_rot = 0, row_title_gp = gpar(fontsize = 8, fontface = "bold"),
  row_gap = unit(1.5, "mm"),
  cluster_columns = TRUE, cluster_column_slices = FALSE, column_split = condition,
  show_row_names = FALSE, column_names_rot = 45, column_names_gp = gpar(fontsize = 8),
  column_title = sprintf("Identified proteome (n=%d quantified)", nrow(all_de)),
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = c(list(at = c(-2, -1, 0, 1, 2)), small_legend),
  use_raster = FALSE)
g_heat <- grid.grabExpr(draw(ht_all, merge_legends = TRUE, heatmap_legend_side = "right",
                             annotation_legend_side = "right"))


# 6. Panel e: GO GSEA bubble ---------------------------------------------------
gse <- go_gsea_top(cn, 8L, CFG)                           # top 8 sets per aspect
p_go <- ggplot(gse, aes(NES, ykey, fill = neglogp, size = setSize, colour = sig)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  bubble_layers("set size", size_range = c(1.8, 6.5), key_size = 3) +
  facet_grid(aspect ~ ., scales = "free_y", space = "free_y") +
  labs(title = "GO GSEA",
       subtitle = sprintf("top 8 sets/aspect; NES>0 = %s, NES<0 = %s", dl_low[["Up"]], dl_low[["Down"]]),
       x = "normalised enrichment score (NES)", y = NULL) +
  theme_mo(8) +
  theme(panel.grid.major.y = element_line(colour = "grey92"),
        strip.text = element_text(face = "bold", size = 7.5),
        axis.text.y = element_text(size = 6.8),
        legend.key.size = unit(3, "mm")) +
  panel_text_theme


# 7. Panel f: COG strip (v1) and COG ORA bubble (v2) ---------------------------
## 7.1 Per-protein COG letters (written by 05c) ----
pmap <- read_result_table(sprintf("prot__enrich__%s__cog_protein_map.tsv", cn), CFG)
pmap <- pmap[!is.na(cog) & nzchar(cog)]
annot_n <- uniqueN(pmap$protein_id)                       # detected proteins with a COG letter

## 7.2 v1: log2FC strip, categories with >= minGSSize detected proteins ----
cexp <- pmap[, .(letter = strsplit(cog, "", fixed = TRUE)[[1]]),
             by = .(protein_id, log2FC, padj, regulation)]
present <- cexp[, .(n = uniqueN(protein_id)), by = letter][n >= CFG$gsea_minGS, letter]
p_cog_strip <- plot_cog_lfc_strip(
  cexp[letter %chin% present], tg, cg,
  title = "log2FC by COG category",
  subtitle = sprintf("each point = one detected protein; | = category median; >0 = %s, <0 = %s",
                     dl_low[["Up"]], dl_low[["Down"]])) +
  panel_text_theme

## 7.3 v2: ORA bubble, category x direction ----
dir_lab <- setNames(unname(dl), c("up", "down"))
ob <- rbindlist(lapply(names(dir_lab), function(dir) {
  x <- read_result_table(sprintf("prot__enrich__%s_%s__cog_ora_all.tsv", cn, dir), CFG)
  if (is.null(x) || !nrow(x)) return(NULL)
  x[, .(letter = ID, FoldEnrichment, Count, pvalue, p.adjust, direction = dir_lab[[dir]])]
}), fill = TRUE)
p_cog_ora <- plot_cog_ora_bubble(
  ob, tg, PADJ,
  title = "COG over-representation",
  subtitle = sprintf("universe = detected proteome with a COG letter (n=%d); outline = BH FDR<=0.05", annot_n),
  size_range = c(2, 7)) +
  panel_text_theme


# 8. Assemble and save ---------------------------------------------------------
## 8.1 Layouts (15-column grid) ----
# With f: a | b | c top band (rows 1-4); d heatmap as a tall left column
# (rows 5-15, cols 1-7); e GO GSEA (top) and f COG (bottom) stacked on the
# right (cols 8-15).
LAYOUT_WITH_COG <- c(
  area(1,  1,  4,  7),   # a PCA
  area(1,  8,  4, 12),   # b volcano
  area(1, 13,  4, 15),   # c Venn
  area(5,  1, 15,  7),   # d heatmap
  area(5,  8,  9, 15),   # e GO GSEA
  area(10, 8, 15, 15))   # f COG
# Without f: 12-row grid on a shorter canvas, so the a | b | c band keeps its
# height (~120 mm) and d + e share the bottom rows side by side.
LAYOUT_NO_COG <- c(
  area(1,  1,  4,  7),   # a PCA
  area(1,  8,  4, 12),   # b volcano
  area(1, 13,  4, 15),   # c Venn
  area(5,  1, 12,  7),   # d heatmap
  area(5,  8, 12, 15))   # e GO GSEA

## 8.2 Builders ----
assemble_panel <- function(plots, layout) {
  wrap_plots(plots, design = layout) +
    plot_annotation(
      tag_levels = "a",
      title = bquote(bolditalic("Aspergillus fumigatus") ~
                     bold("proteome response to lactate (" * .(cn_lab) * ")")),
      theme = theme(plot.title = element_text(face = "bold", size = 15, margin = margin(b = 4)))) &
    theme(plot.tag = element_text(face = "bold", size = 15))
}

PANEL_WIDTH_MM <- 380
save_panel <- function(pnl, suffix, height_mm) {
  base <- file.path(panel_dir, sprintf("prot__panel__%s__%s", cn, suffix))
  svglite::svglite(paste0(base, ".svg"), width = mm2in(PANEL_WIDTH_MM), height = mm2in(height_mm))
  print(pnl); invisible(dev.off())
  ragg::agg_png(paste0(base, ".png"), width = PANEL_WIDTH_MM, height = height_mm, units = "mm", res = 600)
  print(pnl); invisible(dev.off())
  message("DONE: ", base, ".{svg,png}")
}

## 8.3 Render the three versions ----
shared_plots <- list(p_pca, p_volc, wrap_elements(full = g_venn), wrap_elements(full = g_heat), p_go)
save_panel(assemble_panel(c(shared_plots, list(p_cog_strip)), LAYOUT_WITH_COG), "v1_cog_lfc_strip", 440)
save_panel(assemble_panel(c(shared_plots, list(p_cog_ora)),   LAYOUT_WITH_COG), "v2_cog_bubble_ora", 440)
save_panel(assemble_panel(shared_plots,                       LAYOUT_NO_COG),   "v3_no_cog", 360)

message(sprintf(paste0("  a PCA | b volcano | c detected-venn(Lac=%d,Sal=%d) | d heatmap(n=%d quantified) | ",
                       "e GO GSEA(%d sets) | f = v1 COG strip(%d cats) / v2 COG ORA(%d rows)"),
                length(det_sets$Lactate), length(det_sets$Saline), nrow(all_de),
                nrow(gse), length(present), nrow(ob)))
