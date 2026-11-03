#!/usr/bin/env Rscript
# ============================================================================
# 09b_panel_variants.R  -  Publication-ready composite panels for Lactate vs
# Saline, assembled from a user-specified figure set. Native/vector build
# (patchwork + grid grobs, reads cached DE / matrix / enrichment TSVs, NO re-run),
# in TWO versions that differ only in the COG sub-panel (F).
#
# Shared sub-panels (both versions):
#   A  Sample PCA                 (plot_pca, imputed matrix)       [prot__qc__all__pca_limma]
#   B  Volcano                    (plot_volcano, limma DE)         [prot__de__..__volcano_limma]
#   C  Detected proteome Venn     (VennDiagram, scaled = FALSE)    [prot__detected__..__venn]
#   D  Identified-proteome heatmap(ComplexHeatmap, all quantified) [prot__de__..__heatmap_all_proteins]
#   E  GO GSEA bubble             (BP/MF/CC, NES)                  [prot__enrich__..__go_bubble_gsea]
#
# Swapped sub-panel F:
#   v1  COG log2FC strip          (ggbeeswarm, per category)       [prot__enrich__..__cog_lfc_strip]
#   v2  COG over-representation   (bubble, category x direction)   [prot__enrich__..__cog_bubble_ora]
#
# Layout (portrait): A|B|C short top band; D|E medium band; F full-width tall band.
# Output:
#   figures/panel/prot__panel__Lactate_vs_Saline__v1_cog_lfc_strip.{svg,png}
#   figures/panel/prot__panel__Lactate_vs_Saline__v2_cog_bubble_ora.{svg,png}
# (PNG @ 600 dpi; SVG fully vector - colourbars drawn with guide_colourbar(raster
#  = FALSE) and the heatmap kept use_raster = FALSE, so nothing becomes a black
#  <image> box in Affinity/Illustrator.)
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(ComplexHeatmap); library(circlize); library(VennDiagram); library(grid)
  library(ggbeeswarm)
})

args <- commandArgs(FALSE)
self <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
SDIR <- if (length(self)) dirname(normalizePath(self)) else "scripts"
source(file.path(SDIR, "config.R"))
source(file.path(SDIR, "io.R"))
source(file.path(SDIR, "plot_utils.R"))
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

cn     <- CFG$contrasts[[1]]$name
tg     <- CFG$contrasts[[1]]$test
cg     <- CFG$contrasts[[1]]$control
cn_lab <- gsub("_", " ", cn, fixed = TRUE)
PADJ   <- CFG$padj
RES    <- CFG$dir_results
dl     <- dir_labels(tg)                                 # Higher / Lower in <test>
panel_dir <- file.path(CFG$dir_figures, "panel")
if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)
rd <- function(f) { p <- file.path(RES, f); if (file.exists(p)) fread(p) else NULL }

# ---------------------------------------------------------------------------
# Shared inputs (imputed matrix + limma DE), exactly as 08/09 read them
# ---------------------------------------------------------------------------
samples <- unique(fread(CFG$samples_tsv, sep = "\t")[, .(sample_id, group)])
samples[, group := factor(group, levels = CFG$group_levels)]
setorder(samples, group, sample_id)

mdt <- fread(file.path(RES, "prot__imputed__all__matrix_limma.tsv"),
             sep = "\t", na.strings = c("", "NA"))
abundance <- as.matrix(mdt[, setdiff(names(mdt), "protein_id"), with = FALSE])
storage.mode(abundance) <- "double"
rownames(abundance) <- mdt$protein_id
sample_order <- samples$sample_id[samples$sample_id %in% colnames(abundance)]
abundance <- abundance[, sample_order, drop = FALSE]
grp_vec <- samples$group[match(colnames(abundance), samples$sample_id)]

de <- fread(file.path(RES, sprintf("prot__de__%s__results_limma.tsv", cn)),
            sep = "\t", na.strings = c("", "NA"))

# ===========================================================================
# A. PCA  (prot__qc__all__pca_limma)
# Panel build of plot_pca (same prcomp: col-mean impute, drop zero-var, scale):
# legend moved below and labels drawn into the margin (clip = off) + x/y
# expansion, so the long sample names fit in the narrow top-band cell.
# ===========================================================================
.pca <- {
  m <- abundance
  for (j in seq_len(ncol(m))) m[is.na(m[, j]), j] <- mean(m[, j], na.rm = TRUE)
  m <- m[matrixStats::rowVars(m) > 0, , drop = FALSE]
  pc <- prcomp(t(m), scale. = TRUE); ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  list(d = data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                      sample = colnames(abundance), group = as.character(grp_vec)), ve = ve)
}
p_pca <- ggplot(.pca$d, aes(PC1, PC2, colour = group, label = sample)) +
  geom_point(size = 2.8) +
  ggrepel::geom_text_repel(size = 2.3, show.legend = FALSE, max.overlaps = 20,
                           box.padding = 0.5, min.segment.length = 0,
                           segment.size = 0.25, segment.colour = "grey65", seed = 42L) +
  scale_colour_manual(values = GROUP_COLORS) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.08)) +
  coord_cartesian(clip = "off") +
  labs(title = "Sample PCA", x = sprintf("PC1 (%.1f%%)", .pca$ve[1]),
       y = sprintf("PC2 (%.1f%%)", .pca$ve[2])) +
  theme_mo(10) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(face = "bold"), axis.text = element_text(size = 8),
        plot.margin = margin(6, 12, 4, 12))

# ===========================================================================
# B. Volcano  (prot__de__..__volcano_limma) - fewer/smaller labels for the panel
# ===========================================================================
p_volc <- plot_volcano(de, cn, PADJ, REG_COLORS, n_label = 6L, label_size = 2.3) +
  labs(title = "Differential abundance", subtitle = NULL) +
  theme(plot.title = element_text(face = "bold"))

# ===========================================================================
# C. Detected-proteome overlap Venn  (prot__detected__..__venn) as a grob
# ===========================================================================
memb <- fread(file.path(RES, "prot__detected__Lactate_vs_Saline__membership.tsv"), sep = "\t")
det_sets <- list(Lactate = memb[detected_Lactate == TRUE, protein_id],
                 Saline  = memb[detected_Saline  == TRUE, protein_id])

venn_panel_grob <- function(sets, labels, colors, title) {
  colors <- unname(colors)
  vg <- venn.diagram(
    x = sets, filename = NULL, category.names = c("", ""),
    scaled = FALSE, ext.text = FALSE,
    fill = colors, alpha = 0.45, col = colors, lwd = 2.2, lty = "solid",
    cex = 1.25, fontface = "bold", fontfamily = "sans", fontcolour = "grey15",
    cat.cex = 0.01, cat.col = "white", margin = 0.08,
    print.mode = "raw", disable.logging = TRUE)
  grid.grabExpr({
    grid.rect(gp = gpar(fill = "white", col = NA))
    grid.text(title, y = 0.955,
              gp = gpar(fontface = "bold", fontfamily = "sans", fontsize = 11))
    left_idx  <- which.max(lengths(sets))
    right_idx <- setdiff(seq_along(sets), left_idx)
    grid.text(sprintf("%s (n=%d)", labels[left_idx], length(sets[[left_idx]])),
              x = 0.27, y = 0.83,
              gp = gpar(fontface = "bold", fontsize = 9.5, col = colors[left_idx],
                        fontfamily = "sans"))
    grid.text(sprintf("%s (n=%d)", labels[right_idx], length(sets[[right_idx]])),
              x = 0.73, y = 0.83,
              gp = gpar(fontface = "bold", fontsize = 9.5, col = colors[right_idx],
                        fontfamily = "sans"))
    pushViewport(viewport(y = 0.40, height = 0.74, width = 0.96))
    grid.draw(vg); upViewport()
  })
}
g_venn <- venn_panel_grob(
  det_sets, labels = c("Lactate", "Saline"),
  colors = c(GROUP_COLORS[["Lactate"]], GROUP_COLORS[["Saline"]]),
  title = "Detected proteome overlap")

# ===========================================================================
# D. Identified-proteome heatmap of ALL quantified proteins
#    (prot__de__..__heatmap_all_proteins) - ComplexHeatmap grob, same aesthetics
#    as 08 (row names hidden, split Higher/Lower in <test>/Not significant).
# ===========================================================================
zscore_clip <- function(hm) {
  row_mean <- rowMeans(hm, na.rm = TRUE); row_sd <- apply(hm, 1, sd, na.rm = TRUE)
  z <- sweep(hm, 1, row_mean, "-")
  ok <- is.finite(row_sd) & row_sd > 0
  z[ok, ] <- sweep(z[ok, , drop = FALSE], 1, row_sd[ok], "/"); z[!ok, ] <- 0
  z[z < -2] <- -2; z[z > 2] <- 2; z
}
cond <- factor(grp_vec, levels = CFG$group_levels)
heat_col <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
dir_cols <- setNames(c(REG_COLORS[["Up"]], REG_COLORS[["Down"]]), unname(dl))
NS_LAB   <- "Not significant"
status_cols_hm <- c(dir_cols, setNames(REG_COLORS[["NS"]], NS_LAB))

all_de <- de[protein_id %chin% rownames(abundance)]
all_de[, status := factor(
  fifelse(!is.na(padj) & padj <= PADJ & regulation == "Up", dl[["Up"]],
  fifelse(!is.na(padj) & padj <= PADJ & regulation == "Down", dl[["Down"]], NS_LAB)),
  levels = c(unname(dl), NS_LAB))]
setorder(all_de, status, padj)
zall <- zscore_clip(abundance[all_de$protein_id, sample_order, drop = FALSE])
n_status <- table(all_de$status)

top_ha <- HeatmapAnnotation(Condition = cond,
    col = list(Condition = GROUP_COLORS[CFG$group_levels]),
    annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 8),
    annotation_legend_param = list(Condition = list(title_gp = gpar(fontsize = 8),
                                                    labels_gp = gpar(fontsize = 7))))
ht_all <- Heatmap(
  zall, name = "row z-score", col = heat_col, top_annotation = top_ha,
  left_annotation = rowAnnotation(Direction = all_de$status,
      col = list(Direction = status_cols_hm), show_annotation_name = FALSE,
      annotation_legend_param = list(Direction = list(title_gp = gpar(fontsize = 8),
                                                      labels_gp = gpar(fontsize = 7)))),
  right_annotation = rowAnnotation(log2FC = anno_barplot(all_de$log2FC, baseline = 0,
      gp = gpar(fill = unname(status_cols_hm[as.character(all_de$status)]), col = NA),
      width = unit(16, "mm")), annotation_name_rot = 0,
      annotation_name_gp = gpar(fontsize = 8)),
  row_split = all_de$status, cluster_rows = TRUE, cluster_row_slices = FALSE,
  row_title = sprintf("%s\n(n=%d)", names(n_status), as.integer(n_status)),
  row_title_rot = 0, row_title_gp = gpar(fontsize = 8, fontface = "bold"),
  row_gap = unit(1.5, "mm"),
  cluster_columns = TRUE, cluster_column_slices = FALSE, column_split = cond,
  show_row_names = FALSE, column_names_rot = 45, column_names_gp = gpar(fontsize = 8),
  column_title = sprintf("Identified proteome (n=%d quantified)", nrow(all_de)),
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(at = c(-2, -1, 0, 1, 2), title_gp = gpar(fontsize = 8),
                              labels_gp = gpar(fontsize = 7)),
  use_raster = FALSE)
g_heat <- grid.grabExpr(draw(ht_all, merge_legends = TRUE,
                             heatmap_legend_side = "right",
                             annotation_legend_side = "right"))

# ===========================================================================
# E. GO GSEA bubble  (prot__enrich__..__go_bubble_gsea) - rebuilt from tables
# ===========================================================================
ASP_LEVELS <- c(BP = "Biological process", MF = "Molecular function", CC = "Cellular component")
tidy_label_go <- function(x, max = 80, width = 42) {
  x <- ifelse(nchar(x) > max, paste0(substr(x, 1, max - 3), "..."), x)
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
}
gse <- rbindlist(lapply(names(ASP_LEVELS), function(asp) {
  x <- rd(sprintf("prot__enrich__%s__go_%s_gsea_all.tsv", cn, tolower(asp)))
  if (is.null(x) || !nrow(x)) return(NULL)
  x[, aspect := ASP_LEVELS[[asp]]]
  x[order(pvalue)][seq_len(min(8L, .N))]                 # top 8 per aspect
}), fill = TRUE)
gse[, aspect := factor(aspect, levels = ASP_LEVELS)]
gse[, sig := p.adjust <= PADJ]
gse[, neglogp := -log10(pvalue)]
gse[, Description := tidy_label_go(Description)]
setorderv(gse, c("aspect", "NES"))
gse[, ykey := factor(.I, levels = .I, labels = Description)]

p_go <- ggplot(gse, aes(NES, ykey, fill = neglogp, size = setSize, colour = sig)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_point(shape = 21, stroke = 0.5, alpha = 0.9) +
  scale_fill_gradient(low = "grey85", high = "#d73027", name = expression(-log[10]~italic(p))) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey70"),
                      labels = c(`TRUE` = "FDR<=0.05", `FALSE` = "n.s."), name = "BH") +
  scale_size_continuous(range = c(1.8, 6.5), name = "set size") +
  guides(colour = guide_legend(override.aes = list(size = 3)),
         fill = guide_colourbar(raster = FALSE)) +   # vector colourbar (Affinity-safe)
  facet_grid(aspect ~ ., scales = "free_y", space = "free_y") +
  labs(title = "GO GSEA",
       subtitle = sprintf("top 8 sets/aspect; NES>0 = %s, NES<0 = %s",
                          dir_labels(tg, lower = TRUE)[["Up"]], dir_labels(tg, lower = TRUE)[["Down"]]),
       x = "normalised enrichment score (NES)", y = NULL) +
  theme_mo(8) +
  theme(panel.grid.major.y = element_line(colour = "grey92"),
        strip.text = element_text(face = "bold", size = 7.5),
        axis.text.y = element_text(size = 6.8),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 7),
        legend.key.size = unit(3, "mm"))

# ===========================================================================
# COG scaffold (letters, names, supergroups) - shared by F v1 / v2
# ===========================================================================
GROUPS <- c("Information storage\n& processing", "Cellular processes\n& signaling",
            "Metabolism", "Poorly\ncharacterized")
COG <- data.table(
  letter = c("J", "A", "K", "L", "B",
             "D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O",
             "C", "G", "E", "F", "H", "I", "P", "Q",
             "R", "S"),
  name = c("Translation, ribosomal structure & biogenesis", "RNA processing & modification",
           "Transcription", "Replication, recombination & repair",
           "Chromatin structure & dynamics",
           "Cell cycle control & division", "Nuclear structure", "Defense mechanisms",
           "Signal transduction", "Cell wall/membrane/envelope biogenesis", "Cell motility",
           "Cytoskeleton", "Extracellular structures", "Intracellular trafficking & secretion",
           "PTM, protein turnover, chaperones",
           "Energy production & conversion", "Carbohydrate transport & metabolism",
           "Amino acid transport & metabolism", "Nucleotide transport & metabolism",
           "Coenzyme transport & metabolism", "Lipid transport & metabolism",
           "Inorganic ion transport & metabolism", "Secondary metabolite biosynthesis",
           "General function prediction only", "Function unknown"),
  group = rep(GROUPS, c(5, 10, 8, 2)))
COG[, group := factor(group, levels = GROUPS)]
COG[, lab := sprintf("%s  %s", letter, name)]
cog_axis <- function(dt) {
  dt <- merge(dt, COG[, .(letter, lab, group)], by = "letter")
  dt[, lab := factor(lab, levels = rev(COG$lab))]; dt
}
cog_theme <- function() list(theme_mo(8),
  theme(plot.title.position = "plot",
        panel.grid.major.y = element_line(colour = "grey92"),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 7),
        strip.text = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0, hjust = 0)))

# per-protein COG map (written by 05c): protein_id, log2FC, padj, cog(letters)
pmap <- rd(sprintf("prot__enrich__%s__cog_protein_map.tsv", cn))
pmap <- pmap[!is.na(cog) & nzchar(cog)]
annot_n <- uniqueN(pmap$protein_id)                      # detected proteins with a COG letter

# ---- F v1: log2FC strip (beeswarm), one point per detected protein ----------
cexp <- pmap[, .(letter = strsplit(cog, "", fixed = TRUE)[[1]]),
             by = .(protein_id, log2FC, padj)]
present <- cexp[, .(n = uniqueN(protein_id)), by = letter][n >= CFG$gsea_minGS, letter]
st <- cexp[letter %chin% present]
st[, status := fifelse(!is.na(padj) & padj <= PADJ & log2FC > 0, dl[["Up"]],
               fifelse(!is.na(padj) & padj <= PADJ & log2FC < 0, dl[["Down"]], "n.s."))]
st[, status := factor(status, levels = c(unname(dl), "n.s."))]
st <- cog_axis(st)
nlab <- st[, .N, by = lab]
st[, lab := factor(sprintf("%s (n=%d)", lab, nlab$N[match(lab, nlab$lab)]),
                   levels = sprintf("%s (n=%d)", levels(lab), nlab$N[match(levels(lab), nlab$lab)]))]
med <- st[, .(med = median(log2FC)), by = .(lab, group)]
st_cols <- setNames(c(REG_COLORS[["Up"]], REG_COLORS[["Down"]], "#bdbdbd"), levels(st$status))
p_cog_strip <- ggplot(st, aes(log2FC, lab)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_quasirandom(data = st[status == "n.s."], aes(colour = status), orientation = "y",
                   size = 0.45, alpha = 0.55, width = 0.38) +
  geom_quasirandom(data = st[status != "n.s."], aes(colour = status), orientation = "y",
                   size = 1.1, alpha = 0.95, width = 0.38) +
  geom_point(data = med, aes(med, lab), shape = 124, size = 5, colour = "grey10") +
  scale_colour_manual(values = st_cols, name = NULL, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(title = "log2FC by COG category",
       subtitle = sprintf("each point = one detected protein; | = category median; >0 = %s, <0 = %s",
                          dir_labels(tg, lower = TRUE)[["Up"]], dir_labels(tg, lower = TRUE)[["Down"]]),
       x = sprintf("log2FC (%s / %s)", tg, cg), y = NULL) +
  cog_theme() + theme(legend.position = "top")

# ---- F v2: ORA bubble, category x direction ---------------------------------
dir_lab <- setNames(dir_labels(tg), c("up", "down"))
ob <- rbindlist(lapply(names(dir_lab), function(dr) {
  x <- rd(sprintf("prot__enrich__%s_%s__cog_ora_all.tsv", cn, dr))
  if (is.null(x) || !nrow(x)) return(NULL)
  x[, .(letter = ID, FoldEnrichment, Count, pvalue, p.adjust, direction = dir_lab[[dr]])]
}), fill = TRUE)
ob <- cog_axis(ob)
ob[, direction := factor(direction, levels = dir_lab)]
ob[, sig := factor(p.adjust <= PADJ, levels = c(TRUE, FALSE))]
p_cog_ora <- ggplot(ob, aes(FoldEnrichment, lab, fill = -log10(pvalue), size = Count, colour = sig)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
  geom_point(shape = 21, stroke = 0.5, alpha = 0.9, show.legend = TRUE) +
  scale_fill_gradient(low = "grey85", high = "#d73027", name = expression(-log[10]~italic(p))) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey70"),
                      labels = c(`TRUE` = "FDR<=0.05", `FALSE` = "n.s."), name = "BH", drop = FALSE) +
  scale_size_continuous(range = c(2, 7), name = "DEPs in\ncategory") +
  guides(colour = guide_legend(override.aes = list(size = 3.5, fill = "white", alpha = 1)),
         fill = guide_colourbar(raster = FALSE)) +   # vector colourbar (Affinity-safe)
  facet_grid(group ~ direction, scales = "free_y", space = "free_y") +
  labs(title = "COG over-representation",
       subtitle = sprintf("universe = detected proteome with a COG letter (n=%d); outline = BH FDR<=0.05", annot_n),
       x = "fold enrichment", y = NULL) +
  cog_theme()

# ===========================================================================
# Assemble + save
# ===========================================================================
# Heatmap-dominant arrangement v2 (15-col grid): A|B|C top band (rows 1-4); D
# heatmap as a BIG tall left column (rows 5-15, cols 1-7 - one col wider than v1,
# which also uncramps its sample labels) with a taller canvas; E GO GSEA (top)
# and F COG (bottom) stacked in the RIGHT column (cols 8-15). Left column (A over
# D) aligns at cols 1-7; right region (B,C over E,F) at cols 8-15.
design <- c(
  area(1,  1,  4,  7),   # A PCA     (cols 1-7)
  area(1,  8,  4, 12),   # B Volcano (cols 8-12)
  area(1, 13,  4, 15),   # C Venn    (cols 13-15)
  area(5,  1, 15,  7),   # D heatmap (all proteins) - big tall left column, rows 5-15
  area(5,  8,  9, 15),   # E GO GSEA (right column, top)
  area(10, 8, 15, 15))   # F COG (v1 strip / v2 ora) (right column, bottom)

make_panel <- function(f_elem) {
  wrap_plots(
    p_pca, p_volc, wrap_elements(full = g_venn),
    wrap_elements(full = g_heat), p_go, f_elem,
    design = design) +
    plot_annotation(
      tag_levels = "A",
      title = bquote(bolditalic("Aspergillus fumigatus") ~
                     bold("proteome response to lactate (" * .(cn_lab) * ")")),
      theme = theme(plot.title = element_text(face = "bold", size = 15,
                                              margin = margin(b = 4)))) &
    theme(plot.tag = element_text(face = "bold", size = 15))
}

W_mm <- 380; H_mm <- 440   # heatmap-dominant v2: bigger (wider + taller) left column for D
versions <- list(
  list(suffix = "v1_cog_lfc_strip",  elem = p_cog_strip),
  list(suffix = "v2_cog_bubble_ora", elem = p_cog_ora))

for (v in versions) {
  base <- file.path(panel_dir, sprintf("prot__panel__%s__%s", cn, v$suffix))
  pnl  <- make_panel(v$elem)
  svglite::svglite(paste0(base, ".svg"), width = W_mm / 25.4, height = H_mm / 25.4)
  print(pnl); dev.off()
  ragg::agg_png(paste0(base, ".png"), width = W_mm, height = H_mm, units = "mm", res = 600)
  print(pnl); dev.off()
  message("DONE: ", base, ".{svg,png}")
}

message(sprintf("  A PCA | B volcano | C detected-venn(Lac=%d,Sal=%d) | D heatmap(n=%d quantified) | E GO GSEA(%d sets) | F = v1 COG strip(%d cats) / v2 COG ORA(%d rows)",
                length(det_sets$Lactate), length(det_sets$Saline), nrow(all_de),
                nrow(gse), length(present), nrow(ob)))
