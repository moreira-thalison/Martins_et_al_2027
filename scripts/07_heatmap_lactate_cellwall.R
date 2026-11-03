#!/usr/bin/env Rscript
# ============================================================================
# 07_heatmap_lactate_cellwall.R - Curated cell-wall/lactate protein heatmap.
#
# Uses the quantile-normalized, Perseus-imputed matrix and limma results made by
# 01_DE_proteome_limma.R. Rows are kept in the supplied functional-group and
# gene order; columns are ordered Saline then Lactate. Heatmap values are
# per-protein z-scores (relative abundance across the six samples).
#
# Outputs:
#   figures/prot__de__Lactate_vs_Saline__heatmap_cellwall_proteins.png
#   figures/Figures PDF/prot__de__Lactate_vs_Saline__heatmap_cellwall_proteins.pdf
#   logs/07_heatmap_lactate_cellwall.log
#
# Run: Rscript scripts/07_heatmap_lactate_cellwall.R
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

args <- commandArgs(FALSE)
self <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
SDIR <- if (length(self)) dirname(normalizePath(self)) else "scripts"
source(file.path(SDIR, "config.R"))

gene_sets <- list(
  Cell_wall_biosynthesis_remodeling = c(
    "fks1",
    "ags1", "ags2", "ags3",
    "chsA", "chsB", "chsC", "chsD",
    "chsF", "chsG", "csmA", "csmB",
    "gel1", "gel2", "gel4", "gel7",
    "glfA", "uge5",
    "uge3", "sph3", "agd3", "gtb3", "ega3"
  ),
  Cell_wall_integrity_signaling = c(
    "rom2", "rho1",
    "pkcA", "bck1", "mkk2", "mpkA",
    "rlmA"
  ),
  Lactate_metabolism = c(
    "AFUA_5G14800",
    "AFUA_1G17040",
    "AFUA_1G00510",
    "AFUA_1G17520"
  ),
  Pyruvate_ethanol_metabolism = c(
    "AFUA_3G11070",
    "AFUA_6G00750",
    "AFUA_7G01010",
    "AFUA_7G04380"
  )
)

# AF293 aliases from NCBI Gene/UniProt. The input table contains a mixture of
# symbols (e.g. gel1) and AFUA locus tags, so both are tried when matching.
symbol_to_locus <- c(
  fks1 = "AFUA_6G12400",
  ags1 = "AFUA_3G00910", ags2 = "AFUA_2G11270", ags3 = "AFUA_1G15440",
  chsA = "AFUA_2G01870", chsB = "AFUA_4G04180", chsC = "AFUA_5G00760",
  chsD = "AFUA_1G12600", chsF = "AFUA_8G05630", chsG = "AFUA_3G14420",
  csmA = "AFUA_2G13440", csmB = "AFUA_2G13430",
  gel1 = "AFUA_2G01170", gel2 = "AFUA_6G11390", gel4 = "AFUA_2G05340",
  gel7 = "AFUA_6G12410", glfA = "AFUA_3G12690", uge5 = "AFUA_5G10780",
  uge3 = "AFUA_3G07910", sph3 = "AFUA_3G07900", agd3 = "AFUA_3G07870",
  gtb3 = "AFUA_3G07860", ega3 = "AFUA_3G07890",
  rom2 = "AFUA_5G08550", rho1 = "AFUA_6G06900", pkcA = "AFUA_5G11970",
  bck1 = "AFUA_3G11080", mkk2 = "AFUA_1G05800", mpkA = "AFUA_4G13720",
  rlmA = "AFUA_3G08520"
)

contrast <- CFG$contrasts[[1]]$name
matrix_file <- file.path(CFG$dir_results, "prot__imputed__all__matrix_limma.tsv")
de_file <- file.path(
  CFG$dir_results,
  sprintf("prot__de__%s__results_limma.tsv", contrast)
)
required <- c(matrix_file, de_file, CFG$samples_tsv)
missing_files <- required[!file.exists(required)]
if (length(missing_files)) {
  stop(
    "Required input(s) missing; run 00_make_samplesheet.R and ",
    "01_DE_proteome_limma.R first:\n  ",
    paste(missing_files, collapse = "\n  ")
  )
}

# ---- panel definition and ID matching --------------------------------------
panel <- rbindlist(lapply(seq_along(gene_sets), function(i) {
  data.table(
    category = names(gene_sets)[i],
    category_order = i,
    gene_order = seq_along(gene_sets[[i]]),
    requested_gene = gene_sets[[i]]
  )
}))
panel[, locus_tag := requested_gene]
is_symbol <- !grepl(
  "^AFUA_[0-9]+G[0-9]+$",
  panel$requested_gene,
  ignore.case = TRUE
)
panel[is_symbol, locus_tag := unname(symbol_to_locus[requested_gene])]
if (anyNA(panel$locus_tag)) {
  stop("No AF293 locus mapping for: ",
       paste(panel[is.na(locus_tag), requested_gene], collapse = ", "))
}

de <- fread(de_file, sep = "\t", na.strings = c("", "NA"))
needed_de <- c("protein_id", "gene_raw", "description", "log2FC", "padj", "regulation")
if (!all(needed_de %in% names(de))) {
  stop("Unexpected limma result columns; missing: ",
       paste(setdiff(needed_de, names(de)), collapse = ", "))
}
de[, `:=`(
  gene_key = toupper(trimws(gene_raw)),
  protein_key = toupper(trimws(protein_id))
)]

# Resolve each requested member against both its symbol and its AFUA locus. If
# an annotation unexpectedly maps to multiple proteins, retain the most abundant
# limma row and record the number of matches in the companion table.
resolved <- lapply(seq_len(nrow(panel)), function(i) {
  keys <- unique(toupper(c(panel$requested_gene[i], panel$locus_tag[i])))
  hit <- de[gene_key %chin% keys | protein_key %chin% keys]
  if (!nrow(hit)) {
    return(data.table(
      protein_id = NA_character_, gene_raw = NA_character_,
      description = NA_character_, log2FC = NA_real_, padj = NA_real_,
      regulation = NA_character_, n_matches = 0L
    ))
  }
  if ("AveExpr" %in% names(hit)) setorder(hit, -AveExpr, padj, na.last = TRUE)
  hit[1, .(
    protein_id, gene_raw, description, log2FC, padj,
    regulation = as.character(regulation), n_matches = nrow(hit)
  )]
})
panel <- cbind(panel, rbindlist(resolved))
panel[, detected := !is.na(protein_id)]
panel[, status := fifelse(detected, "detected", "not detected after limma filtering")]
panel[, display_label := fifelse(
  is_symbol,
  sprintf("%s (%s)", requested_gene, locus_tag),
  requested_gene
)]
setorder(panel, category_order, gene_order)

# ---- abundance matrix and sample ordering ----------------------------------
mdt <- fread(matrix_file, sep = "\t")
if (!"protein_id" %in% names(mdt)) stop("Matrix lacks a protein_id column: ", matrix_file)
abundance <- as.matrix(mdt[, setdiff(names(mdt), "protein_id"), with = FALSE])
storage.mode(abundance) <- "double"
rownames(abundance) <- mdt$protein_id

samples <- unique(fread(CFG$samples_tsv, sep = "\t")[, .(sample_id, group)])
samples[, group := factor(group, levels = CFG$group_levels)]
setorder(samples, group, sample_id)
sample_order <- samples$sample_id[samples$sample_id %in% colnames(abundance)]
if (!length(sample_order)) stop("No sample-sheet sample IDs found in the limma matrix")

plot_panel <- panel[detected & protein_id %in% rownames(abundance)]
if (nrow(plot_panel) < 2L) {
  stop("Only ", nrow(plot_panel), " requested protein(s) detected; need at least 2 for a heatmap")
}
mat <- abundance[plot_panel$protein_id, sample_order, drop = FALSE]
rownames(mat) <- make.unique(plot_panel$display_label)

# Per-protein z-score; imputed limma values make this complete, but retain a
# defensive zero-variance treatment so a constant row cannot create NaN values.
row_mean <- rowMeans(mat, na.rm = TRUE)
row_sd <- apply(mat, 1, sd, na.rm = TRUE)
zmat <- sweep(mat, 1, row_mean, "-")
valid_sd <- is.finite(row_sd) & row_sd > 0
zmat[valid_sd, ] <- sweep(zmat[valid_sd, , drop = FALSE], 1, row_sd[valid_sd], "/")
zmat[!valid_sd, ] <- 0

# ---- annotations and heatmap ----------------------------------------------
sample_group <- factor(
  samples$group[match(colnames(zmat), samples$sample_id)],
  levels = CFG$group_levels
)
top_ann <- HeatmapAnnotation(
  Condition = sample_group,
  col = list(Condition = GROUP_COLORS[CFG$group_levels]),
  annotation_name_side = "left"
)

reg <- factor(plot_panel$regulation, levels = c("Up", "NS", "Down"))
lfc <- plot_panel$log2FC
lfc_limit <- max(1, ceiling(max(abs(lfc), na.rm = TRUE)))
sig_label <- fifelse(
  !is.na(plot_panel$padj) & plot_panel$padj <= CFG$padj,
  "*", ""
)
right_ann <- rowAnnotation(
  `log2FC\nLactate/Saline` = anno_barplot(
    lfc,
    baseline = 0,
    ylim = c(-lfc_limit, lfc_limit),
    gp = gpar(fill = unname(REG_COLORS[as.character(reg)]), col = NA),
    axis_param = list(side = "bottom", labels_rot = 0),
    width = unit(32, "mm")
  ),
  `padj<=0.05` = anno_text(sig_label, gp = gpar(fontsize = 12, fontface = "bold")),
  annotation_name_rot = 0,
  gap = unit(2, "mm")
)

category_labels <- c(
  Cell_wall_biosynthesis_remodeling = "Cell wall biosynthesis & remodeling",
  Cell_wall_integrity_signaling = "Cell wall integrity signaling",
  Lactate_metabolism = "Lactate metabolism",
  Pyruvate_ethanol_metabolism = "Pyruvate/ethanol metabolism"
)
row_split <- factor(
  unname(category_labels[plot_panel$category]),
  levels = unname(category_labels[names(gene_sets)])
)
heat_col <- colorRamp2(c(-2, 0, 2), c("#4575b4", "white", "#d73027"))
ht <- Heatmap(
  zmat,
  name = "row\nz-score",
  col = heat_col,
  na_col = "grey85",
  top_annotation = top_ann,
  right_annotation = right_ann,
  row_split = row_split,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 8),
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  column_title = "Cell-wall and lactate-related proteins",
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  heatmap_legend_param = list(at = c(-2, -1, 0, 1, 2))
)

# ---- save figure and run log ------------------------------------------------
fig_name <- sprintf("prot__de__%s__heatmap_cellwall_proteins", contrast)
height_mm <- max(125, 8 * nrow(plot_panel) + 35)
width_mm <- 235
png(file.path(CFG$dir_figures, paste0(fig_name, ".png")), width = width_mm, height = height_mm,
    units = "mm", res = CFG$dpi)
draw(ht, merge_legends = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()
pdf(file.path(CFG$dir_figures_pdf, paste0(fig_name, ".pdf")), width = width_mm / 25.4,
    height = height_mm / 25.4)
draw(ht, merge_legends = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()

missing_message <- if (any(!panel$detected)) {
  paste(
    panel[detected == FALSE, sprintf("%s [%s]", requested_gene, locus_tag)],
    collapse = "; "
  )
} else {
  "none"
}
writeLines(c(
  "step\theatmap_lactate_cellwall",
  sprintf("time\t%s", Sys.time()),
  sprintf("contrast\t%s", contrast),
  sprintf("panel_members\t%d", nrow(panel)),
  sprintf("detected\t%d", sum(panel$detected)),
  sprintf("not_detected\t%d", sum(!panel$detected)),
  sprintf("not_detected_members\t%s", missing_message),
  "matrix\tquantile-normalized, Perseus-imputed limma matrix",
  "heatmap_scale\tper-protein z-score"
), file.path(CFG$dir_logs, "07_heatmap_lactate_cellwall.log"))

message(sprintf(
  "[heatmap] detected %d/%d requested proteins; not detected: %d",
  sum(panel$detected), nrow(panel), sum(!panel$detected)
))
if (any(!panel$detected)) message("  ", missing_message)
message("[heatmap] figure: ", fig_name, ".{png,pdf}")
