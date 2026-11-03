#!/usr/bin/env Rscript
# ============================================================================
# 08_global_heatmap_venn.R - DEP heatmap, all-protein heatmap, detection Venn.
#
# Uses the significant Lactate-vs-Saline DEPs from 01_DE_proteome_limma.R.
# The DEP ComplexHeatmap ("top180") shows every significant DEP, split by
# direction; a second heatmap shows every quantified protein. One VennDiagram
# shows proteins detected in >=2 Lactate vs >=2 Saline biological replicates.
#
# Outputs (PNG in figures/, PDF in figures/Figures PDF/):
#   prot__de__Lactate_vs_Saline__heatmap_top180.{png,pdf}
#   prot__de__Lactate_vs_Saline__heatmap_all_proteins.{png,pdf}
#   prot__detected__Lactate_vs_Saline__venn.{png,pdf}
#   results/prot__de__Lactate_vs_Saline__deps.tsv
#   results/prot__detected__Lactate_vs_Saline__membership.tsv
#   results/prot__detected__Lactate_vs_Saline__venn_regions.tsv
#   logs/08_global_heatmap_venn.log
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(VennDiagram)
  library(grid)
})

args <- commandArgs(FALSE)
self <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
SDIR <- if (length(self)) dirname(normalizePath(self)) else "scripts"
source(file.path(SDIR, "config.R"))
source(file.path(SDIR, "io.R"))

contrast <- CFG$contrasts[[1]]$name
contrast_label <- gsub("_", " ", contrast, fixed = TRUE)
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

de <- fread(de_file, sep = "\t", na.strings = c("", "NA"))
deps <- de[!is.na(padj) & padj <= CFG$padj & regulation %chin% c("Up", "Down")]
deps[, abs_log2FC_sort := abs(log2FC)]
setorder(deps, regulation, padj, -abs_log2FC_sort)
deps[, abs_log2FC_sort := NULL]
if (nrow(deps) < 2L) stop("Fewer than two DEPs at padj <= ", CFG$padj)

# A compact, Excel-ready DEP table retained separately from the full DE table.
deps_file <- file.path(
  CFG$dir_results,
  sprintf("prot__de__%s__deps.tsv", contrast)
)
fwrite(deps, deps_file, sep = "\t", na = "NA")

# ---- reusable two-set Venn helpers ----------------------------------------
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
pair_regions <- function(sets) {
  stopifnot(length(sets) == 2L, !is.null(names(sets)))
  shared <- intersect(sets[[1]], sets[[2]])
  data.table(
    venn_section = c(paste0(names(sets)[1], "_only"),
                     paste0(names(sets)[2], "_only"), "Shared"),
    size = c(length(setdiff(sets[[1]], sets[[2]])),
             length(setdiff(sets[[2]], sets[[1]])), length(shared)),
    protein_ids = c(
      paste(sort(setdiff(sets[[1]], sets[[2]])), collapse = "; "),
      paste(sort(setdiff(sets[[2]], sets[[1]])), collapse = "; "),
      paste(sort(shared), collapse = "; ")
    )
  )
}

save_pair_venn <- function(sets, labels, colors, title, subtitle, name) {
  stopifnot(length(sets) == 2L, length(labels) == 2L, length(colors) == 2L)
  colors <- unname(colors)
  # Equal-size circles (scaled = FALSE): keeps all three region counts legible
  # inside the diagram even when one set nearly contains the other, so no count
  # is pushed outside on a leader line. VennDiagram still draws the LARGER set on
  # the left, so the category labels below are positioned by set size, not order.
  vg <- venn.diagram(
    x = sets, filename = NULL, category.names = c("", ""),
    scaled = FALSE, ext.text = FALSE,
    fill = colors, alpha = 0.45, col = colors, lwd = 2.4, lty = "solid",
    cex = 1.4, fontface = "bold", fontfamily = "sans", fontcolour = "grey15",
    cat.cex = 0.01, cat.col = "white", margin = 0.12,
    print.mode = "raw", disable.logging = TRUE
  )
  draw_page <- function() {
    grid.newpage()
    grid.rect(gp = gpar(fill = "white", col = NA))
    grid.text(title, y = 0.965,
              gp = gpar(fontfamily = "sans", fontface = "bold", fontsize = 16))
    grid.text(subtitle, y = 0.922,
              gp = gpar(fontfamily = "sans", col = "grey35", fontsize = 10))
    # VennDiagram draws the larger set on the left regardless of list order.
    left_idx  <- which.max(lengths(sets))
    right_idx <- setdiff(seq_along(sets), left_idx)
    grid.text(sprintf("%s\n(n=%d)", labels[left_idx], length(sets[[left_idx]])),
              x = 0.30, y = 0.74,
              gp = gpar(fontfamily = "sans", fontface = "bold", fontsize = 11.5,
                        col = colors[left_idx]))
    grid.text(sprintf("%s\n(n=%d)", labels[right_idx], length(sets[[right_idx]])),
              x = 0.70, y = 0.74,
              gp = gpar(fontfamily = "sans", fontface = "bold", fontsize = 11.5,
                        col = colors[right_idx]))
    pushViewport(viewport(y = 0.42, width = 0.90, height = 0.74))
    grid.draw(vg)
    upViewport()
  }
  png(file.path(CFG$dir_figures, paste0(name, ".png")),
      width = 180, height = 150, units = "mm", res = CFG$dpi)
  draw_page(); dev.off()
  pdf(file.path(CFG$dir_figures_pdf, paste0(name, ".pdf")), width = 180 / 25.4, height = 150 / 25.4)
  draw_page(); dev.off()
}

# ---- Venn: detected proteins in Lactate vs Saline -------------------------
ss_input <- fread(CFG$samples_tsv, sep = "\t")
pm <- load_biological_matrix(CFG, ss_input)
det_mat <- pm$log2mat
det_group <- pm$group
n_rep <- min(table(det_group))   # biological reps per group (3 in this 3v3 design)
n_lactate <- rowSums(!is.na(det_mat[, det_group == "Lactate", drop = FALSE]))
n_saline <- rowSums(!is.na(det_mat[, det_group == "Saline", drop = FALSE]))
detected <- data.table(
  protein_id = rownames(det_mat),
  n_Lactate = n_lactate,
  n_Saline = n_saline,
  detected_Lactate = n_lactate >= CFG$min_valid_group,
  detected_Saline = n_saline >= CFG$min_valid_group
)
detected <- merge(
  pm$meta[, .(protein_id, gene_raw, entry_name, description)],
  detected, by = "protein_id", all.y = TRUE
)
detected[, venn_section := fifelse(
  detected_Lactate & detected_Saline, "Shared",
  fifelse(detected_Lactate, "Lactate_only", "Saline_only")
)]
setorder(detected, venn_section, protein_id)
fwrite(
  detected,
  file.path(CFG$dir_results,
            "prot__detected__Lactate_vs_Saline__membership.tsv"),
  sep = "\t"
)
det_sets <- list(
  Lactate = detected[detected_Lactate == TRUE, protein_id],
  Saline = detected[detected_Saline == TRUE, protein_id]
)
det_regions <- pair_regions(det_sets)
fwrite(
  det_regions,
  file.path(CFG$dir_results,
            "prot__detected__Lactate_vs_Saline__venn_regions.tsv"),
  sep = "\t"
)
save_pair_venn(
  det_sets,
  labels = c("Detected in Lactate", "Detected in Saline"),
  colors = c(GROUP_COLORS[["Lactate"]], GROUP_COLORS[["Saline"]]),
  title = "Detected proteome overlap - Lactate vs Saline",
  subtitle = sprintf("Detected in >=%d of %d biological replicates | filtered DIA proteome",
                     CFG$min_valid_group, n_rep),
  name = "prot__detected__Lactate_vs_Saline__venn"
)

# ---- global ComplexHeatmap of all DEPs ------------------------------------
mdt <- fread(matrix_file, sep = "\t", na.strings = c("", "NA"))
if (!"protein_id" %in% names(mdt)) stop("Matrix lacks protein_id: ", matrix_file)
abundance <- as.matrix(mdt[, setdiff(names(mdt), "protein_id"), with = FALSE])
storage.mode(abundance) <- "double"
rownames(abundance) <- mdt$protein_id

samples <- unique(fread(CFG$samples_tsv, sep = "\t")[, .(sample_id, group)])
samples[, group := factor(group, levels = CFG$group_levels)]
setorder(samples, group, sample_id)
sample_order <- samples$sample_id[samples$sample_id %in% colnames(abundance)]
if (!length(sample_order)) stop("No sample-sheet IDs found in the abundance matrix")

# per-protein z-score, clipped to [-2, 2]; constant rows -> 0 (no NaN)
zscore_clip <- function(hm) {
  row_mean <- rowMeans(hm, na.rm = TRUE)
  row_sd <- apply(hm, 1, sd, na.rm = TRUE)
  z <- sweep(hm, 1, row_mean, "-")
  valid_sd <- is.finite(row_sd) & row_sd > 0
  z[valid_sd, ] <- sweep(z[valid_sd, , drop = FALSE], 1, row_sd[valid_sd], "/")
  z[!valid_sd, ] <- 0
  z[z < -2] <- -2
  z[z > 2] <- 2
  z
}
# draw a ComplexHeatmap to PNG (figures/) + PDF (figures/Figures PDF/), shared legend layout
save_ht <- function(ht, name, width_mm, height_mm) {
  png(file.path(CFG$dir_figures, paste0(name, ".png")), width = width_mm, height = height_mm, units = "mm", res = CFG$dpi)
  draw(ht, merge_legends = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  pdf(file.path(CFG$dir_figures_pdf, paste0(name, ".pdf")), width = width_mm / 25.4, height = height_mm / 25.4)
  draw(ht, merge_legends = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

zmat <- zscore_clip(abundance[deps$protein_id, sample_order, drop = FALSE])

row_labels <- deps$label
row_labels[is.na(row_labels) | !nzchar(row_labels)] <- deps$protein_id[
  is.na(row_labels) | !nzchar(row_labels)
]
rownames(zmat) <- make.unique(row_labels)

condition_ordered <- factor(
  samples$group[match(colnames(zmat), samples$sample_id)],
  levels = CFG$group_levels
)
top_ha <- HeatmapAnnotation(
  Condition = condition_ordered,
  col = list(Condition = GROUP_COLORS[CFG$group_levels]),
  annotation_name_side = "left"
)
dl <- dir_labels()                                    # Higher / Lower in <test>
direction <- factor(
  ifelse(deps$regulation == "Up", dl[["Up"]], dl[["Down"]]),
  levels = unname(dl)
)
direction_colors <- setNames(c(REG_COLORS[["Up"]], REG_COLORS[["Down"]]), unname(dl))
ha_left <- rowAnnotation(
  Direction = direction,
  col = list(Direction = direction_colors),
  show_annotation_name = FALSE
)
ha_right <- rowAnnotation(
  log2FC = anno_barplot(
    deps$log2FC,
    baseline = 0,
    gp = gpar(fill = unname(direction_colors[as.character(direction)]), col = NA),
    width = unit(27, "mm")
  ),
  annotation_name_rot = 0
)

heat_col <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
ht <- Heatmap(
  zmat,
  name = "row z-score",
  col = heat_col,
  top_annotation = top_ha,
  left_annotation = ha_left,
  right_annotation = ha_right,
  row_split = direction,
  cluster_rows = TRUE,
  cluster_row_slices = TRUE,
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  column_split = condition_ordered,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_max_width = unit(32, "mm"),
  row_names_gp = gpar(fontsize = 5),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 8),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title = sprintf(
    "Global differential proteome: %s (n=%d)", contrast_label, nrow(deps)
  ),
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  heatmap_legend_param = list(at = c(-2, -1, 0, 1, 2))
)

height_mm <- max(240, 1.65 * nrow(deps) + 75)
width_mm <- 220
save_ht(ht, sprintf("prot__de__%s__heatmap_top180", contrast), width_mm, height_mm)

# ---- heatmap of ALL quantified proteins (same aesthetics) ------------------
# Every protein in the limma matrix (identified AND quantified: >=2 valid values
# in >=1 group), split Higher / Lower in <test> / not significant. Row names are
# hidden (too many to read); slices keep a fixed order, rows cluster within them.
NS_LAB <- "Not significant"
all_de <- de[protein_id %chin% rownames(abundance)]
all_de[, status := factor(
  fifelse(!is.na(padj) & padj <= CFG$padj & regulation == "Up", dl[["Up"]],
  fifelse(!is.na(padj) & padj <= CFG$padj & regulation == "Down", dl[["Down"]], NS_LAB)),
  levels = c(unname(dl), NS_LAB))]
setorder(all_de, status, padj)
zall <- zscore_clip(abundance[all_de$protein_id, sample_order, drop = FALSE])
status_colors <- c(direction_colors, setNames(REG_COLORS[["NS"]], NS_LAB))
n_status <- table(all_de$status)

ht_all <- Heatmap(
  zall,
  name = "row z-score",
  col = heat_col,
  top_annotation = top_ha,
  left_annotation = rowAnnotation(
    Direction = all_de$status,
    col = list(Direction = status_colors),
    show_annotation_name = FALSE
  ),
  right_annotation = rowAnnotation(
    log2FC = anno_barplot(
      all_de$log2FC,
      baseline = 0,
      gp = gpar(fill = unname(status_colors[as.character(all_de$status)]), col = NA),
      width = unit(27, "mm")
    ),
    annotation_name_rot = 0
  ),
  row_split = all_de$status,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  row_title = sprintf("%s\n(n=%d)", names(n_status), as.integer(n_status)),
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_gap = unit(1.5, "mm"),
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  column_split = condition_ordered,
  show_row_names = FALSE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 8),
  column_title = sprintf(
    "Identified proteome: %s (n=%d quantified proteins)", contrast_label, nrow(all_de)
  ),
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  heatmap_legend_param = list(at = c(-2, -1, 0, 1, 2)),
  use_raster = FALSE
)
save_ht(ht_all, sprintf("prot__de__%s__heatmap_all_proteins", contrast),
        width_mm = 220, height_mm = 260)

writeLines(c(
  "step\tglobal_heatmap_venn",
  sprintf("time\t%s", Sys.time()),
  sprintf("contrast\t%s", contrast),
  sprintf("padj\t%.3g", CFG$padj),
  "heatmap_matrix\tquantile-normalized_Perseus-imputed_row-zscore",
  "heatmap_set\tall significant DEPs from 01_DE_proteome_limma.R",
  sprintf("deps_up\t%d", sum(deps$regulation == "Up")),
  sprintf("deps_down\t%d", sum(deps$regulation == "Down")),
  sprintf("deps_total\t%d", nrow(deps)),
  sprintf("all_protein_heatmap\t%d quantified proteins (%s)", nrow(all_de),
          paste(sprintf("%s=%d", names(n_status), as.integer(n_status)), collapse = ", ")),
  sprintf("detected_lactate\t%d", length(det_sets$Lactate)),
  sprintf("detected_saline\t%d", length(det_sets$Saline)),
  sprintf("detected_shared\t%d", length(intersect(det_sets$Lactate, det_sets$Saline))),
  sprintf("detection_rule\t>= %d of %d biological replicates", CFG$min_valid_group, n_rep)
), file.path(CFG$dir_logs, "08_global_heatmap_venn.log"))

message("DONE: DEP heatmap + all-protein heatmap + detection Venn")
message(sprintf("  DEPs=%d", nrow(deps)))
message(sprintf("  detected Lactate=%d, Saline=%d, shared=%d",
                length(det_sets$Lactate), length(det_sets$Saline),
                length(intersect(det_sets$Lactate, det_sets$Saline))))
