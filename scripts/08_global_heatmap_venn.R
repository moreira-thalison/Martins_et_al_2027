#!/usr/bin/env Rscript
# =============================================================================
# 08_global_heatmap_venn.R
# DEP heatmap, all-protein heatmap and detected-proteome Venn for the contrast.
#   DEP heatmap ("top180")  : every significant DEP from 01, split by direction
#   All-protein heatmap     : every quantified protein, split Higher / Lower /
#                             not significant (row names hidden)
#   Venn                    : proteins detected in >= 2 Lactate vs >= 2 Saline
#                             biological replicates
#
# Inputs : results/prot__imputed__all__matrix_limma.tsv,
#          results/prot__de__<contrast>__results_limma.tsv, config/samples.tsv,
#          data/combined_protein.tsv (detection counts, before imputation)
# Outputs (PNG in figures/, PDF in figures/Figures PDF/):
#          prot__de__<contrast>__heatmap_top180.{png,pdf}
#          prot__de__<contrast>__heatmap_all_proteins.{png,pdf}
#          prot__detected__Lactate_vs_Saline__venn.{png,pdf}
#          results/prot__de__<contrast>__deps.tsv
#          results/prot__detected__Lactate_vs_Saline__membership.tsv
#          results/prot__detected__Lactate_vs_Saline__venn_regions.tsv
#          logs/08_global_heatmap_venn.log
# Run    : Rscript scripts/08_global_heatmap_venn.R
#
# Steps
#   0. Setup
#   1. DEP table
#   2. Detected-proteome Venn
#      2.1 Venn helpers
#      2.2 Detection per group
#      2.3 Tables and figure
#   3. DEP heatmap
#   4. All-protein heatmap
#   5. Run log
# =============================================================================


# 0. Setup ---------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(ComplexHeatmap); library(circlize); library(VennDiagram); library(grid)
})
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
SDIR <- if (length(script_file)) dirname(normalizePath(script_file)) else "scripts"
source(file.path(SDIR, "config.R"))
source(file.path(SDIR, "io.R"))
source(file.path(SDIR, "plot_utils.R"))
invisible(futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))

contrast       <- CFG$contrasts[[1]]$name
contrast_label <- gsub("_", " ", contrast, fixed = TRUE)
require_inputs(c(imputed_matrix_path(CFG), de_results_path(CFG, contrast), CFG$samples_tsv))

dl       <- dir_labels()                                   # Higher / Lower in <test>
NS_LAB   <- "Not significant"
HEAT_COL <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
DIRECTION_COLORS <- setNames(c(REG_COLORS[["Up"]], REG_COLORS[["Down"]]), unname(dl))


# 1. DEP table -----------------------------------------------------------------
# Compact, Excel-ready DEP table, sorted by direction, padj, then |log2FC|.
de   <- read_de_results(CFG, contrast, na.strings = c("", "NA"))
deps <- de[!is.na(padj) & padj <= CFG$padj & regulation %chin% c("Up", "Down")]
deps <- deps[order(regulation, padj, -abs(log2FC))]
if (nrow(deps) < 2L) stop("Fewer than two DEPs at padj <= ", CFG$padj)
fwrite(deps, file.path(CFG$dir_results, sprintf("prot__de__%s__deps.tsv", contrast)),
       sep = "\t", na = "NA")


# 2. Detected-proteome Venn ----------------------------------------------------
## 2.1 Venn helpers ----
# Size and members of the three regions of a two-set Venn.
pair_regions <- function(sets) {
  stopifnot(length(sets) == 2L, !is.null(names(sets)))
  only1  <- setdiff(sets[[1]], sets[[2]])
  only2  <- setdiff(sets[[2]], sets[[1]])
  shared <- intersect(sets[[1]], sets[[2]])
  data.table(
    venn_section = c(paste0(names(sets)[1], "_only"), paste0(names(sets)[2], "_only"), "Shared"),
    size         = c(length(only1), length(only2), length(shared)),
    protein_ids  = c(paste(sort(only1), collapse = "; "), paste(sort(only2), collapse = "; "),
                     paste(sort(shared), collapse = "; "))
  )
}

# Equal-size circles (scaled = FALSE) keep all three counts legible inside the
# diagram even when one set nearly contains the other. VennDiagram always draws
# the LARGER set on the left, so category labels are placed by set size.
save_pair_venn <- function(sets, labels, colors, title, subtitle, name) {
  stopifnot(length(sets) == 2L, length(labels) == 2L, length(colors) == 2L)
  colors <- unname(colors)
  vg <- venn.diagram(
    x = sets, filename = NULL, category.names = c("", ""),
    scaled = FALSE, ext.text = FALSE,
    fill = colors, alpha = 0.45, col = colors, lwd = 2.4, lty = "solid",
    cex = 1.4, fontface = "bold", fontfamily = "sans", fontcolour = "grey15",
    cat.cex = 0.01, cat.col = "white", margin = 0.12,
    print.mode = "raw", disable.logging = TRUE
  )
  left_idx  <- which.max(lengths(sets))
  right_idx <- setdiff(seq_along(sets), left_idx)
  set_label <- function(idx, x)
    grid.text(sprintf("%s\n(n=%d)", labels[idx], length(sets[[idx]])), x = x, y = 0.74,
              gp = gpar(fontfamily = "sans", fontface = "bold", fontsize = 11.5, col = colors[idx]))
  draw_page <- function() {
    grid.newpage()
    grid.rect(gp = gpar(fill = "white", col = NA))
    grid.text(title, y = 0.965, gp = gpar(fontfamily = "sans", fontface = "bold", fontsize = 16))
    grid.text(subtitle, y = 0.922, gp = gpar(fontfamily = "sans", col = "grey35", fontsize = 10))
    set_label(left_idx, 0.30)
    set_label(right_idx, 0.70)
    pushViewport(viewport(y = 0.42, width = 0.90, height = 0.74))
    grid.draw(vg)
    upViewport()
  }
  png(fig_path(name, "png", CFG), width = 180, height = 150, units = "mm", res = CFG$dpi)
  draw_page(); invisible(dev.off())
  pdf(fig_path(name, "pdf", CFG), width = mm2in(180), height = mm2in(150))
  draw_page(); invisible(dev.off())
}

## 2.2 Detection per group (filtered matrix, before imputation) ----
pm        <- load_biological_matrix(CFG, fread(CFG$samples_tsv, sep = "\t"))
det_mat   <- pm$log2mat
det_group <- pm$group
n_rep     <- min(table(det_group))                         # biological replicates per group
n_lactate <- rowSums(!is.na(det_mat[, det_group == "Lactate", drop = FALSE]))
n_saline  <- rowSums(!is.na(det_mat[, det_group == "Saline", drop = FALSE]))
detected <- data.table(
  protein_id       = rownames(det_mat),
  n_Lactate        = n_lactate,
  n_Saline         = n_saline,
  detected_Lactate = n_lactate >= CFG$min_valid_group,
  detected_Saline  = n_saline >= CFG$min_valid_group
)
detected <- merge(pm$meta[, .(protein_id, gene_raw, entry_name, description)],
                  detected, by = "protein_id", all.y = TRUE)
detected[, venn_section := fifelse(detected_Lactate & detected_Saline, "Shared",
                           fifelse(detected_Lactate, "Lactate_only", "Saline_only"))]
setorder(detected, venn_section, protein_id)
det_sets <- list(Lactate = detected[detected_Lactate == TRUE, protein_id],
                 Saline  = detected[detected_Saline == TRUE, protein_id])

## 2.3 Tables and figure ----
fwrite(detected, file.path(CFG$dir_results, "prot__detected__Lactate_vs_Saline__membership.tsv"), sep = "\t")
fwrite(pair_regions(det_sets),
       file.path(CFG$dir_results, "prot__detected__Lactate_vs_Saline__venn_regions.tsv"), sep = "\t")
save_pair_venn(
  det_sets,
  labels   = c("Detected in Lactate", "Detected in Saline"),
  colors   = c(GROUP_COLORS[["Lactate"]], GROUP_COLORS[["Saline"]]),
  title    = "Detected proteome overlap - Lactate vs Saline",
  subtitle = sprintf("Detected in >=%d of %d biological replicates | filtered DIA proteome",
                     CFG$min_valid_group, n_rep),
  name     = "prot__detected__Lactate_vs_Saline__venn"
)


# 3. DEP heatmap ---------------------------------------------------------------
imp       <- read_imputed_matrix(CFG)
abundance <- imp$abundance
condition <- imp$samples$group
top_ha <- HeatmapAnnotation(
  Condition = condition,
  col = list(Condition = GROUP_COLORS[CFG$group_levels]),
  annotation_name_side = "left"
)

zmat <- row_zscore(abundance[deps$protein_id, , drop = FALSE], clip = 2)
rownames(zmat) <- make.unique(display_label(deps$label, deps$protein_id))
direction <- factor(ifelse(deps$regulation == "Up", dl[["Up"]], dl[["Down"]]), levels = unname(dl))

ht <- Heatmap(
  zmat,
  name = "row z-score",
  col = HEAT_COL,
  top_annotation = top_ha,
  left_annotation = rowAnnotation(Direction = direction,
                                  col = list(Direction = DIRECTION_COLORS),
                                  show_annotation_name = FALSE),
  right_annotation = rowAnnotation(
    log2FC = anno_barplot(deps$log2FC, baseline = 0,
                          gp = gpar(fill = unname(DIRECTION_COLORS[as.character(direction)]), col = NA),
                          width = unit(27, "mm")),
    annotation_name_rot = 0),
  row_split = direction,
  cluster_rows = TRUE,
  cluster_row_slices = TRUE,
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  column_split = condition,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_max_width = unit(32, "mm"),
  row_names_gp = gpar(fontsize = 5),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 8),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title = sprintf("Global differential proteome: %s (n=%d)", contrast_label, nrow(deps)),
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  heatmap_legend_param = list(at = c(-2, -1, 0, 1, 2))
)
save_complex_heatmap(ht, sprintf("prot__de__%s__heatmap_top180", contrast),
                     width_mm = 220, height_mm = max(240, 1.65 * nrow(deps) + 75), CFG)


# 4. All-protein heatmap -------------------------------------------------------
# Every protein in the limma matrix (identified and quantified: >= 2 valid
# values in >= 1 group). Slices keep a fixed order; rows cluster within them.
all_de   <- status_ordered_proteins(de, rownames(abundance), NS_LAB)
zall     <- row_zscore(abundance[all_de$protein_id, , drop = FALSE], clip = 2)
n_status <- table(all_de$status)
status_colors <- c(DIRECTION_COLORS, setNames(REG_COLORS[["NS"]], NS_LAB))

ht_all <- Heatmap(
  zall,
  name = "row z-score",
  col = HEAT_COL,
  top_annotation = top_ha,
  left_annotation = rowAnnotation(Direction = all_de$status,
                                  col = list(Direction = status_colors),
                                  show_annotation_name = FALSE),
  right_annotation = rowAnnotation(
    log2FC = anno_barplot(all_de$log2FC, baseline = 0,
                          gp = gpar(fill = unname(status_colors[as.character(all_de$status)]), col = NA),
                          width = unit(27, "mm")),
    annotation_name_rot = 0),
  row_split = all_de$status,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  row_title = sprintf("%s\n(n=%d)", names(n_status), as.integer(n_status)),
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_gap = unit(1.5, "mm"),
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  column_split = condition,
  show_row_names = FALSE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 8),
  column_title = sprintf("Identified proteome: %s (n=%d quantified proteins)", contrast_label, nrow(all_de)),
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  heatmap_legend_param = list(at = c(-2, -1, 0, 1, 2)),
  use_raster = FALSE
)
save_complex_heatmap(ht_all, sprintf("prot__de__%s__heatmap_all_proteins", contrast),
                     width_mm = 220, height_mm = 260, CFG)


# 5. Run log -------------------------------------------------------------------
n_shared <- length(intersect(det_sets$Lactate, det_sets$Saline))
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
  sprintf("detected_shared\t%d", n_shared),
  sprintf("detection_rule\t>= %d of %d biological replicates", CFG$min_valid_group, n_rep)
), file.path(CFG$dir_logs, "08_global_heatmap_venn.log"))

message("DONE: DEP heatmap + all-protein heatmap + detection Venn")
message(sprintf("  DEPs=%d", nrow(deps)))
message(sprintf("  detected Lactate=%d, Saline=%d, shared=%d",
                length(det_sets$Lactate), length(det_sets$Saline), n_shared))
