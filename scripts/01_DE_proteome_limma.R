#!/usr/bin/env Rscript
# =============================================================================
# 01_DE_proteome_limma.R
# FragPipe MaxLFQ -> limma moderated-t differential abundance.
# Quantile normalisation + Perseus left-shifted Gaussian imputation, then limma
# with arrayWeights. Significance: BH padj <= 0.05, no log2FC cutoff.
#
# Inputs : data/combined_protein.tsv, config/samples.tsv
# Outputs: results/prot__de__<contrast>__results_limma.tsv     (all proteins)
#          results/prot__de__<contrast>__significant_limma.tsv (DEPs)
#          results/prot__de__all__summary_limma.tsv
#          results/prot__qc__all__array_weights_limma.tsv
#          results/prot__norm__all__matrix_limma.tsv, prot__imputed__all__matrix_limma.tsv
#          figures: PCA, sample correlation, DE counts, volcano, MA, p-value histogram
#          objects/prot_de_limma_objects.rds, logs/01_limma_run.log
# Run    : Rscript scripts/01_DE_proteome_limma.R
#
# Steps
#   0. Setup
#   1. Load and filter the protein matrix
#   2. Quantile normalisation
#   3. Perseus imputation
#   4. limma model
#   5. DE tables per contrast
#   6. Figures
#   7. Matrices, objects and run log
# =============================================================================


# 0. Setup ---------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table); library(limma) })
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
SDIR <- if (length(script_file)) dirname(normalizePath(script_file)) else "scripts"
source(file.path(SDIR, "config.R"))
source(file.path(SDIR, "io.R"))
source(file.path(SDIR, "plot_utils.R"))


# 1. Load and filter the protein matrix ----------------------------------------
ss <- fread(CFG$samples_tsv, sep = "\t")
ss[, group := factor(group, levels = CFG$group_levels)]
pm <- load_biological_matrix(CFG, ss)
log2mat <- pm$log2mat
meta    <- pm$meta
sg      <- pm$sample_group[match(colnames(log2mat), sample_id)]
message(sprintf("[limma] input=%d proteins -> kept=%d (removed %d) x %d biol. samples",
                pm$n_input, nrow(log2mat), pm$n_removed, ncol(log2mat)))


# 2. Quantile normalisation ----------------------------------------------------
norm_mat <- normalizeBetweenArrays(log2mat, method = "quantile")


# 3. Perseus imputation --------------------------------------------------------
# Per sample: missing values drawn from N(mean - downshift*SD, width*SD).
impute_perseus <- function(x, downshift, width) {
  miss <- is.na(x)
  if (!any(miss)) return(x)
  mu <- mean(x, na.rm = TRUE)
  sdv <- sd(x, na.rm = TRUE)
  x[miss] <- rnorm(sum(miss), mean = mu - downshift * sdv, sd = width * sdv)
  x
}
set.seed(CFG$seed_limma_impute)
imp_mat <- apply(norm_mat, 2, impute_perseus,
                 downshift = CFG$impute_downshift, width = CFG$impute_width)
rownames(imp_mat) <- rownames(norm_mat)


# 4. limma model ---------------------------------------------------------------
## 4.1 Design (batch covariate only if estimable) ----
# Batch is used only when it is not aliased with group (design rank check);
# otherwise fall back to ~0 + group with a warning. In this 3v3 design it is
# always aliased (see config.R).
condition <- factor(sg$group, levels = CFG$group_levels)
batch     <- ss$batch[match(colnames(imp_mat), ss$sample_id)]
batch_f   <- factor(make.names(as.character(batch)))
group_only_design <- function() {
  d <- model.matrix(~ 0 + condition)
  colnames(d) <- levels(condition)
  d
}
use_batch <- isTRUE(CFG$batch_adjust) && nlevels(batch_f) > 1
if (use_batch) {
  design <- model.matrix(~ 0 + condition + batch_f)
  colnames(design) <- sub("^condition", "", colnames(design))   # bare group levels
  if (qr(design)$rank < ncol(design)) {
    warning("batch is confounded with group -> cannot adjust; using ~0+group")
    design <- group_only_design()
    use_batch <- FALSE
  } else {
    message("[limma] BATCH-ADJUSTED design (~0 + group + batch): batches = ",
            paste(levels(batch_f), collapse = ", "))
  }
} else {
  design <- group_only_design()
}

## 4.2 Sample weights (down-weight noisier samples, none removed) ----
aw <- if (isTRUE(CFG$arrayweights)) arrayWeights(imp_mat, design) else NULL
if (!is.null(aw)) {
  awt <- data.table(sample_id = colnames(imp_mat), group = as.character(sg$group),
                    array_weight = round(as.numeric(aw), 3))[order(array_weight)]
  fwrite(awt, file.path(CFG$dir_results, "prot__qc__all__array_weights_limma.tsv"), sep = "\t")
  message("[limma] arrayWeights range ", paste(round(range(aw), 2), collapse = "-"),
          " (lowest: ", awt$sample_id[1], "=", awt$array_weight[1], ")")
}

## 4.3 Fit and contrasts (test - control, from CFG$contrasts) ----
fit <- lmFit(imp_mat, design, weights = aw)
contr_str <- setNames(
  vapply(CFG$contrasts, function(ct) paste0(ct$test, " - ", ct$control), ""),
  vapply(CFG$contrasts, `[[`, "", "name"))
cm <- makeContrasts(contrasts = contr_str, levels = design)
colnames(cm) <- names(contr_str)
fit2 <- eBayes(contrasts.fit(fit, cm))


# 5. DE tables per contrast ----------------------------------------------------
de_list <- list()
summ    <- list()
for (cn in colnames(cm)) {
  tt <- as.data.table(topTable(fit2, coef = cn, number = Inf, sort.by = "none"),
                      keep.rownames = "protein_id")
  setnames(tt, c("logFC", "P.Value", "adj.P.Val"), c("log2FC", "pvalue", "padj"))
  tt <- merge(tt, meta[, .(protein_id, gene_raw, entry_name, description)],
              by = "protein_id", all.x = TRUE)
  tt[, label := display_label(gene_raw, protein_id)]
  tt <- add_regulation(tt, CFG$padj)
  setcolorder(tt, c("protein_id", "gene_raw", "label", "log2FC", "AveExpr", "pvalue", "padj",
                    "regulation", "description"))
  setorder(tt, padj)
  de_list[[cn]] <- tt

  fwrite(tt, de_results_path(CFG, cn), sep = "\t")
  fwrite(tt[regulation != "NS"],
         file.path(CFG$dir_results, sprintf("prot__de__%s__significant_limma.tsv", cn)), sep = "\t")
  summ[[cn]] <- data.table(contrast = cn, engine = "limma", tested = nrow(tt),
                           up = sum(tt$regulation == "Up"), down = sum(tt$regulation == "Down"))
  message(sprintf("  %s [limma]: tested=%d up=%d down=%d",
                  cn, nrow(tt), summ[[cn]]$up, summ[[cn]]$down))
}
summ_dt <- rbindlist(summ)
fwrite(summ_dt, file.path(CFG$dir_results, "prot__de__all__summary_limma.tsv"), sep = "\t")


# 6. Figures -------------------------------------------------------------------
## 6.1 QC: PCA and sample correlation ----
save_fig(plot_pca(imp_mat, sg$group, GROUP_COLORS, "Proteome PCA"), "qc", "all", "pca_limma", CFG)
if (use_batch) {
  # batch-corrected PCA for visualisation only (DE uses batch in the design)
  corr_mat <- removeBatchEffect(imp_mat, batch = batch_f, design = model.matrix(~ condition))
  save_fig(plot_pca(corr_mat, sg$group, GROUP_COLORS, "Proteome PCA (batch-corrected, viz only)"),
           "qc", "all", "pca_limma_batchcorr", CFG)
}
save_heatmap(plot_corr_heatmap(imp_mat, sg$group, GROUP_COLORS, "Sample correlation (limma)"),
             "qc", "all", "corr_heatmap_limma", CFG)

## 6.2 DE: counts, volcano, MA, p-value histogram ----
W2 <- CFG$width_double_mm
save_fig(plot_count_bar(summ_dt, "Proteome DE counts (limma)"),
         "de", "all", "count_bar_limma", CFG, width_mm = W2)
for (cn in names(de_list)) {
  d <- de_list[[cn]]
  save_fig(plot_volcano(d, cn, CFG$padj, REG_COLORS), "de", cn, "volcano_limma", CFG,
           width_mm = W2, height_mm = CFG$height_mm * 1.35)
  save_fig(plot_ma(d, cn, CFG$padj, REG_COLORS), "de", cn, "ma_limma", CFG, width_mm = W2)
  save_fig(plot_pval_hist(d, cn), "de", cn, "pval_hist_limma", CFG, width_mm = W2)
}

## 6.3 Overlap between contrasts (only with >= 2 contrasts) ----
sets <- lapply(de_list, function(d) d[regulation != "NS"]$protein_id)
if (.have("UpSetR") && sum(vapply(sets, length, 0L) > 0) >= 2) {
  png(file.path(CFG$dir_figures, "prot__de__all__overlap_limma.png"),
      width = 120, height = 90, units = "mm", res = CFG$dpi)
  print(plot_overlap(sets, "DE overlap (limma)"))
  invisible(dev.off())
}


# 7. Matrices, objects and run log ---------------------------------------------
write_matrix <- function(m, f) fwrite(data.table(protein_id = rownames(m), m), f, sep = "\t")
write_matrix(norm_mat, file.path(CFG$dir_results, "prot__norm__all__matrix_limma.tsv"))
write_matrix(imp_mat,  imputed_matrix_path(CFG))
saveRDS(list(fit = fit2, de = de_list, imp = imp_mat, norm = norm_mat, samples = ss),
        file.path(CFG$dir_objects, "prot_de_limma_objects.rds"))

writeLines(c(
  "engine\tlimma",
  sprintf("time\t%s", Sys.time()),
  sprintf("n_input\t%d", pm$n_input),
  sprintf("n_kept\t%d", nrow(log2mat)),
  sprintf("n_removed\t%d", pm$n_removed),
  "normalization\tquantile",
  sprintf("imputation\tperseus_%.1fSD_%.1fw_seed%d",
          CFG$impute_downshift, CFG$impute_width, CFG$seed_limma_impute),
  sprintf("padj\t%.3g\tlfc_cutoff\t%d", CFG$padj, CFG$lfc_cutoff)),
  file.path(CFG$dir_logs, "01_limma_run.log"))
message("[limma] DONE. Results in ", CFG$dir_results)
