#!/usr/bin/env Rscript
# ============================================================================
# 01_DE_proteome_limma.R  -  FragPipe MaxLFQ -> limma moderated-t DE.
# Strategy identical to reference Wt-Vs-Wt_MultiOmics/proteomics 01_DE_proteome.R:
#   quantile normalization + Perseus left-shifted Gaussian imputation, then limma.
# Input: biological matrix from io.R (DIA, N=4 per group).
# Significance (locked): padj <= 0.05, NO log2FC cutoff.
# LIMMA-ONLY analysis (proDA removed for this DIA folder per user request).
#
# Run:  Rscript scripts/01_DE_proteome_limma.R
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(limma) })
args <- commandArgs(FALSE); self <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
SDIR <- dirname(self); source(file.path(SDIR, "config.R"))
source(file.path(SDIR, "io.R")); source(file.path(SDIR, "plot_utils.R"))
ENGINE <- "limma"

ss <- fread(CFG$samples_tsv, sep = "\t")
ss[, group := factor(group, levels = CFG$group_levels)]
pm <- load_biological_matrix(CFG, ss)
log2mat <- pm$log2mat; meta <- pm$meta; grp <- pm$group
sg <- pm$sample_group[match(colnames(log2mat), sample_id)]
message(sprintf("[limma] input=%d proteins -> kept=%d (removed %d) x %d biol. samples",
                pm$n_input, nrow(log2mat), pm$n_removed, ncol(log2mat)))

# ---- STEP 2: quantile normalization ---------------------------------------
norm_mat <- normalizeBetweenArrays(log2mat, method = "quantile")

# ---- STEP 3: Perseus left-shifted Gaussian imputation (per sample column) --
set.seed(CFG$seed_limma_impute)
impute_perseus <- function(x, downshift, width) {
  miss <- is.na(x); if (!any(miss)) return(x)
  mu <- mean(x, na.rm = TRUE); sdv <- sd(x, na.rm = TRUE)
  x[miss] <- rnorm(sum(miss), mean = mu - downshift * sdv, sd = width * sdv); x
}
imp_mat <- apply(norm_mat, 2, impute_perseus,
                 downshift = CFG$impute_downshift, width = CFG$impute_width)
rownames(imp_mat) <- rownames(norm_mat)

# ---- STEP 4: DE (limma), padj-only significance ---------------------------
condition <- factor(sg$group, levels = CFG$group_levels)
# Optional batch adjustment: include acquisition batch as a covariate. Estimable
# only when batch is not aliased with group (checked via design rank); otherwise
# fall back to ~0+group with a warning (e.g. fully confounded subsets).
batch   <- ss$batch[match(colnames(imp_mat), ss$sample_id)]
batch_f <- factor(make.names(as.character(batch)))
use_batch <- isTRUE(CFG$batch_adjust) && nlevels(batch_f) > 1
if (use_batch) {
  design <- model.matrix(~ 0 + condition + batch_f)
  colnames(design) <- sub("^condition", "", colnames(design))  # bare group levels
  if (qr(design)$rank < ncol(design)) {
    warning("batch is confounded with group -> cannot adjust; using ~0+group")
    design <- model.matrix(~ 0 + condition); colnames(design) <- levels(condition)
    use_batch <- FALSE
  } else message("[limma] BATCH-ADJUSTED design (~0 + group + batch): batches = ",
                 paste(levels(batch_f), collapse = ", "))
} else {
  design <- model.matrix(~ 0 + condition); colnames(design) <- levels(condition)
}
# limma arrayWeights: down-weight lower-quality samples instead of deleting them.
aw <- if (isTRUE(CFG$arrayweights)) arrayWeights(imp_mat, design) else NULL
if (!is.null(aw)) {
  awt <- data.table(sample_id = colnames(imp_mat), group = as.character(sg$group),
                    array_weight = round(as.numeric(aw), 3))[order(array_weight)]
  fwrite(awt, file.path(CFG$dir_results, "prot__qc__all__array_weights_limma.tsv"), sep = "\t")
  message("[limma] arrayWeights range ", paste(round(range(aw), 2), collapse = "-"),
          " (lowest: ", awt$sample_id[1], "=", awt$array_weight[1], ")")
}
fit <- lmFit(imp_mat, design, weights = aw)
# contrasts built from config (test - control), any number
contr_str <- setNames(
  vapply(CFG$contrasts, function(ct) paste0(ct$test, " - ", ct$control), ""),
  vapply(CFG$contrasts, `[[`, "", "name"))
cm <- makeContrasts(contrasts = contr_str, levels = design)
colnames(cm) <- names(contr_str)
fit2 <- eBayes(contrasts.fit(fit, cm))

de_list <- list(); summ <- list()
for (cn in colnames(cm)) {
  tt <- as.data.table(topTable(fit2, coef = cn, number = Inf, sort.by = "none"),
                      keep.rownames = "protein_id")
  setnames(tt, c("logFC", "P.Value", "adj.P.Val"), c("log2FC", "pvalue", "padj"))
  tt <- merge(tt, meta[, .(protein_id, gene_raw, entry_name, description)],
              by = "protein_id", all.x = TRUE)
  tt[, label := display_label(gene_raw, protein_id)]
  tt <- add_regulation(tt, CFG$padj)
  setcolorder(tt, c("protein_id","gene_raw","label","log2FC","AveExpr","pvalue","padj",
                    "regulation","description"))
  setorder(tt, padj)
  de_list[[cn]] <- tt
  fwrite(tt, file.path(CFG$dir_results, sprintf("prot__de__%s__results_limma.tsv", cn)), sep = "\t")
  fwrite(tt[regulation != "NS"],
         file.path(CFG$dir_results, sprintf("prot__de__%s__significant_limma.tsv", cn)), sep = "\t")
  summ[[cn]] <- data.table(contrast = cn, engine = ENGINE, tested = nrow(tt),
                           up = sum(tt$regulation == "Up"), down = sum(tt$regulation == "Down"))
  message(sprintf("  %s [limma]: tested=%d up=%d down=%d",
                  cn, nrow(tt), summ[[cn]]$up, summ[[cn]]$down))
}
summ_dt <- rbindlist(summ)
fwrite(summ_dt, file.path(CFG$dir_results, "prot__de__all__summary_limma.tsv"), sep = "\t")

# ---- STEP 5: figures ------------------------------------------------------
save_fig(plot_pca(imp_mat, sg$group, GROUP_COLORS, "Proteome PCA"),
         "qc", "all", "pca_limma", CFG)
# batch-corrected PCA (visualization ONLY; DE uses batch in the design, not this)
if (use_batch) {
  corr_mat <- removeBatchEffect(imp_mat, batch = batch_f,
                                design = model.matrix(~ condition))
  save_fig(plot_pca(corr_mat, sg$group, GROUP_COLORS,
                    "Proteome PCA (batch-corrected, viz only)"),
           "qc", "all", "pca_limma_batchcorr", CFG)
}
save_heatmap(plot_corr_heatmap(imp_mat, sg$group, GROUP_COLORS, "Sample correlation (limma)"),
             "qc", "all", "corr_heatmap_limma", CFG)
save_fig(plot_count_bar(summ_dt, "Proteome DE counts (limma)"),
         "de", "all", "count_bar_limma", CFG, width_mm = CFG$width_double_mm)
for (cn in names(de_list)) {
  d <- de_list[[cn]]
  W2 <- CFG$width_double_mm
  save_fig(plot_volcano(d, cn, CFG$padj, REG_COLORS), "de", cn, "volcano_limma", CFG,
           width_mm = W2, height_mm = CFG$height_mm * 1.35)
  save_fig(plot_ma(d, cn, CFG$padj, REG_COLORS),       "de", cn, "ma_limma", CFG, width_mm = W2)
  save_fig(plot_pval_hist(d, cn),                       "de", cn, "pval_hist_limma", CFG, width_mm = W2)
}
sets <- lapply(de_list, function(d) d[regulation != "NS"]$protein_id)
if (.have("UpSetR") && sum(vapply(sets, length, 0L) > 0) >= 2) {
  b <- file.path(CFG$dir_figures, "prot__de__all__overlap_limma")
  png(paste0(b, ".png"), width = 120, height = 90, units = "mm", res = CFG$dpi)
  print(plot_overlap(sets, "DE overlap (limma)")); dev.off()
}

# ---- STEP 6: matrices + objects + log -------------------------------------
wmat <- function(m, f) fwrite(data.table(protein_id = rownames(m), m), f, sep = "\t")
wmat(norm_mat, file.path(CFG$dir_results, "prot__norm__all__matrix_limma.tsv"))
wmat(imp_mat,  file.path(CFG$dir_results, "prot__imputed__all__matrix_limma.tsv"))
saveRDS(list(fit = fit2, de = de_list, imp = imp_mat, norm = norm_mat, samples = ss),
        file.path(CFG$dir_objects, "prot_de_limma_objects.rds"))
writeLines(c(
  sprintf("engine\tlimma"),
  sprintf("time\t%s", Sys.time()),
  sprintf("n_input\t%d", pm$n_input), sprintf("n_kept\t%d", nrow(log2mat)),
  sprintf("n_removed\t%d", pm$n_removed),
  sprintf("normalization\tquantile"),
  sprintf("imputation\tperseus_%.1fSD_%.1fw_seed%d", CFG$impute_downshift, CFG$impute_width, CFG$seed_limma_impute),
  sprintf("padj\t%.3g\tlfc_cutoff\t%d", CFG$padj, CFG$lfc_cutoff)),
  file.path(CFG$dir_logs, "01_limma_run.log"))
message("[limma] DONE. Results in ", CFG$dir_results)
