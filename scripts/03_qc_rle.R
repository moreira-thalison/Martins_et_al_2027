#!/usr/bin/env Rscript
# ============================================================================
# 03_qc_rle.R  -  Relative log expression (RLE) per sample: each protein's
# median-normalized log2 abundance minus its across-sample median. Boxes centred
# on 0 with similar spread = comparable samples.
#
# Produces:
#   figures/prot__qc__all__rle_box.png  (+ PDF in figures/Figures PDF/)
#
# Run:  Rscript scripts/03_qc_rle.R
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(matrixStats); library(ggplot2)
})
args <- commandArgs(FALSE); self <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
SDIR <- dirname(self); source(file.path(SDIR, "config.R"))
source(file.path(SDIR, "io.R")); source(file.path(SDIR, "plot_utils.R"))

median_norm <- function(m) sweep(m, 2, apply(m, 2, median, na.rm=TRUE))

ss <- fread(CFG$samples_tsv); ss[, group := factor(group, levels = CFG$group_levels)]
pm <- load_biological_matrix(CFG, ss)
M  <- pm$log2mat; grp <- pm$group

Mn  <- median_norm(M)                                  # loading-normalized for RLE
rle <- Mn - rowMedians(Mn, na.rm = TRUE)

rle_dt <- melt(as.data.table(rle, keep.rownames="protein_id"), id.vars="protein_id",
               variable.name="sample_id", value.name="rle")
rle_dt[, group := grp[match(sample_id, colnames(M))]]
rle_dt[, sample_id := factor(sample_id, levels = colnames(M))]
save_fig(ggplot(rle_dt, aes(sample_id, rle, fill=group)) +
           geom_hline(yintercept=0, linetype=2, colour="grey50") +
           geom_boxplot(outlier.size=0.3, linewidth=0.3) +
           scale_fill_manual(values=GROUP_COLORS) +
           labs(title="Relative log expression (RLE) per sample",
                x=NULL, y="RLE (log2)") + theme_mo(8) +
           theme(axis.text.x=element_text(angle=45, hjust=1)),
         "qc","all","rle_box", CFG, width_mm=CFG$width_double_mm)
message("[QC] DONE. figure figures/prot__qc__all__rle_box")
