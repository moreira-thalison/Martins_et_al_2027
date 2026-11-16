#!/usr/bin/env Rscript
# =============================================================================
# 03_qc_rle.R
# Relative log expression (RLE) per sample: each protein's median-normalised
# log2 abundance minus its across-sample median. Boxes centred on 0 with a
# similar spread mean comparable samples.
#
# Inputs : data/combined_protein.tsv, config/samples.tsv
# Output : figures/prot__qc__all__rle_box.png (+ PDF in figures/Figures PDF/)
# Run    : Rscript scripts/03_qc_rle.R
#
# Steps
#   0. Setup
#   1. Load the filtered matrix
#   2. Compute RLE
#   3. Plot
# =============================================================================


# 0. Setup ---------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table); library(matrixStats); library(ggplot2) })
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
SDIR <- if (length(script_file)) dirname(normalizePath(script_file)) else "scripts"
source(file.path(SDIR, "config.R"))
source(file.path(SDIR, "io.R"))
source(file.path(SDIR, "plot_utils.R"))


# 1. Load the filtered matrix --------------------------------------------------
ss <- fread(CFG$samples_tsv)
ss[, group := factor(group, levels = CFG$group_levels)]
pm <- load_biological_matrix(CFG, ss)
M  <- pm$log2mat


# 2. Compute RLE ---------------------------------------------------------------
median_norm <- function(m) sweep(m, 2, apply(m, 2, median, na.rm = TRUE))
Mn  <- median_norm(M)                                  # loading-normalised
rle <- Mn - rowMedians(Mn, na.rm = TRUE)

rle_dt <- melt(as.data.table(rle, keep.rownames = "protein_id"), id.vars = "protein_id",
               variable.name = "sample_id", value.name = "rle")
rle_dt[, group := pm$group[match(sample_id, colnames(M))]]
rle_dt[, sample_id := factor(sample_id, levels = colnames(M))]


# 3. Plot ----------------------------------------------------------------------
p <- ggplot(rle_dt, aes(sample_id, rle, fill = group)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
  scale_fill_manual(values = GROUP_COLORS) +
  labs(title = "Relative log expression (RLE) per sample", x = NULL, y = "RLE (log2)") +
  theme_mo(8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p, "qc", "all", "rle_box", CFG, width_mm = CFG$width_double_mm)
message("[QC] DONE. figure figures/prot__qc__all__rle_box")
