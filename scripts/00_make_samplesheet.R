#!/usr/bin/env Rscript
# =============================================================================
# 00_make_samplesheet.R
# Build the sample sheet from the FragPipe MaxLFQ column names
# ("Saline_R1" ... "Lactate_R3"; paper names, see config/sample_name_map.csv).
# Each injection is one biological replicate (its own sample_id). The
# acquisition batch comes from CFG$acq_batch and is recorded for transparency;
# 01 tries to use it but it is confounded with group (see config.R).
#
# Input : data/combined_protein.tsv
# Output: config/samples.tsv
# Run   : Rscript scripts/00_make_samplesheet.R
#
# Steps
#   0. Setup
#   1. Parse MaxLFQ columns
#   2. Select the analysed samples
#   3. Write and report
# =============================================================================


# 0. Setup ---------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
SDIR <- if (length(script_file)) dirname(normalizePath(script_file)) else "scripts"
source(file.path(SDIR, "config.R"))


# 1. Parse MaxLFQ columns ------------------------------------------------------
header   <- names(fread(CFG$input_protein, sep = "\t", nrows = 0, check.names = FALSE))
maxlfq   <- grep(" MaxLFQ Intensity$", header, value = TRUE)
raw_name <- sub(" MaxLFQ Intensity$", "", maxlfq)     # e.g. Saline_R1
group    <- sub("_R[0-9].*$", "", raw_name)           # text before the first "_R<digit>"

ss <- data.table(maxlfq_col = maxlfq, raw_name = raw_name, group_raw = group,
                 batch = unname(CFG$acq_batch[group]), sample_id = raw_name)


# 2. Select the analysed samples -----------------------------------------------
ss <- ss[group_raw %in% CFG$group_levels]
if (!is.null(CFG$keep_samples)) {
  not_found <- setdiff(CFG$keep_samples, ss$sample_id)
  if (length(not_found)) stop("keep_samples not found: ", paste(not_found, collapse = ", "))
  ss <- ss[sample_id %in% CFG$keep_samples]
}
ss[, group := factor(group_raw, levels = CFG$group_levels)][, group_raw := NULL]
setorder(ss, group, batch, sample_id)
stopifnot(uniqueN(ss$sample_id) == nrow(ss))          # one injection per sample_id


# 3. Write and report ----------------------------------------------------------
fwrite(ss, CFG$samples_tsv, sep = "\t")

cat(sprintf("Groups: %s\n", paste(CFG$group_levels, collapse = ", ")))
cat("\n== sample map ==\n")
print(ss[, .(injection = raw_name, group, batch, sample_id)])
cat("\n== samples per group x batch ==\n")
print(dcast(ss[, .N, by = .(group, batch)], group ~ batch, value.var = "N", fill = 0))
cat("\nWrote:", CFG$samples_tsv, "\n")
