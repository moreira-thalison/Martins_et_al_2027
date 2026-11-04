#!/usr/bin/env Rscript
# ============================================================================
# 00_make_samplesheet.R  (DIA new-reps)  -  Derive the sample sheet from the
# FragPipe MaxLFQ column names ("Saline_R1" ... "Lactate_R3"; paper names, see
# config/sample_name_map.csv). Acquisition batch comes from CFG$acq_batch.
#
# DESIGN: 2 groups, biological replicates. Saline expanded with a 2nd acquisition
# batch (2026-08-19, "rep2"). Each injection = one biological replicate (its own
# sample_id; no collapsing). Acquisition batch is recorded for transparency but
# is NOT used in the model (per user: no batch adjustment).
# ============================================================================
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(FALSE)
self <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
source(file.path(dirname(self), "config.R"))

hdr <- names(fread(CFG$input_protein, sep = "\t", nrows = 0, check.names = FALSE))
maxlfq <- grep(" MaxLFQ Intensity$", hdr, value = TRUE)
raw_name <- sub(" MaxLFQ Intensity$", "", maxlfq)             # e.g. Saline_R1
# group = everything before the first "_R<digit>" token (robust to _rep2_ suffixes)
group    <- sub("_R[0-9].*$", "", raw_name)                   # Lactate / Saline
batch    <- unname(CFG$acq_batch[group])
sample_id <- raw_name                                          # unique per injection

ss_all <- data.table(maxlfq_col = maxlfq, raw_name = raw_name, group_raw = group,
                     batch = batch, sample_id = sample_id)
ss <- ss_all[group_raw %in% CFG$group_levels]
# Optional whitelist: keep only the requested biological samples.
if (!is.null(CFG$keep_samples)) {
  miss <- setdiff(CFG$keep_samples, ss$sample_id)
  if (length(miss)) stop("keep_samples not found: ", paste(miss, collapse = ", "))
  ss <- ss[sample_id %in% CFG$keep_samples]
}
ss[, group := factor(group_raw, levels = CFG$group_levels)][, group_raw := NULL]
setorder(ss, group, batch, sample_id)
stopifnot(uniqueN(ss$sample_id) == nrow(ss))                  # 1 injection per sample_id

fwrite(ss, CFG$samples_tsv, sep = "\t")

cat(sprintf("Groups: %s\n", paste(CFG$group_levels, collapse = ", ")))
cat("\n== sample map ==\n")
print(ss[, .(injection = raw_name, group, batch, sample_id)])
cat("\n== samples per group x batch ==\n")
print(dcast(ss[, .N, by = .(group, batch)], group ~ batch, value.var = "N", fill = 0))
cat("\nWrote:", CFG$samples_tsv, "\n")
