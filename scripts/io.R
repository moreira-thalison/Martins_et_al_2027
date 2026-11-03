#!/usr/bin/env Rscript
# ============================================================================
# io.R  -  Load + clean the FragPipe MaxLFQ matrix, collapse technical reps to
# the biological level, then filter.  Shared by BOTH DE scripts (limma + proDA)
# so the protein set, identifiers, collapse, and filtering are GUARANTEED
# identical between engines (mirrors reference shared/proteome_io.R).
#
# Normalization is intentionally NOT done here: limma uses quantile,
# proDA uses median-shift, so each engine normalizes its own way.
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(matrixStats) })

META_COLS  <- c("Protein ID", "Gene", "Entry Name", "Description", "Organism", "Protein Length")
META_NAMES <- c("protein_id", "gene_raw", "entry_name", "description", "organism", "protein_length")

# load_biological_matrix(CFG, ss)
#   ss : sample sheet (from 00_make_samplesheet.R) with maxlfq_col, group, sample_id
# Returns list(log2mat, meta, group, sample_group, n_removed, n_input):
#   * log2mat : filtered, technical-rep-collapsed log2 MaxLFQ (rows=proteins,
#               cols=biological sample_id). FragPipe 0 ("not quantified") -> NA.
#   * meta    : per-protein annotation for kept rows.
#   * group   : group label per matrix column (biological sample).
#   * sample_group : data.table(sample_id, group) unique, ordered.
load_biological_matrix <- function(CFG, ss) {
  stopifnot(all(META_COLS %in% names(fread(CFG$input_protein, nrows = 0))))
  prot <- fread(CFG$input_protein, sep = "\t", header = TRUE, check.names = FALSE)

  meta <- prot[, ..META_COLS]; setnames(meta, META_NAMES)

  # --- injection-level MaxLFQ matrix (0 -> NA -> log2) ---------------------
  mlcols <- ss$maxlfq_col                    # plain var: data.table `..` needs it
  missing <- setdiff(mlcols, names(prot))
  if (length(missing)) stop("MaxLFQ columns not found:\n  ", paste(missing, collapse = "\n  "))
  inj <- as.matrix(prot[, ..mlcols])
  inj[inj == 0] <- NA                       # FragPipe writes 0 for "not quantified"
  log2inj <- log2(inj)
  colnames(log2inj) <- ss$maxlfq_col
  rownames(log2inj) <- meta$protein_id

  sid <- ss$sample_id[match(colnames(log2inj), ss$maxlfq_col)]

  # --- collapse technical (injection) reps -> biological sample -----------
  # Statistically correct fix for pseudoreplication: average injections of the
  # SAME biological sample. Mean on the log2 scale (geometric mean of intensity),
  # robust to a single high-intensity injection. A protein detected in >=1
  # injection is retained for that biological sample (na.rm = TRUE); NA only if
  # absent in ALL injections of that sample.
  bio_ids <- unique(ss$sample_id)                      # keep sheet order
  bio <- sapply(bio_ids, function(s) {
    cols <- which(sid == s)
    m <- log2inj[, cols, drop = FALSE]
    v <- rowMeans(m, na.rm = TRUE)
    v[is.nan(v)] <- NA_real_                            # all-NA -> NA (not NaN)
    v
  })
  colnames(bio) <- bio_ids
  rownames(bio) <- meta$protein_id
  grp <- as.character(ss$group[match(bio_ids, ss$sample_id)])

  # --- filter: drop contaminants/decoys; keep if >=N valid in >=1 group ----
  keep_rows <- !grepl("^contam_|reverse", meta$protein_id, ignore.case = TRUE) &
               !grepl("^contam_|reverse", meta$entry_name,  ignore.case = TRUE)
  valid_ok <- apply(bio, 1, function(x)
    any(tapply(!is.na(x), grp, sum) >= CFG$min_valid_group))
  keep <- keep_rows & valid_ok

  list(log2mat = bio[keep, , drop = FALSE],
       meta    = meta[keep],
       group   = grp,
       sample_group = unique(data.table(sample_id = bio_ids, group = grp)),
       n_removed = sum(!keep), n_input = nrow(bio))
}

# add_regulation: Up / Down / NS by padj only (locked rule, no LFC cutoff).
add_regulation <- function(dt, padj_cut) {
  dt <- copy(dt)
  dt[, regulation := "NS"]
  dt[!is.na(padj) & padj <= padj_cut & log2FC > 0, regulation := "Up"]
  dt[!is.na(padj) & padj <= padj_cut & log2FC < 0, regulation := "Down"]
  dt[, regulation := factor(regulation, levels = c("Up", "NS", "Down"))]
  dt[]
}

# gene label for plots: prefer a real gene symbol, else the locus/UniProt id.
display_label <- function(meta_row_gene, protein_id) {
  lab <- meta_row_gene
  lab[is.na(lab) | lab == ""] <- protein_id[is.na(lab) | lab == ""]
  lab
}
