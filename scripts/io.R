#!/usr/bin/env Rscript
# =============================================================================
# io.R
# Reading and writing the pipeline's core tables. Sourced after config.R.
#
# Normalisation is deliberately NOT done here: 01_DE_proteome_limma.R applies
# quantile normalisation itself.
#
# Sections
#   1. Raw FragPipe matrix         load_biological_matrix()
#   2. DE result columns           add_regulation(), display_label()
#   3. Files written by step 01    de_results_path(), read_de_results(),
#                                  imputed_matrix_path(), read_imputed_matrix(),
#                                  require_inputs()
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(matrixStats) })


# 1. Raw FragPipe matrix -------------------------------------------------------
META_COLS  <- c("Protein ID", "Gene", "Entry Name", "Description", "Organism", "Protein Length")
META_NAMES <- c("protein_id", "gene_raw", "entry_name", "description", "organism", "protein_length")

# Load the FragPipe MaxLFQ matrix, collapse injections to biological samples and
# filter proteins.
#   ss : sample sheet from 00_make_samplesheet.R (maxlfq_col, group, sample_id)
# Returns list(log2mat, meta, group, sample_group, n_removed, n_input):
#   log2mat      filtered log2 MaxLFQ (rows = proteins, cols = biological
#                sample_id); FragPipe 0 ("not quantified") -> NA
#   meta         per-protein annotation of the kept rows
#   group        group label per matrix column
#   sample_group data.table(sample_id, group), unique, sheet order
load_biological_matrix <- function(CFG, ss) {
  ## 1.1 Read the report and its annotation columns ----
  stopifnot(all(META_COLS %in% names(fread(CFG$input_protein, nrows = 0))))
  prot <- fread(CFG$input_protein, sep = "\t", header = TRUE, check.names = FALSE)
  meta <- prot[, ..META_COLS]
  setnames(meta, META_NAMES)

  ## 1.2 Injection-level MaxLFQ matrix (0 -> NA -> log2) ----
  mlcols  <- ss$maxlfq_col                       # plain variable: data.table `..` needs it
  missing <- setdiff(mlcols, names(prot))
  if (length(missing)) stop("MaxLFQ columns not found:\n  ", paste(missing, collapse = "\n  "))
  inj <- as.matrix(prot[, ..mlcols])
  inj[inj == 0] <- NA
  log2inj <- log2(inj)
  colnames(log2inj) <- ss$maxlfq_col
  rownames(log2inj) <- meta$protein_id
  sid <- ss$sample_id[match(colnames(log2inj), ss$maxlfq_col)]

  ## 1.3 Collapse injections of the same biological sample ----
  # Mean on the log2 scale (geometric mean of intensity). A protein seen in
  # >= 1 injection keeps a value; NA only if absent in all injections.
  bio_ids <- unique(ss$sample_id)                # keep sheet order
  bio <- sapply(bio_ids, function(s) {
    v <- rowMeans(log2inj[, sid == s, drop = FALSE], na.rm = TRUE)
    v[is.nan(v)] <- NA_real_                     # all-NA -> NA (not NaN)
    v
  })
  colnames(bio) <- bio_ids
  rownames(bio) <- meta$protein_id
  grp <- as.character(ss$group[match(bio_ids, ss$sample_id)])

  ## 1.4 Filter: drop contaminants / decoys; keep >= N valid in >= 1 group ----
  not_contam <- !grepl("^contam_|reverse", meta$protein_id, ignore.case = TRUE) &
                !grepl("^contam_|reverse", meta$entry_name, ignore.case = TRUE)
  enough_valid <- apply(bio, 1, function(x) any(tapply(!is.na(x), grp, sum) >= CFG$min_valid_group))
  keep <- not_contam & enough_valid

  list(log2mat      = bio[keep, , drop = FALSE],
       meta         = meta[keep],
       group        = grp,
       sample_group = unique(data.table(sample_id = bio_ids, group = grp)),
       n_removed    = sum(!keep),
       n_input      = nrow(bio))
}


# 2. DE result columns ---------------------------------------------------------
# Up / Down / NS by adjusted p only (no log2FC cutoff).
add_regulation <- function(dt, padj_cut) {
  dt <- copy(dt)
  dt[, regulation := "NS"]
  dt[!is.na(padj) & padj <= padj_cut & log2FC > 0, regulation := "Up"]
  dt[!is.na(padj) & padj <= padj_cut & log2FC < 0, regulation := "Down"]
  dt[, regulation := factor(regulation, levels = c("Up", "NS", "Down"))]
  dt[]
}

# Plot label: the gene symbol when there is one, else the UniProt id.
display_label <- function(gene, protein_id) {
  missing <- is.na(gene) | gene == ""
  gene[missing] <- protein_id[missing]
  gene
}


# 3. Files written by step 01 --------------------------------------------------
de_results_path <- function(CFG, contrast)
  file.path(CFG$dir_results, sprintf("prot__de__%s__results_limma.tsv", contrast))

# Full limma table for one contrast. Extra arguments go to fread().
read_de_results <- function(CFG, contrast, ...)
  fread(de_results_path(CFG, contrast), sep = "\t", ...)

imputed_matrix_path <- function(CFG)
  file.path(CFG$dir_results, "prot__imputed__all__matrix_limma.tsv")

# Quantile-normalised, Perseus-imputed matrix with columns ordered Saline then
# Lactate (sample_id within group). Returns list(abundance, samples), where
# samples$group is aligned with the abundance columns.
read_imputed_matrix <- function(CFG) {
  path <- imputed_matrix_path(CFG)
  mdt  <- fread(path, sep = "\t", na.strings = c("", "NA"))
  if (!"protein_id" %in% names(mdt)) stop("Matrix lacks a protein_id column: ", path)
  abundance <- as.matrix(mdt[, setdiff(names(mdt), "protein_id"), with = FALSE])
  storage.mode(abundance) <- "double"
  rownames(abundance) <- mdt$protein_id

  samples <- unique(fread(CFG$samples_tsv, sep = "\t")[, .(sample_id, group)])
  samples[, group := factor(group, levels = CFG$group_levels)]
  setorder(samples, group, sample_id)
  samples <- samples[sample_id %in% colnames(abundance)]
  if (!nrow(samples)) stop("No sample-sheet sample IDs found in ", path)

  list(abundance = abundance[, samples$sample_id, drop = FALSE], samples = samples)
}

# Stop early, with the list of missing files, when steps 00 / 01 have not run.
require_inputs <- function(files) {
  missing <- files[!file.exists(files)]
  if (length(missing))
    stop("Required input(s) missing; run 00_make_samplesheet.R and ",
         "01_DE_proteome_limma.R first:\n  ", paste(missing, collapse = "\n  "))
}
