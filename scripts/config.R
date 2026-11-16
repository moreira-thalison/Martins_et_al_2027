#!/usr/bin/env Rscript
# =============================================================================
# config.R
# Single source of every path, design choice, threshold, colour and direction
# label used by the pipeline. Sourced first by every script; nothing else
# hard-codes paths or thresholds.
#
# Analysis: A. fumigatus Af293 conidia, DIA (FragPipe / IonQuant MaxLFQ
# combined_protein.tsv), limma only.
#
# Design: 3 Saline vs 3 Lactate biological replicates, renamed for the paper
# (2026-11-04) to Saline_R1-R3 / Lactate_R1-R3; old -> paper names are in
# config/sample_name_map.csv.
#   Saline  = {R1_rep2_1, R2_rep2_2, R3_rep2_3}  -> all acquired 2026-08-19
#   Lactate = {R2_2, R3_3, R4_4}                 -> all acquired 2026-08-07
#
# Batch correction is impossible: group is perfectly confounded with
# acquisition day. batch_adjust = TRUE is requested, but the rank check in
# 01_DE_proteome_limma.R detects the confound, warns, and falls back to
# ~0 + group. Reported DEPs are therefore unadjusted and exploratory. Enabling
# batch correction needs >= 2 Saline replicates from the 2026-08-07 batch.
#
# Sections
#   1. Project root
#   2. Analysis settings (CFG)
#   3. Output folders
#   4. Colours
#   5. Direction wording
# =============================================================================


# 1. Project root --------------------------------------------------------------
# Resolved from this file's location, so the project runs after `git clone`
# anywhere. Order: $DIA_PROJECT_ROOT, parent of the running script's folder,
# parent of this sourced file, working directory.
.resolve_root <- function() {
  env <- Sys.getenv("DIA_PROJECT_ROOT", unset = NA)
  if (!is.na(env) && nzchar(env)) return(normalizePath(env))
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(script)) return(normalizePath(file.path(dirname(script[1]), "..")))
  for (i in rev(seq_len(sys.nframe()))) {
    sourced <- sys.frame(i)$ofile
    if (!is.null(sourced)) return(normalizePath(file.path(dirname(sourced), "..")))
  }
  normalizePath(getwd())
}
ROOT <- .resolve_root()


# 2. Analysis settings (CFG) ---------------------------------------------------
CFG <- list(
  ## 2.1 Input ----
  # data/combined_protein.tsv is the FragPipe report already reduced to the six
  # analysed samples; keep_samples below only checks that they are present.
  # Falls back to the original FragPipe run folder if data/ is absent.
  root          = ROOT,
  input_protein = {
    p <- file.path(ROOT, "data", "combined_protein.tsv")
    if (file.exists(p)) p
    else "/Volumes/SSDthalis/Proteomes/Fred_Aspergillus/Af293_DIA_new_reps/combined_protein.tsv"
  },

  ## 2.2 Experimental design ----
  control      = "Saline",
  group_levels = c("Saline", "Lactate"),
  batch_adjust = TRUE,   # requested; auto-disabled by the rank check (confounded)
  keep_samples = c("Saline_R1", "Saline_R2", "Saline_R3",
                   "Lactate_R1", "Lactate_R2", "Lactate_R3"),
  acq_batch    = c(Saline = "2026-08-19", Lactate = "2026-08-07"),  # acquisition day per group
  contrasts    = list(
    list(name = "Lactate_vs_Saline", test = "Lactate", control = "Saline")
  ),

  ## 2.3 Filtering and differential abundance (limma) ----
  padj              = 0.05,   # significance by BH-adjusted p only
  lfc_cutoff        = 0,      # no effect-size filter
  min_valid_group   = 2L,     # keep a protein with >= 2 valid values in >= 1 group
  impute_downshift  = 1.8,    # Perseus imputation: left shift (SD units)
  impute_width      = 0.3,    # Perseus imputation: width (SD units)
  seed_limma_impute = 123L,
  arrayweights      = TRUE,   # limma arrayWeights: down-weight noisier samples

  ## 2.4 Enrichment ----
  # KEGG, A. fumigatus Af293
  kegg_org     = "afm",
  kegg_keytype = "uniprot",
  gsea_minGS   = 5L,
  gsea_maxGS   = 500L,
  seed_gsea    = 42L,
  # GO: curated UniProt/GOA annotations of the Af293 reference proteome, per
  # aspect. Universe = detected proteome (same as KEGG).
  go_proteome       = "UP000002530",
  go_aspects        = c("BP", "MF", "CC"),
  go_force_download = FALSE,   # TRUE = re-fetch even if the cache exists
  # COG letters: NCBI COG covers prokaryotes only, so use UniProt eggNOG
  # cross-references (same proteome as GO) -> eggNOG 5.0 categories at the
  # Eukaryota level (2759).
  cog_eggnog_level   = 2759L,
  cog_force_download = FALSE,

  ## 2.5 Figures ----
  dpi = 600, width_single_mm = 90, width_double_mm = 180, height_mm = 90
)


# 3. Output folders ------------------------------------------------------------
CFG$dir_results          <- file.path(ROOT, "results")
CFG$dir_figures          <- file.path(ROOT, "figures")                      # PNG figures
CFG$dir_figures_pdf      <- file.path(CFG$dir_figures, "Figures PDF")       # PDF copy of each figure
CFG$dir_figures_kegg     <- file.path(CFG$dir_figures, "KEGG Pathways")     # ggkegg maps, PNG (06)
CFG$dir_figures_kegg_pdf <- file.path(CFG$dir_figures_kegg, "KEGG Figures PDF")  # ggkegg maps, PDF (06)
CFG$dir_objects          <- file.path(ROOT, "objects")
CFG$dir_logs             <- file.path(ROOT, "logs")
CFG$dir_config           <- file.path(ROOT, "config")
CFG$dir_annotation       <- file.path(ROOT, "annotation")                   # cached GO / eggNOG tables
CFG$dir_cache            <- file.path(ROOT, "cache", "kegg_kgml")           # KGML downloads (06)
CFG$samples_tsv          <- file.path(CFG$dir_config, "samples.tsv")

for (d in c(CFG$dir_results, CFG$dir_figures, CFG$dir_figures_pdf, CFG$dir_figures_kegg,
            CFG$dir_figures_kegg_pdf, CFG$dir_objects, CFG$dir_logs, CFG$dir_config,
            CFG$dir_annotation, CFG$dir_cache))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)


# 4. Colours -------------------------------------------------------------------
GROUP_COLORS <- c(Saline = "#7f7f7f", Lactate = "#1b9e77")
REG_COLORS   <- c(Up = "#d73027", NS = "#bdbdbd", Down = "#4575b4")


# 5. Direction wording ---------------------------------------------------------
# Always from the TEST group's perspective: "Higher in Lactate" / "Lower in
# Lactate", never "Higher in Saline". Every figure label, facet, legend and
# subtitle takes its wording from here.

# lower = TRUE for mid-sentence use ("higher in Lactate"; the group keeps its case).
dir_labels <- function(test = CFG$contrasts[[1]]$test, lower = FALSE)
  c(Up   = sprintf("%s in %s", if (lower) "higher" else "Higher", test),
    Down = sprintf("%s in %s", if (lower) "lower"  else "Lower",  test))

# Regulation (Up / NS / Down) -> factor of display labels, in legend order
# Higher, Lower, <ns_label>.
direction_status <- function(regulation, ns_label = "NS", test = CFG$contrasts[[1]]$test) {
  dl <- dir_labels(test)
  label <- ifelse(regulation == "Up", dl[["Up"]],
           ifelse(regulation == "Down", dl[["Down"]], ns_label))
  factor(label, levels = c(dl[["Up"]], dl[["Down"]], ns_label))
}
