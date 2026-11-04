#!/usr/bin/env Rscript
# ============================================================================
# config.R  -  AF293 (A. fumigatus) DIA re-analysis. Same workflow as
# New_analysis_August (DDA) but LIMMA-ONLY (proDA removed per user request),
# applied to the FragPipe (IonQuant MaxLFQ) combined_protein.tsv.
# Sourced by every script. NOTHING else hard-codes paths or thresholds.
#
# DESIGN (new-reps run, 2026-08-25): R1-R4 = biological replicates ("Pool" was a
# technician mislabel). Saline was EXPANDED with 6 more biological replicates from
# a later acquisition date (2026-08-19, "rep2") -> Saline N=10, Lactate N=4.
#
# *** USER-SPECIFIED 3v3 run (2026-08-25) ***: exact replicate list requested.
#   Saline  N=3 = {R1_rep2_1, R2_rep2_2, R3_rep2_3}  -> ALL batch-2 (2026-08-19)
#   Lactate N=3 = {R2_2, R3_3, R4_4}                 -> ALL batch-1 (2026-08-07)
# Renamed for the paper (2026-11-04) to Saline_R1-R3 / Lactate_R1-R3 in data/;
# old -> paper names in config/sample_name_map.csv.
#
# *** BATCH CORRECTION IS IMPOSSIBLE HERE ***: group is PERFECTLY confounded with
# batch (Saline=only batch-2, Lactate=only batch-1; no shared batch). The design
# ~0+group+batch is rank-deficient, so batch_adjust=TRUE below is CAUGHT by the
# runtime rank check in 01, which warns and FALLS BACK to ~0+group (no correction).
# The reported DEPs are therefore UNADJUSTED and fully confounded with acquisition
# day -> exploratory only. To actually enable batch correction, the selection must
# include >=2 batch-1 Saline reps (so batch-1 holds both groups).
# ============================================================================

# ---- Portable project root ------------------------------------------------
# Resolve ROOT from this file's own location (not a hard-coded absolute path),
# so the project is reproducible after `git clone` anywhere. Order: (1)
# $DIA_PROJECT_ROOT override, (2) parent of the running script's dir, (3) parent
# of this sourced file, (4) cwd.
.resolve_root <- function() {
  env <- Sys.getenv("DIA_PROJECT_ROOT", unset = NA)
  if (!is.na(env) && nzchar(env)) return(normalizePath(env))
  a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  if (length(f)) return(normalizePath(file.path(dirname(f[1]), "..")))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(file.path(dirname(of), "..")))
  }
  normalizePath(getwd())
}
ROOT <- .resolve_root()

CFG <- list(
  # ---- Paths (all relative to the project ROOT) ----------------------------
  root          = ROOT,
  # In-project input: FragPipe combined_protein.tsv already reduced to the 6
  # analysed samples (full 14-run report in Fragpipe_results/); keep_samples
  # below now acts only as a check that these 6 are present.
  # Falls back to the external FragPipe run if data/ is absent.
  input_protein = {
    p <- file.path(ROOT, "data", "combined_protein.tsv")
    if (file.exists(p)) p
    else "/Volumes/SSDthalis/Proteomes/Fred_Aspergillus/Af293_DIA_new_reps/combined_protein.tsv"
  },
  base          = ROOT,

  # ---- Experimental design -------------------------------------------------
  # 2 groups, control = Saline (N=3). Lactate N=3. batch_adjust requested but NOT
  # estimable (perfect confound) -> pipeline falls back to ~0+group.
  control      = "Saline",
  group_levels = c("Saline", "Lactate"),
  batch_adjust = TRUE,   # requested; auto-disabled by rank check (confounded)
  keep_samples = c("Saline_R1", "Saline_R2", "Saline_R3",
                   "Lactate_R1", "Lactate_R2", "Lactate_R3"),
  # Acquisition date per group (sample names no longer encode it).
  acq_batch    = c(Saline = "2026-08-19", Lactate = "2026-08-07"),
  contrasts = list(
    list(name = "Lactate_vs_Saline", test = "Lactate", control = "Saline")
  ),

  # ---- Locked statistical parameters (identical to DDA workflow) -----------
  padj              = 0.05,   # significance by adjusted p only (padj <= 0.05)
  lfc_cutoff        = 0,      # NO effect-size filter
  min_valid_group   = 2L,     # keep protein if >=2 valid in >=1 group
  # limma path
  impute_downshift  = 1.8,    # Perseus left-shift (SD units)
  impute_width      = 0.3,
  seed_limma_impute = 123L,
  # proDA path
  seed_proda        = 42L,
  quant_min_both    = 2L,     # quantitative: >=2 valid in BOTH groups
  # KEGG (A. fumigatus Af293)
  kegg_org          = "afm",
  kegg_keytype      = "uniprot",
  gsea_minGS        = 5L,
  gsea_maxGS        = 500L,
  seed_gsea         = 42L,
  # GO (A. fumigatus Af293) - curated UniProt/GOA annotations, per aspect (BP/MF/CC).
  # Term->protein map built from the UniProt reference proteome; universe = detected
  # proteome (same as KEGG). No org.Af.eg.db / GO.db needed (names come from UniProt).
  go_proteome       = "UP000002530",  # UniProt reference proteome for Af293
  go_aspects        = c("BP", "MF", "CC"),
  go_force_download = FALSE,           # TRUE = re-fetch even if cache exists
  # COG functional categories (A-Z). NCBI COG covers prokaryotes only, so for this
  # fungus use the eukaryotic extension: UniProt eggNOG cross-refs (KOG / ENOG50
  # orthologous groups, same proteome as GO) -> eggNOG 5.0 OG category letters at
  # the Eukaryota level (2759). Universe = detected proteome (same as GO/KEGG).
  cog_eggnog_level   = 2759L,
  cog_force_download = FALSE,

  # ---- Sample weighting ----------------------------------------------------
  # limma arrayWeights: down-weight noisier samples (no sample is removed).
  arrayweights  = TRUE,

  # ---- Figures -------------------------------------------------------------
  dpi = 600, width_single_mm = 90, width_double_mm = 180, height_mm = 90
)

CFG$dir_data       <- file.path(CFG$base, "data")
CFG$dir_results    <- file.path(CFG$base, "results")
CFG$dir_figures    <- file.path(CFG$base, "figures")                 # PNG figures
CFG$dir_figures_pdf  <- file.path(CFG$dir_figures, "Figures PDF")     # PDF copies of every figure
CFG$dir_figures_kegg <- file.path(CFG$dir_figures, "KEGG Pathways")   # ggkegg maps, PNG (06)
CFG$dir_figures_kegg_pdf <- file.path(CFG$dir_figures_kegg, "KEGG Figures PDF")  # ggkegg maps, PDF (06)
CFG$dir_objects    <- file.path(CFG$base, "objects")
CFG$dir_logs       <- file.path(CFG$base, "logs")
CFG$dir_config     <- file.path(CFG$base, "config")
CFG$dir_annotation <- file.path(CFG$base, "annotation")      # cached GO map (UniProt)
CFG$dir_cache      <- file.path(CFG$base, "cache", "kegg_kgml")  # KGML download cache (06)
CFG$samples_tsv    <- file.path(CFG$dir_config, "samples.tsv")
for (d in c(CFG$dir_results, CFG$dir_figures, CFG$dir_figures_pdf, CFG$dir_figures_kegg,
            CFG$dir_figures_kegg_pdf,
            CFG$dir_objects, CFG$dir_logs, CFG$dir_config, CFG$dir_annotation, CFG$dir_cache))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# Group colours: saline grey, lactate green.
GROUP_COLORS <- c(Saline = "#7f7f7f", Lactate = "#1b9e77")
REG_COLORS <- c(Up = "#d73027", NS = "#bdbdbd", Down = "#4575b4")

# Direction wording, always from the TEST group's perspective (user, 2026-10-27):
# "Higher in <test>" / "Lower in <test>" - never "Higher in <control>". Every
# figure label, facet, legend and subtitle takes its wording from here.
# lower = TRUE for mid-sentence use ("higher in Lactate"; group name keeps its case).
dir_labels <- function(test = CFG$contrasts[[1]]$test, lower = FALSE)
  c(Up   = sprintf("%s in %s", if (lower) "higher" else "Higher", test),
    Down = sprintf("%s in %s", if (lower) "lower"  else "Lower",  test))
