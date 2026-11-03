#!/usr/bin/env bash
# ============================================================================
# run_all.sh  -  Reproduce the AF293 DIA "3v3" subset analysis (exploratory).
# Lactate {R2,R3,R4} vs Saline {R1,R2,R3}_rep2, selected via CFG$keep_samples
# from the full new-reps matrix in data/.
#
#   bash scripts/run_all.sh
#   DIA_PROJECT_ROOT=/path bash scripts/run_all.sh   # override project root
#
# Dependencies: Rscript env/dependencies.R  (CRAN + Bioconductor).
# Core steps (00,01,03) run offline and fail hard. Enrichment/figure steps need
# internet (KEGG, and UniProt/eggNOG on first run) and are non-fatal: a failure
# is warned and the pipeline continues.
# Outputs land in results/ figures/ objects/ logs/ (git-ignored). Figures: PNG in
# figures/, PDF in "figures/Figures PDF/", ggkegg maps in "figures/KEGG Pathways/",
# composite panels in figures/panel/.
#
# NOTE: this is a CONFOUNDED, exploratory design — Saline is all batch-2 and
# Lactate all batch-1, so batch correction is impossible (see README / config.R).
# ============================================================================
set -uo pipefail
cd "$(dirname "$0")/.."
ROOT="$(pwd)"
mkdir -p logs cache/kegg_kgml

run_core() { echo ">> $1"; shift; Rscript "$@" || { echo "FATAL: core step failed"; exit 1; }; }
run_opt()  { echo ">> $1 (optional)"; shift; Rscript "$@" || echo "WARN: optional step failed/skipped -> $*"; }

echo "=== AF293 DIA 3v3 pipeline | $(date) ==="

# ---- core (offline) --------------------------------------------------------
run_core "00 sample sheet (keep_samples: 3 Saline rep2 + 3 Lactate)" scripts/00_make_samplesheet.R
run_core "01 limma DE (batch_adjust requested -> auto-disabled; ~0+group)" scripts/01_DE_proteome_limma.R
run_core "03 QC: RLE per sample"                                     scripts/03_qc_rle.R

# ---- enrichment (need internet: KEGG/GO) -----------------------------------
run_opt  "02 KEGG ORA + GSEA (afm)"             scripts/02_kegg_enrichment.R
run_opt  "02b KEGG bubble plots"                scripts/02b_kegg_bubbleplots.R
run_opt  "05 GO enrichment (UniProt GOA)"       scripts/05_go_enrichment.R
run_opt  "05b GO bubble plots"                  scripts/05b_go_bubbleplots.R
run_opt  "05c COG enrichment + plots (eggNOG)"  scripts/05c_cog_enrichment.R
# ggkegg reads/writes <pid>.xml in the working directory -> run inside the KGML cache
(cd cache/kegg_kgml && run_opt "06 ggkegg pathway maps (KGML)" "$ROOT/scripts/06_ggkegg_pathways.R")

# ---- downstream figures ----------------------------------------------------
run_opt  "07 cell-wall protein heatmap"         scripts/07_heatmap_lactate_cellwall.R
run_opt  "08 DEP heatmaps + detection Venn"     scripts/08_global_heatmap_venn.R
run_opt  "09b panel variants (COG strip / ORA)" scripts/09b_panel_variants.R

# ---- provenance: capture the exact environment of THIS run -----------------
Rscript -e 'writeLines(capture.output(sessionInfo()), file.path("logs","sessionInfo.txt"))' \
  && echo ">> wrote logs/sessionInfo.txt"
echo "=== DONE. See results/, figures/, logs/. ==="
