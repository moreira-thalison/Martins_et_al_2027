#!/usr/bin/env Rscript
# =============================================================================
# 02_kegg_enrichment.R
# KEGG ORA + GSEA for A. fumigatus Af293 (KEGG code "afm") on the limma DE.
#   ORA  : up   = padj <= 0.05 & log2FC > 0 (higher in test)
#          down = padj <= 0.05 & log2FC < 0 (lower in test)
#          universe = every tested protein (detected proteome)
#   GSEA : every tested protein ranked by limma log2FC
# BH p.adjust <= 0.05, minGSSize 5, maxGSSize 500. IDs are UniProt.
# Needs internet (KEGG REST).
#
# Input  : results/prot__de__<contrast>__results_limma.tsv
# Outputs: results/prot__enrich__<contrast>_<up|down>__kegg_ora{_all,}.tsv
#          results/prot__enrich__<contrast>__kegg_gsea{_all,}.tsv
#          dotplots (only when a term passes FDR), logs/02_kegg_enrichment.log
#          (_all = full table, read by 02b; no suffix = FDR-passing terms)
# Run    : Rscript scripts/02_kegg_enrichment.R
#
# Steps
#   0. Setup
#   1. Enrichment per contrast
#      1.1 DE sets
#      1.2 ORA (up / down)
#      1.3 GSEA
#   2. Run log
# =============================================================================


# 0. Setup ---------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(enrichplot); library(ggplot2)
})
options(timeout = 300)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
SDIR <- if (length(script_file)) dirname(normalizePath(script_file)) else "scripts"
source(file.path(SDIR, "config.R"))
source(file.path(SDIR, "io.R"))
source(file.path(SDIR, "plot_utils.R"))
source(file.path(SDIR, "enrich_utils.R"))
PADJ <- CFG$padj


# 1. Enrichment per contrast ---------------------------------------------------
for (ct in CFG$contrasts) {
  cn <- ct$name
  message("\n=== KEGG enrichment: ", cn, " ===")

  ## 1.1 DE sets ----
  d    <- read_de_results(CFG, cn)
  sets <- enrichment_sets(d, PADJ)
  message(sprintf("  DE set (limma): up(=higher in %s)=%d  down(=lower in %s)=%d  universe=%d",
                  ct$test, length(sets$up), ct$test, length(sets$down), length(sets$universe)))

  ## 1.2 ORA (up / down) ----
  for (dir in c("up", "down")) {
    r <- run_ora("KEGG", sprintf("%s/%s", cn, dir), sets[[dir]], padj_cut = PADJ,
                 enrich_fn = function(genes)
                   enrichKEGG(gene = genes, organism = CFG$kegg_org, keyType = CFG$kegg_keytype,
                              universe = sets$universe, pvalueCutoff = 1, qvalueCutoff = 1,
                              pAdjustMethod = "BH",
                              minGSSize = CFG$gsea_minGS, maxGSSize = CFG$gsea_maxGS))
    if (!has_terms(r)) next
    save_enrichment(r, sprintf("%s_%s", cn, dir), "kegg_ora",
                    sprintf("KEGG ORA - %s (%s)", cn, dir), CFG)
  }

  ## 1.3 GSEA ----
  ranks <- ranked_log2fc(d)
  set.seed(CFG$seed_gsea)
  gse <- run_gsea("KEGG", cn, padj_cut = PADJ, gsea_fn = function()
    gseKEGG(geneList = ranks, organism = CFG$kegg_org, keyType = CFG$kegg_keytype,
            minGSSize = CFG$gsea_minGS, maxGSSize = CFG$gsea_maxGS,
            pvalueCutoff = 1, pAdjustMethod = "BH", eps = 0, verbose = FALSE))
  if (!is.null(gse))
    save_enrichment(gse, cn, "kegg_gsea", sprintf("KEGG GSEA - %s", cn), CFG, gsea = TRUE)
}


# 2. Run log -------------------------------------------------------------------
writeLines(c(
  "step\tkegg_enrichment",
  sprintf("time\t%s", Sys.time()),
  sprintf("organism\t%s\tkeytype\t%s", CFG$kegg_org, CFG$kegg_keytype),
  sprintf("padj\t%.3g\tminGSSize\t%d\tmaxGSSize\t%d", PADJ, CFG$gsea_minGS, CFG$gsea_maxGS),
  "de_set\tlimma_significant up/down"),
  file.path(CFG$dir_logs, "02_kegg_enrichment.log"))
message("\nDONE: KEGG ORA + GSEA")
