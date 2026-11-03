#!/usr/bin/env Rscript
# ============================================================================
# 02_kegg_enrichment.R  -  KEGG ORA + GSEA for the A. fumigatus (Af293, KEGG
# code "afm") proteome, on the LIMMA DE (limma-only analysis). Adapted from the
# reference 02_kegg_enrichment.R (ORA up/down + GSEA).
#
# Differential gene set for ORA = limma-significant (padj <= 0.05):
#   up   = padj <= 0.05 & log2FC > 0   (higher in test)
#   down = padj <= 0.05 & log2FC < 0   (lower in test)
# GSEA ranks ALL tested proteins by limma log2FC (captures diffuse pathway signal).
# Thresholds: p.adjust<=0.05 (BH), minGSSize 5, maxGSSize 500; set.seed for GSEA.
# Requires internet (KEGG REST). IDs = UniProt (keyType="uniprot").
#
# Run:  Rscript scripts/02_kegg_enrichment.R
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(enrichplot); library(ggplot2)
})
options(timeout = 300)
args <- commandArgs(FALSE); self <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
SDIR <- dirname(self); source(file.path(SDIR, "config.R")); source(file.path(SDIR, "plot_utils.R"))
ORG <- CFG$kegg_org; KT <- CFG$kegg_keytype; PADJ <- CFG$padj

# Figures: dotplot only when terms pass FDR. No nominal-trend barplot (removed
# per user, 2026-10-27: bar charts are not used for enrichment) - the full tables
# are always written, so nominal trends stay available in results/.

ora <- function(genes, universe, label) {
  genes <- unique(genes[!is.na(genes)]); universe <- unique(universe[!is.na(universe)])
  if (length(genes) < 3) { message(sprintf("  [skip] %-32s <3 genes", label)); return(NULL) }
  # pvalueCutoff = 1 -> keep the FULL table; we filter/report FDR ourselves.
  r <- tryCatch(enrichKEGG(gene = genes, organism = ORG, keyType = KT, universe = universe,
                           pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
                           minGSSize = CFG$gsea_minGS, maxGSSize = CFG$gsea_maxGS),
                error = function(e) { message("  [error] ORA ", label, ": ", conditionMessage(e)); NULL })
  n_sig <- if (is.null(r)) 0 else sum(as.data.frame(r)$p.adjust <= PADJ)
  message(sprintf("  KEGG ORA %-32s : %d terms tested, %d pass FDR",
                  label, if (is.null(r)) 0 else nrow(as.data.frame(r)), n_sig)); r
}

for (ct in CFG$contrasts) {
  cn <- ct$name; tg <- ct$test; cg <- ct$control
  message("\n=== KEGG enrichment: ", cn, " ===")
  # LIMMA-only: DE set = limma-significant up/down (no proDA / no on/off layer).
  d  <- fread(file.path(CFG$dir_results, sprintf("prot__de__%s__results_limma.tsv", cn)))

  universe <- unique(d$protein_id)                       # all tested proteins
  up   <- d[!is.na(padj) & padj <= PADJ & log2FC > 0, protein_id]
  down <- d[!is.na(padj) & padj <= PADJ & log2FC < 0, protein_id]
  message(sprintf("  DE set (limma): up(=higher in %s)=%d  down(=lower in %s)=%d  universe=%d",
                  tg, length(up), tg, length(down), length(universe)))

  # ---- ORA: up / down -----------------------------------------------------
  for (dir in c("up","down")) {
    genes <- switch(dir, up = up, down = down)
    r <- ora(genes, universe, sprintf("%s/%s", cn, dir))
    if (is.null(r) || nrow(as.data.frame(r)) == 0) next
    rdf <- as.data.frame(r)
    # full (unfiltered) table: read by 02b_kegg_bubbleplots.R
    fwrite(rdf, file.path(CFG$dir_results,
           sprintf("prot__enrich__%s_%s__kegg_ora_all.tsv", cn, dir)), sep = "\t")
    sig <- rdf[rdf$p.adjust <= PADJ, , drop = FALSE]
    if (nrow(sig)) {                       # FDR-significant subset (+ dotplot)
      fwrite(sig, file.path(CFG$dir_results,
             sprintf("prot__enrich__%s_%s__kegg_ora.tsv", cn, dir)), sep = "\t")
      p <- enrichplot::dotplot(r, showCategory = 20) +
             ggtitle(sprintf("KEGG ORA - %s (%s)", cn, dir)) + theme_mo(8)
      save_fig(p, "enrich", sprintf("%s_%s", cn, dir), "kegg_ora", CFG,
               width_mm = CFG$width_double_mm, height_mm = CFG$height_mm * 1.3)
    }                                       # else: table only (no barplot)
  }

  # ---- GSEA on ranked quantitative log2FC ---------------------------------
  dd <- d[is.finite(log2FC) & !is.na(protein_id)]
  ranks <- dd$log2FC; names(ranks) <- dd$protein_id
  ranks <- sort(ranks[!duplicated(names(ranks))], decreasing = TRUE)
  set.seed(CFG$seed_gsea)
  # pvalueCutoff = 1 -> keep the FULL ranked pathway table; report FDR ourselves.
  gse <- tryCatch(gseKEGG(geneList = ranks, organism = ORG, keyType = KT,
                          minGSSize = CFG$gsea_minGS, maxGSSize = CFG$gsea_maxGS,
                          pvalueCutoff = 1, pAdjustMethod = "BH", eps = 0, verbose = FALSE),
                  error = function(e) { message("  [error] GSEA: ", conditionMessage(e)); NULL })
  if (!is.null(gse) && nrow(as.data.frame(gse)) > 0) {
    g <- as.data.frame(gse)
    n_sig <- sum(g$p.adjust <= PADJ)
    message(sprintf("  KEGG GSEA %s : %d pathways tested, %d pass FDR (top nominal: %s, NES=%.2f, p=%.3g)",
                    cn, nrow(g), n_sig, g$Description[which.min(g$pvalue)],
                    g$NES[which.min(g$pvalue)], min(g$pvalue)))
    # full table: read by 02b_kegg_bubbleplots.R
    fwrite(g, file.path(CFG$dir_results, sprintf("prot__enrich__%s__kegg_gsea_all.tsv", cn)), sep = "\t")
    if (n_sig) {
      fwrite(g[g$p.adjust <= PADJ, ],
             file.path(CFG$dir_results, sprintf("prot__enrich__%s__kegg_gsea.tsv", cn)), sep = "\t")
      p <- enrichplot::dotplot(gse, showCategory = 20, split = ".sign") +
             facet_grid(. ~ .sign) + ggtitle(sprintf("KEGG GSEA - %s", cn)) + theme_mo(8)
      save_fig(p, "enrich", cn, "kegg_gsea", CFG,
               width_mm = CFG$width_double_mm, height_mm = CFG$height_mm * 1.3)
    }                                       # else: table only (no barplot)
  } else message("  [empty] KEGG GSEA ", cn)
}

writeLines(c("step\tkegg_enrichment", sprintf("time\t%s", Sys.time()),
  sprintf("organism\t%s\tkeytype\t%s", ORG, KT),
  sprintf("padj\t%.3g\tminGSSize\t%d\tmaxGSSize\t%d", PADJ, CFG$gsea_minGS, CFG$gsea_maxGS),
  "de_set\tlimma_significant up/down"),
  file.path(CFG$dir_logs, "02_kegg_enrichment.log"))
message("\nDONE: KEGG ORA + GSEA")
