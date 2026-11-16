#!/usr/bin/env Rscript
# =============================================================================
# 05_go_enrichment.R
# GO ORA + GSEA (BP / MF / CC) for A. fumigatus Af293 on the limma DE.
# Sibling of 02_kegg_enrichment.R: same DE sets, universe, thresholds and file
# naming. The term -> protein map is curated GO from the UniProt reference
# proteome, so no org.Af.eg.db / GO.db is needed (term names come from UniProt).
#
# Why curated UniProt GO rather than Blast2GO: Af293 is a curated reference
# organism; UniProt/GOA already carries InterPro2GO evidence plus manual
# curation, is versioned and needs no BLAST step.
#
#   ORA  : up / down = limma-significant (padj <= 0.05) with log2FC > 0 / < 0
#          universe = every tested protein, intersected by enricher() with
#          each aspect's annotated set
#   GSEA : every tested protein ranked by limma log2FC
# BH p.adjust <= 0.05, minGSSize 5, maxGSSize 500. Internet needed once, to
# build the GO cache in annotation/.
#
# Input  : results/prot__de__<contrast>__results_limma.tsv
# Outputs: results/prot__enrich__go_map_summary.tsv
#          results/prot__enrich__<contrast>_<up|down>__go_<bp|mf|cc>_ora{_all,}.tsv
#          results/prot__enrich__<contrast>__go_<bp|mf|cc>_gsea{_all,}.tsv
#          dotplots (only when a term passes FDR), logs/05_go_enrichment.log
# Run    : Rscript scripts/05_go_enrichment.R
#
# Steps
#   0. Setup
#   1. GO annotation map
#   2. Enrichment per contrast and aspect
#      2.1 DE sets
#      2.2 ORA (up / down)
#      2.3 GSEA
#   3. Run log
# =============================================================================


# 0. Setup ---------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(enrichplot); library(ggplot2)
})
options(timeout = 600)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
SDIR <- if (length(script_file)) dirname(normalizePath(script_file)) else "scripts"
source(file.path(SDIR, "config.R"))
source(file.path(SDIR, "io.R"))
source(file.path(SDIR, "plot_utils.R"))
source(file.path(SDIR, "enrich_utils.R"))
PADJ <- CFG$padj


# 1. GO annotation map ---------------------------------------------------------
# Cache = UniProt TSV with one GO column per aspect ("name [GO:id]; ...").
go_cache <- file.path(CFG$dir_annotation, sprintf("uniprot_go_%s.tsv", CFG$go_proteome))
fetch_cached(sprintf(paste0("https://rest.uniprot.org/uniprotkb/stream?",
                            "query=proteome:%s&fields=accession,go_p,go_f,go_c&format=tsv"),
                     CFG$go_proteome),
             go_cache, force = CFG$go_force_download)

# Wide UniProt TSV -> long: one (accession, go_id, go_name, aspect) per row.
parse_go_long <- function(cache) {
  dt <- fread(cache, sep = "\t", quote = "", colClasses = "character")
  stopifnot(ncol(dt) == 4L)
  setnames(dt, c("Entry", "bp", "mf", "cc"))
  aspect_cols <- c(BP = "bp", MF = "mf", CC = "cc")
  rbindlist(lapply(names(aspect_cols), function(a) {
    col <- aspect_cols[[a]]
    sub <- dt[nzchar(dt[[col]]), .(accession = Entry, cell = get(col))]
    if (!nrow(sub)) return(NULL)
    parts <- strsplit(sub$cell, "; ", fixed = TRUE)        # split multi-term cells
    long  <- data.table(accession = rep(sub$accession, lengths(parts)),
                        term_raw  = unlist(parts, use.names = FALSE))
    long[, go_id   := sub("^.*\\[(GO:[0-9]+)\\]$", "\\1", term_raw)]
    long[, go_name := sub("\\s*\\[GO:[0-9]+\\]$", "", term_raw)]
    long[grepl("^GO:[0-9]+$", go_id), .(accession, go_id, go_name, aspect = a)]
  }))
}

GO  <- parse_go_long(go_cache)
T2N <- unique(GO[, .(term = go_id, name = go_name)])[!duplicated(term)]   # term names
message(sprintf("GO map: %d annotations, %d proteins, %d terms | BP %d / MF %d / CC %d terms",
                nrow(GO), uniqueN(GO$accession), uniqueN(GO$go_id),
                uniqueN(GO[aspect == "BP", go_id]), uniqueN(GO[aspect == "MF", go_id]),
                uniqueN(GO[aspect == "CC", go_id])))
fwrite(GO[, .(n_annot = .N, n_prot = uniqueN(accession), n_term = uniqueN(go_id)), by = aspect],
       file.path(CFG$dir_results, "prot__enrich__go_map_summary.tsv"), sep = "\t")


# 2. Enrichment per contrast and aspect ----------------------------------------
for (ct in CFG$contrasts) {
  cn <- ct$name
  message("\n=== GO enrichment: ", cn, " ===")

  ## 2.1 DE sets ----
  d     <- read_de_results(CFG, cn)
  sets  <- enrichment_sets(d, PADJ)
  ranks <- ranked_log2fc(d)
  n_univ_go <- length(intersect(sets$universe, GO$accession))
  message(sprintf("  DE set (limma): up=%d down=%d | universe=%d (%d with any GO, %.0f%%)",
                  length(sets$up), length(sets$down), length(sets$universe),
                  n_univ_go, 100 * n_univ_go / length(sets$universe)))

  for (asp in CFG$go_aspects) {
    t2g <- GO[aspect == asp, .(term = go_id, gene = accession)]
    ak  <- tolower(asp)                                     # bp / mf / cc

    ## 2.2 ORA (up / down) ----
    for (dir in c("up", "down")) {
      r <- run_ora("GO", sprintf("%s/%s/%s", cn, asp, dir), sets[[dir]], padj_cut = PADJ,
                   enrich_fn = function(genes)
                     enricher(gene = genes, universe = sets$universe,
                              TERM2GENE = t2g, TERM2NAME = T2N,
                              pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
                              minGSSize = CFG$gsea_minGS, maxGSSize = CFG$gsea_maxGS))
      if (!has_terms(r)) next
      save_enrichment(r, sprintf("%s_%s", cn, dir), sprintf("go_%s_ora", ak),
                      sprintf("GO %s ORA - %s (%s)", asp, cn, dir), CFG)
    }

    ## 2.3 GSEA ----
    set.seed(CFG$seed_gsea)
    gse <- run_gsea(sprintf("GO %s", asp), cn, padj_cut = PADJ, gsea_fn = function()
      clusterProfiler::GSEA(geneList = ranks, TERM2GENE = t2g, TERM2NAME = T2N,
                            minGSSize = CFG$gsea_minGS, maxGSSize = CFG$gsea_maxGS,
                            pvalueCutoff = 1, pAdjustMethod = "BH", eps = 0,
                            seed = TRUE, verbose = FALSE))
    if (!is.null(gse))
      save_enrichment(gse, cn, sprintf("go_%s_gsea", ak),
                      sprintf("GO %s GSEA - %s", asp, cn), CFG, gsea = TRUE)
  }
}


# 3. Run log -------------------------------------------------------------------
writeLines(c(
  "step\tgo_enrichment",
  sprintf("time\t%s", Sys.time()),
  sprintf("proteome\t%s\tcache\t%s", CFG$go_proteome, go_cache),
  sprintf("aspects\t%s", paste(CFG$go_aspects, collapse = ",")),
  sprintf("annotations\t%d\tproteins\t%d\tterms\t%d",
          nrow(GO), uniqueN(GO$accession), uniqueN(GO$go_id)),
  sprintf("padj\t%.3g\tminGSSize\t%d\tmaxGSSize\t%d", PADJ, CFG$gsea_minGS, CFG$gsea_maxGS),
  "de_set\tlimma_significant up/down\tuniverse\tdetected proteome"),
  file.path(CFG$dir_logs, "05_go_enrichment.log"))
message("\nDONE: GO ORA + GSEA (BP/MF/CC)")
