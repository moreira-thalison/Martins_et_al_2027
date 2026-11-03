#!/usr/bin/env Rscript
# ============================================================================
# 05_go_enrichment.R  -  GO ORA + GSEA (BP / MF / CC) for the A. fumigatus
# (Af293) proteome, on the LIMMA DE (limma-only analysis). Sibling of
# 02_kegg_enrichment.R: same DE sets, same universe, same thresholds and output
# conventions, but the term->protein map is CURATED GO from the UniProt
# reference proteome (per aspect), so no org.Af.eg.db / GO.db is required
# (term names come from UniProt).
#
# Why curated UniProt GO (not Blast2GO): Af293 is a curated reference organism;
# UniProt/GOA already carries InterPro2GO evidence PLUS manual curation, is
# versioned/reproducible, and needs no BLAST step. (See project notes.)
#
# Differential gene set for ORA = limma-significant (padj <= 0.05):
#   up   = padj <= 0.05 & log2FC > 0   (higher in test)
#   down = padj <= 0.05 & log2FC < 0   (lower in test)
# GSEA ranks ALL tested proteins by limma log2FC.
# Universe = all tested proteins (the DETECTED proteome), intersected by
# enricher() with each aspect's annotated set.
# Thresholds: p.adjust<=0.05 (BH), minGSSize 5, maxGSSize 500; set.seed for GSEA.
# Requires internet ONCE to build the GO cache (then offline). IDs = UniProt.
#
# Run:  Rscript scripts/05_go_enrichment.R
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(enrichplot); library(ggplot2)
})
options(timeout = 600)
args <- commandArgs(FALSE); self <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
SDIR <- dirname(self); source(file.path(SDIR, "config.R")); source(file.path(SDIR, "plot_utils.R"))
PADJ <- CFG$padj

# ---- GO annotation: fetch (once) + parse to a long map ----------------------
# Cache = UniProt TSV with per-aspect GO columns ("name [GO:id]; ...").
go_cache <- file.path(CFG$dir_annotation,
                      sprintf("uniprot_go_%s.tsv", CFG$go_proteome))

fetch_go <- function(cache, proteome, force = FALSE) {
  if (!force && file.exists(cache) && file.info(cache)$size > 0) {
    message("  GO cache: ", cache); return(invisible(cache))
  }
  url <- sprintf(paste0("https://rest.uniprot.org/uniprotkb/stream?",
                        "query=proteome:%s&fields=accession,go_p,go_f,go_c&format=tsv"),
                 proteome)
  message("  downloading UniProt GO for ", proteome, " (release-versioned) ...")
  utils::download.file(url, cache, quiet = TRUE, mode = "wb")
  invisible(cache)
}

# Wide UniProt TSV -> long: one (accession, go_id, go_name, aspect) per row.
parse_go_long <- function(cache) {
  dt <- fread(cache, sep = "\t", quote = "", colClasses = "character")
  stopifnot(ncol(dt) == 4L)
  setnames(dt, c("Entry", "bp", "mf", "cc"))
  asp <- c(BP = "bp", MF = "mf", CC = "cc")
  out <- lapply(names(asp), function(a) {
    col <- asp[[a]]
    sub <- dt[nzchar(dt[[col]]), .(accession = Entry, cell = get(col))]
    if (!nrow(sub)) return(NULL)
    parts <- strsplit(sub$cell, "; ", fixed = TRUE)      # split multi-term cells
    long  <- data.table(accession = rep(sub$accession, lengths(parts)),
                        term_raw  = unlist(parts, use.names = FALSE))
    long[, go_id   := sub("^.*\\[(GO:[0-9]+)\\]$", "\\1", term_raw)]
    long[, go_name := sub("\\s*\\[GO:[0-9]+\\]$", "",     term_raw)]
    long[grepl("^GO:[0-9]+$", go_id),
         .(accession, go_id, go_name, aspect = a)]
  })
  rbindlist(out)
}

fetch_go(go_cache, CFG$go_proteome, force = CFG$go_force_download)
GO <- parse_go_long(go_cache)
T2N <- unique(GO[, .(term = go_id, name = go_name)])[!duplicated(term)]  # global names
message(sprintf("GO map: %d annotations, %d proteins, %d terms | BP %d / MF %d / CC %d terms",
                nrow(GO), uniqueN(GO$accession), uniqueN(GO$go_id),
                uniqueN(GO[aspect == "BP", go_id]), uniqueN(GO[aspect == "MF", go_id]),
                uniqueN(GO[aspect == "CC", go_id])))
fwrite(GO[, .(n_annot = .N, n_prot = uniqueN(accession), n_term = uniqueN(go_id)), by = aspect],
       file.path(CFG$dir_results, "prot__enrich__go_map_summary.tsv"), sep = "\t")

# Figures: dotplot only when terms pass FDR. No nominal-trend barplot (removed
# per user, 2026-10-27) - nominal trends are shown by the 05b bubble plots.

# ---- ORA over a custom (per-aspect) GO map ---------------------------------
ora_go <- function(genes, universe, t2g, label) {
  genes <- unique(genes[!is.na(genes)]); universe <- unique(universe[!is.na(universe)])
  if (length(genes) < 3) { message(sprintf("  [skip] %-40s <3 genes", label)); return(NULL) }
  r <- tryCatch(enricher(gene = genes, universe = universe,
                         TERM2GENE = t2g, TERM2NAME = T2N,
                         pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
                         minGSSize = CFG$gsea_minGS, maxGSSize = CFG$gsea_maxGS),
                error = function(e) { message("  [error] ORA ", label, ": ", conditionMessage(e)); NULL })
  n_sig <- if (is.null(r)) 0 else sum(as.data.frame(r)$p.adjust <= PADJ)
  message(sprintf("  GO ORA %-40s : %d terms tested, %d pass FDR",
                  label, if (is.null(r)) 0 else nrow(as.data.frame(r)), n_sig)); r
}

for (ct in CFG$contrasts) {
  cn <- ct$name; tg <- ct$test; cg <- ct$control
  message("\n=== GO enrichment: ", cn, " ===")
  # LIMMA-only: DE set = limma-significant up/down (no proDA / no on/off layer).
  d  <- fread(file.path(CFG$dir_results, sprintf("prot__de__%s__results_limma.tsv", cn)))

  universe <- unique(d$protein_id)                       # all tested proteins
  up   <- d[!is.na(padj) & padj <= PADJ & log2FC > 0, protein_id]
  down <- d[!is.na(padj) & padj <= PADJ & log2FC < 0, protein_id]
  n_univ_go <- length(intersect(universe, GO$accession))
  message(sprintf("  DE set (limma): up=%d down=%d | universe=%d (%d with any GO, %.0f%%)",
                  length(up), length(down), length(universe),
                  n_univ_go, 100 * n_univ_go / length(universe)))

  # ---- ORA: aspect (BP/MF/CC) x direction (up/down) ------------------------
  for (asp in CFG$go_aspects) {
    t2g <- GO[aspect == asp, .(term = go_id, gene = accession)]
    ak  <- tolower(asp)                                  # bp / mf / cc
    for (dir in c("up", "down")) {
      genes <- switch(dir, up = up, down = down)
      r <- ora_go(genes, universe, t2g, sprintf("%s/%s/%s", cn, asp, dir))
      if (is.null(r) || nrow(as.data.frame(r)) == 0) next
      rdf <- as.data.frame(r)
      fwrite(rdf, file.path(CFG$dir_results,                     # full table (read by 05b)
             sprintf("prot__enrich__%s_%s__go_%s_ora_all.tsv", cn, dir, ak)), sep = "\t")
      sig <- rdf[rdf$p.adjust <= PADJ, , drop = FALSE]
      if (nrow(sig)) {                                                     # FDR subset
        fwrite(sig, file.path(CFG$dir_results,
               sprintf("prot__enrich__%s_%s__go_%s_ora.tsv", cn, dir, ak)), sep = "\t")
        p <- enrichplot::dotplot(r, showCategory = 20) +
               ggtitle(sprintf("GO %s ORA - %s (%s)", asp, cn, dir)) + theme_mo(8)
        save_fig(p, "enrich", sprintf("%s_%s", cn, dir), sprintf("go_%s_ora", ak), CFG,
                 width_mm = CFG$width_double_mm, height_mm = CFG$height_mm * 1.3)
      }                                                    # else: table only
    }
  }

  # ---- GSEA on ranked limma log2FC, per aspect ----------------------------
  dd <- d[is.finite(log2FC) & !is.na(protein_id)]
  ranks <- dd$log2FC; names(ranks) <- dd$protein_id
  ranks <- sort(ranks[!duplicated(names(ranks))], decreasing = TRUE)
  for (asp in CFG$go_aspects) {
    t2g <- GO[aspect == asp, .(term = go_id, gene = accession)]
    ak  <- tolower(asp)
    set.seed(CFG$seed_gsea)
    gse <- tryCatch(clusterProfiler::GSEA(geneList = ranks, TERM2GENE = t2g, TERM2NAME = T2N,
                          minGSSize = CFG$gsea_minGS, maxGSSize = CFG$gsea_maxGS,
                          pvalueCutoff = 1, pAdjustMethod = "BH", eps = 0,
                          seed = TRUE, verbose = FALSE),
                    error = function(e) { message("  [error] GSEA ", asp, ": ", conditionMessage(e)); NULL })
    if (is.null(gse) || nrow(as.data.frame(gse)) == 0) { message("  [empty] GO ", asp, " GSEA ", cn); next }
    g <- as.data.frame(gse); n_sig <- sum(g$p.adjust <= PADJ)
    message(sprintf("  GO %s GSEA %s : %d sets tested, %d pass FDR (top nominal: %s, NES=%.2f, p=%.3g)",
                    asp, cn, nrow(g), n_sig, g$Description[which.min(g$pvalue)],
                    g$NES[which.min(g$pvalue)], min(g$pvalue)))
    fwrite(g, file.path(CFG$dir_results,
           sprintf("prot__enrich__%s__go_%s_gsea_all.tsv", cn, ak)), sep = "\t")
    if (n_sig) {
      fwrite(g[g$p.adjust <= PADJ, ], file.path(CFG$dir_results,
             sprintf("prot__enrich__%s__go_%s_gsea.tsv", cn, ak)), sep = "\t")
      p <- enrichplot::dotplot(gse, showCategory = 20, split = ".sign") +
             facet_grid(. ~ .sign) + ggtitle(sprintf("GO %s GSEA - %s", asp, cn)) + theme_mo(8)
      save_fig(p, "enrich", cn, sprintf("go_%s_gsea", ak), CFG,
               width_mm = CFG$width_double_mm, height_mm = CFG$height_mm * 1.3)
    }                                                      # else: table only
  }
}

writeLines(c("step\tgo_enrichment",
  sprintf("time\t%s", Sys.time()),
  sprintf("proteome\t%s\tcache\t%s", CFG$go_proteome, go_cache),
  sprintf("aspects\t%s", paste(CFG$go_aspects, collapse = ",")),
  sprintf("annotations\t%d\tproteins\t%d\tterms\t%d",
          nrow(GO), uniqueN(GO$accession), uniqueN(GO$go_id)),
  sprintf("padj\t%.3g\tminGSSize\t%d\tmaxGSSize\t%d", PADJ, CFG$gsea_minGS, CFG$gsea_maxGS),
  "de_set\tlimma_significant up/down\tuniverse\tdetected proteome"),
  file.path(CFG$dir_logs, "05_go_enrichment.log"))
message("\nDONE: GO ORA + GSEA (BP/MF/CC)")
