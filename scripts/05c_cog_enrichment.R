#!/usr/bin/env Rscript
# =============================================================================
# 05c_cog_enrichment.R
# COG functional-category over-representation (ORA) for A. fumigatus Af293 on
# the limma DE. Sibling of 02 (KEGG) and 05 (GO): same DE sets, universe
# (detected proteome) and thresholds.
#
# Annotation: NCBI COG covers prokaryotes only, so the COG letters come from the
# eukaryotic extension (KOG / eggNOG), versioned and BLAST-free:
#   UniProt reference proteome (CFG$go_proteome) -> eggNOG cross-refs (KOG#### or
#   ENOG50#####) -> eggNOG 5.0 OG category at the Eukaryota level
#   (CFG$cog_eggnog_level = 2759). A multi-letter OG (e.g. "IQ") counts toward
#   each letter (standard COG convention). Both files are cached in annotation/.
#
# ORA: up / down = limma-significant (padj <= 0.05) with log2FC > 0 / < 0.
# COG categories are deliberately broad, so there is no maxGSSize cap.
#
# Input  : results/prot__de__<contrast>__results_limma.tsv
# Outputs: results/prot__enrich__cog_map_summary.tsv
#          results/prot__enrich__<contrast>__cog_protein_map.tsv
#          results/prot__enrich__<contrast>_<up|down>__cog_ora{_all,}.tsv
#          results/prot__enrich__<contrast>__cog_profile.tsv
#            (% of COG-annotated proteins per category: detected vs DEPs up / down)
#          figures/prot__enrich__<contrast>__cog_bubble_ora   ORA bubble, category x direction
#          figures/prot__enrich__<contrast>__cog_lfc_strip    log2FC per category (beeswarm)
#          logs/05c_cog_enrichment.log
# Run    : Rscript scripts/05c_cog_enrichment.R   (internet once, for the caches)
#
# Steps
#   0. Setup
#   1. Protein -> COG letter map
#   2. Enrichment per contrast
#      2.1 DE sets
#      2.2 Per-protein COG table
#      2.3 ORA (up / down)
#      2.4 Category profile table
#      2.5 ORA bubble
#      2.6 log2FC strip
#   3. Run log
# =============================================================================


# 0. Setup ---------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(ggplot2); library(ggbeeswarm)
})
options(timeout = 600)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
SDIR <- if (length(script_file)) dirname(normalizePath(script_file)) else "scripts"
source(file.path(SDIR, "config.R"))
source(file.path(SDIR, "io.R"))
source(file.path(SDIR, "plot_utils.R"))
source(file.path(SDIR, "enrich_utils.R"))
PADJ <- CFG$padj
LVL  <- CFG$cog_eggnog_level
T2N  <- COG_CATEGORIES[, .(term = letter, name)]


# 1. Protein -> COG letter map -------------------------------------------------
## 1.1 Download (once) ----
xref_cache <- file.path(CFG$dir_annotation, sprintf("uniprot_eggnog_%s.tsv", CFG$go_proteome))
og_cache   <- file.path(CFG$dir_annotation, sprintf("eggnog5_%d_annotations.tsv.gz", LVL))
fetch_cached(sprintf(paste0("https://rest.uniprot.org/uniprotkb/stream?",
                            "query=proteome:%s&fields=accession,xref_eggnog&format=tsv"),
                     CFG$go_proteome),
             xref_cache, CFG$cog_force_download)
fetch_cached(sprintf("http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/%d/%d_annotations.tsv.gz",
                     LVL, LVL),
             og_cache, CFG$cog_force_download)

## 1.2 UniProt accession -> eggNOG OG -> COG letters ----
xr <- fread(xref_cache, sep = "\t", quote = "", colClasses = "character")
setnames(xr, c("accession", "eggnog"))
xl <- xr[nzchar(eggnog), .(og = unlist(strsplit(eggnog, ";", fixed = TRUE))), by = accession]
xl <- unique(xl[nzchar(og)])
og <- fread(cmd = paste("gzip -dc", shQuote(og_cache)), header = FALSE, sep = "\t", quote = "",
            colClasses = "character", select = 2:3, col.names = c("og", "cats"))
og[!grepl("^KOG", og), og := paste0("ENOG50", og)]     # eggNOG 5 stores ENOG50 ids unprefixed
map <- merge(xl, og, by = "og")
C2G <- unique(map[, .(term = unlist(strsplit(cats, "", fixed = TRUE))), by = .(gene = accession)])
C2G <- C2G[term %chin% COG_CATEGORIES$letter, .(term, gene)]
message(sprintf("COG map: %d proteins with eggNOG xref, %d with a COG letter, %d letters used",
                uniqueN(xl$accession), uniqueN(C2G$gene), uniqueN(C2G$term)))

## 1.3 Map summary ----
fwrite(C2G[, .(n_prot = uniqueN(gene)), by = .(letter = term)][
         COG_CATEGORIES, on = "letter",
         .(letter, name, group = gsub("\n", " ", group), n_prot = fifelse(is.na(n_prot), 0L, n_prot))],
       file.path(CFG$dir_results, "prot__enrich__cog_map_summary.tsv"), sep = "\t")


# 2. Enrichment per contrast ---------------------------------------------------
for (ct in CFG$contrasts) {
  cn <- ct$name
  tg <- ct$test
  dl <- dir_labels(tg)
  dl_low <- dir_labels(tg, lower = TRUE)
  message("\n=== COG enrichment: ", cn, " ===")

  ## 2.1 DE sets ----
  d     <- read_de_results(CFG, cn)
  sets  <- enrichment_sets(d, PADJ)
  annot <- intersect(sets$universe, C2G$gene)               # detected proteins with a COG letter
  message(sprintf("  DE set (limma): up=%d down=%d | universe=%d (%d with a COG letter, %.0f%%)",
                  length(sets$up), length(sets$down), length(sets$universe),
                  length(annot), 100 * length(annot) / length(sets$universe)))

  ## 2.2 Per-protein COG table (detected proteome, Excel-ready) ----
  pmap <- C2G[gene %chin% sets$universe,
              .(cog = paste(sort(unique(term)), collapse = "")), by = .(protein_id = gene)]
  pmap <- merge(d[, .(protein_id, gene_raw, description, log2FC, padj, regulation)], pmap,
                by = "protein_id", all.x = TRUE)
  fwrite(pmap[order(padj)],
         file.path(CFG$dir_results, sprintf("prot__enrich__%s__cog_protein_map.tsv", cn)),
         sep = "\t", na = "NA")

  ## 2.3 ORA (up / down) ----
  ora_tabs <- list()
  for (dir in c("up", "down")) {
    r <- run_ora("COG", sprintf("%s/%s", cn, dir), sets[[dir]], padj_cut = PADJ,
                 enrich_fn = function(genes)
                   enricher(gene = genes, universe = sets$universe,
                            TERM2GENE = C2G, TERM2NAME = T2N,
                            pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
                            minGSSize = CFG$gsea_minGS, maxGSSize = length(sets$universe)))
    if (!has_terms(r)) next
    rdf <- as.data.table(as.data.frame(r))
    write_enrichment_tables(rdf, sprintf("%s_%s", cn, dir), "cog_ora", CFG)
    ora_tabs[[dir]] <- rdf
  }

  ## 2.4 Category profile table (% of COG-annotated proteins per set) ----
  profile_sets <- setNames(list(annot, intersect(sets$up, annot), intersect(sets$down, annot)),
                           c("Detected proteome", unname(dl)))
  prof <- rbindlist(lapply(names(profile_sets), function(s) {
    members <- profile_sets[[s]]
    C2G[gene %chin% members, .(n = uniqueN(gene)), by = .(letter = term)][
      , `:=`(set = s, n_set = length(members), pct = 100 * n / length(members))]
  }))
  fwrite(dcast(prof, letter ~ set, value.var = c("n", "pct"))[
           COG_CATEGORIES[, .(letter, name)], on = "letter", nomatch = 0],
         file.path(CFG$dir_results, sprintf("prot__enrich__%s__cog_profile.tsv", cn)), sep = "\t")

  # categories testable in the detected proteome (>= minGSSize, as in ORA)
  present <- prof[set == "Detected proteome" & n >= CFG$gsea_minGS, letter]
  fig_h   <- 30 + 6.2 * length(present)

  ## 2.5 ORA bubble ----
  dir_lab <- setNames(unname(dl), c("up", "down"))
  ob <- rbindlist(lapply(names(dir_lab), function(dir)
    if (!is.null(ora_tabs[[dir]]))
      ora_tabs[[dir]][, .(letter = ID, FoldEnrichment, Count, pvalue, p.adjust, direction = dir_lab[[dir]])]))
  if (nrow(ob)) {
    p <- plot_cog_ora_bubble(
      ob, tg, PADJ,
      title = sprintf("COG over-representation - %s", cn),
      subtitle = sprintf(paste0("categories with >=1 DEP; universe = detected proteome with a COG letter (n=%d)\n",
                                "fill = -log10 nominal p; black outline = BH FDR<=0.05"), length(annot)))
    save_fig(p, "enrich", cn, "cog_bubble_ora", CFG, width_mm = CFG$width_double_mm, height_mm = fig_h)
  }

  ## 2.6 log2FC strip ----
  st <- merge(C2G[gene %chin% sets$universe & term %chin% present, .(protein_id = gene, letter = term)],
              d[, .(protein_id, log2FC, padj, regulation)], by = "protein_id")
  p <- plot_cog_lfc_strip(
    st, tg, ct$control,
    title = sprintf("log2FC by COG category - %s", cn),
    subtitle = sprintf("each point = one detected protein; | = category median; log2FC > 0 = %s, < 0 = %s",
                       dl_low[["Up"]], dl_low[["Down"]]))
  save_fig(p, "enrich", cn, "cog_lfc_strip", CFG, width_mm = CFG$width_double_mm, height_mm = fig_h)
}


# 3. Run log -------------------------------------------------------------------
writeLines(c(
  "step\tcog_enrichment",
  sprintf("time\t%s", Sys.time()),
  sprintf("proteome\t%s\txref_cache\t%s", CFG$go_proteome, xref_cache),
  sprintf("eggnog\tv5.0 level %d\tcache\t%s", LVL, og_cache),
  sprintf("map\tproteins_with_cog\t%d\tletters\t%d", uniqueN(C2G$gene), uniqueN(C2G$term)),
  sprintf("padj\t%.3g\tminGSSize\t%d\tmaxGSSize\tnone", PADJ, CFG$gsea_minGS),
  "de_set\tlimma_significant up/down\tuniverse\tdetected proteome"),
  file.path(CFG$dir_logs, "05c_cog_enrichment.log"))
message("\nDONE: COG ORA + figures")
