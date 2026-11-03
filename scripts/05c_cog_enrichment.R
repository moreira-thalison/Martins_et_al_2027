#!/usr/bin/env Rscript
# ============================================================================
# 05c_cog_enrichment.R  -  COG functional-category enrichment (ORA) for
# the A. fumigatus (Af293) proteome, on the LIMMA DE. Sibling of 02 (KEGG) and
# 05 (GO): same DE sets, same universe (detected proteome), same thresholds.
#
# Annotation: NCBI COG covers prokaryotes only, so for this fungus the COG letters
# come from the eukaryotic extension (KOG / eggNOG). Route, fully versioned and
# BLAST-free (same philosophy as the curated UniProt GO in 05):
#   UniProt reference proteome (CFG$go_proteome) -> eggNOG cross-refs (KOG#### or
#   ENOG50#####)  ->  eggNOG 5.0 OG functional category at the Eukaryota level
#   (CFG$cog_eggnog_level = 2759). Multi-letter OGs (e.g. "IQ") count toward
#   EACH letter (standard COG convention). Both files cached in annotation/.
#
# DE set for ORA = limma-significant (padj <= 0.05): up = log2FC > 0 (higher in
# test), down = log2FC < 0 (lower in test). COG categories are deliberately
# broad, so no maxGSSize cap.
#
# Tables: cog_protein_map, cog_ora_all (up / down), cog_profile (% of annotated
#         proteins per category in the detected proteome vs DEPs up / down)
# Figures (no bar charts):
#   cog_bubble_ora   - ORA bubble, category x direction (fold enrichment)
#   cog_lfc_strip    - log2FC of every detected protein per category (beeswarm)
#
# Run:  Rscript scripts/05c_cog_enrichment.R      (internet once, for the caches)
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(ggplot2); library(ggbeeswarm)
})
options(timeout = 600)
args <- commandArgs(FALSE); self <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
SDIR <- dirname(self); source(file.path(SDIR, "config.R")); source(file.path(SDIR, "plot_utils.R"))
PADJ <- CFG$padj
LVL  <- CFG$cog_eggnog_level

# ---- COG functional categories (letter, short name, supergroup) -------------
GROUPS <- c("Information storage\n& processing", "Cellular processes\n& signaling",
            "Metabolism", "Poorly\ncharacterized")
COG <- data.table(
  letter = c("J", "A", "K", "L", "B",
             "D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O",
             "C", "G", "E", "F", "H", "I", "P", "Q",
             "R", "S"),
  name = c("Translation, ribosomal structure & biogenesis", "RNA processing & modification",
           "Transcription", "Replication, recombination & repair",
           "Chromatin structure & dynamics",
           "Cell cycle control & division", "Nuclear structure", "Defense mechanisms",
           "Signal transduction", "Cell wall/membrane/envelope biogenesis", "Cell motility",
           "Cytoskeleton", "Extracellular structures", "Intracellular trafficking & secretion",
           "PTM, protein turnover, chaperones",
           "Energy production & conversion", "Carbohydrate transport & metabolism",
           "Amino acid transport & metabolism", "Nucleotide transport & metabolism",
           "Coenzyme transport & metabolism", "Lipid transport & metabolism",
           "Inorganic ion transport & metabolism", "Secondary metabolite biosynthesis",
           "General function prediction only", "Function unknown"),
  group = rep(GROUPS, c(5, 10, 8, 2))
)
COG[, group := factor(group, levels = GROUPS)]
COG[, lab := sprintf("%s  %s", letter, name)]
T2N <- COG[, .(term = letter, name)]

# ---- annotation: fetch (once) + build protein -> COG letter map -------------
xref_cache <- file.path(CFG$dir_annotation, sprintf("uniprot_eggnog_%s.tsv", CFG$go_proteome))
og_cache   <- file.path(CFG$dir_annotation, sprintf("eggnog5_%d_annotations.tsv.gz", LVL))
fetch <- function(url, cache, force = FALSE) {
  if (!force && file.exists(cache) && file.info(cache)$size > 0) {
    message("  cache: ", cache); return(invisible(cache))
  }
  message("  downloading ", url)
  utils::download.file(url, cache, quiet = TRUE, mode = "wb")
  invisible(cache)
}
fetch(sprintf(paste0("https://rest.uniprot.org/uniprotkb/stream?",
                     "query=proteome:%s&fields=accession,xref_eggnog&format=tsv"),
              CFG$go_proteome), xref_cache, CFG$cog_force_download)
fetch(sprintf("http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/%d/%d_annotations.tsv.gz",
              LVL, LVL), og_cache, CFG$cog_force_download)

xr <- fread(xref_cache, sep = "\t", quote = "", colClasses = "character")
setnames(xr, c("accession", "eggnog"))
xl <- xr[nzchar(eggnog), .(og = unlist(strsplit(eggnog, ";", fixed = TRUE))), by = accession]
xl <- unique(xl[nzchar(og)])
og <- fread(cmd = paste("gzip -dc", shQuote(og_cache)), header = FALSE, sep = "\t", quote = "",
            colClasses = "character", select = 2:3, col.names = c("og", "cats"))
og[!grepl("^KOG", og), og := paste0("ENOG50", og)]     # eggNOG 5 stores ENOG50 ids unprefixed
map <- merge(xl, og, by = "og")
C2G <- unique(map[, .(term = unlist(strsplit(cats, "", fixed = TRUE))), by = .(gene = accession)])
C2G <- C2G[term %chin% COG$letter, .(term, gene)]
message(sprintf("COG map: %d proteins with eggNOG xref, %d with a COG letter, %d letters used",
                uniqueN(xl$accession), uniqueN(C2G$gene), uniqueN(C2G$term)))
fwrite(C2G[, .(n_prot = uniqueN(gene)), by = .(letter = term)][COG, on = "letter",
           .(letter, name, group = gsub("\n", " ", group), n_prot = fifelse(is.na(n_prot), 0L, n_prot))],
       file.path(CFG$dir_results, "prot__enrich__cog_map_summary.tsv"), sep = "\t")

# ---- shared plot pieces (fill = -log10 nominal p, outline = FDR pass) -------
bubble_layers <- function(size_name) {
  list(
    # show.legend = TRUE: draw the FDR key even when no category passes (ggplot2 >= 3.5)
    geom_point(shape = 21, stroke = 0.5, alpha = 0.9, show.legend = TRUE),
    scale_fill_gradient(low = "grey85", high = "#d73027",
                        name = expression(-log[10]~italic(p))),
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey70"),
                        labels = c(`TRUE` = "FDR<=0.05", `FALSE` = "n.s."),
                        name = "BH", drop = FALSE),
    scale_size_continuous(range = c(2, 8), name = size_name),
    guides(colour = guide_legend(override.aes = list(size = 3.5, fill = "white", alpha = 1)),
           fill = guide_colourbar(raster = FALSE)),
    cog_theme()
  )
}
cog_theme <- function() {
  list(theme_mo(8),
       theme(plot.title.position = "plot",
             panel.grid.major.y = element_line(colour = "grey92"),
             strip.text = element_text(face = "bold"),
             strip.text.y = element_text(angle = 0, hjust = 0)))
}
# category axis: COG order top -> bottom inside each supergroup facet
cog_axis <- function(dt) {
  dt <- merge(dt, COG[, .(letter, lab, group)], by = "letter")
  dt[, lab := factor(lab, levels = rev(COG$lab))]
  dt
}

for (ct in CFG$contrasts) {
  cn <- ct$name; tg <- ct$test; cg <- ct$control
  message("\n=== COG enrichment: ", cn, " ===")
  d <- fread(file.path(CFG$dir_results, sprintf("prot__de__%s__results_limma.tsv", cn)))

  universe <- unique(d$protein_id)                       # all tested proteins
  up   <- d[!is.na(padj) & padj <= PADJ & log2FC > 0, protein_id]
  down <- d[!is.na(padj) & padj <= PADJ & log2FC < 0, protein_id]
  annot <- intersect(universe, C2G$gene)
  message(sprintf("  DE set (limma): up=%d down=%d | universe=%d (%d with a COG letter, %.0f%%)",
                  length(up), length(down), length(universe),
                  length(annot), 100 * length(annot) / length(universe)))

  # per-protein COG table for the detected proteome (Excel-ready)
  pm <- C2G[gene %chin% universe, .(cog = paste(sort(unique(term)), collapse = "")), by = .(protein_id = gene)]
  pm <- merge(d[, .(protein_id, gene_raw, description, log2FC, padj, regulation)], pm,
              by = "protein_id", all.x = TRUE)
  fwrite(pm[order(padj)], file.path(CFG$dir_results, sprintf("prot__enrich__%s__cog_protein_map.tsv", cn)),
         sep = "\t", na = "NA")

  # ---- ORA: up / down -------------------------------------------------------
  ora_tabs <- list()
  for (dir in c("up", "down")) {
    genes <- switch(dir, up = up, down = down)
    if (length(genes) < 3) { message(sprintf("  [skip] COG ORA %s <3 genes", dir)); next }
    r <- tryCatch(enricher(gene = genes, universe = universe, TERM2GENE = C2G, TERM2NAME = T2N,
                           pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
                           minGSSize = CFG$gsea_minGS, maxGSSize = length(universe)),
                  error = function(e) { message("  [error] COG ORA ", dir, ": ", conditionMessage(e)); NULL })
    if (is.null(r) || nrow(as.data.frame(r)) == 0) next
    rdf <- as.data.table(as.data.frame(r))
    message(sprintf("  COG ORA %-5s : %d categories tested, %d pass FDR (top nominal: %s %s, p=%.3g)",
                    dir, nrow(rdf), sum(rdf$p.adjust <= PADJ), rdf$ID[which.min(rdf$pvalue)],
                    rdf$Description[which.min(rdf$pvalue)], min(rdf$pvalue)))
    fwrite(rdf, file.path(CFG$dir_results, sprintf("prot__enrich__%s_%s__cog_ora_all.tsv", cn, dir)), sep = "\t")
    if (any(rdf$p.adjust <= PADJ))
      fwrite(rdf[p.adjust <= PADJ], file.path(CFG$dir_results,
             sprintf("prot__enrich__%s_%s__cog_ora.tsv", cn, dir)), sep = "\t")
    ora_tabs[[dir]] <- rdf
  }

  # ---- category profile table (% of COG-annotated proteins per set) ---------
  sets <- list(`Detected proteome` = annot,
               up   = intersect(up, annot),
               down = intersect(down, annot))
  names(sets)[2:3] <- unname(dir_labels(tg))
  prof <- rbindlist(lapply(names(sets), function(s) {
    C2G[gene %chin% sets[[s]], .(n = uniqueN(gene)), by = .(letter = term)][
      , `:=`(set = s, n_set = length(sets[[s]]), pct = 100 * n / length(sets[[s]]))]
  }))
  fwrite(dcast(prof, letter ~ set, value.var = c("n", "pct"))[COG[, .(letter, name)], on = "letter", nomatch = 0],
         file.path(CFG$dir_results, sprintf("prot__enrich__%s__cog_profile.tsv", cn)), sep = "\t")

  # categories testable in the detected proteome (>= minGSSize, as in ORA)
  present <- prof[set == "Detected proteome" & n >= CFG$gsea_minGS, letter]
  fig_h <- 30 + 6.2 * length(present)

  # ---- (1) ORA bubble: category x direction ---------------------------------
  dir_lab <- setNames(dir_labels(tg), c("up", "down"))
  ob <- rbindlist(lapply(names(dir_lab), function(dr)
    if (!is.null(ora_tabs[[dr]])) ora_tabs[[dr]][, .(letter = ID, FoldEnrichment, Count, pvalue,
                                                    p.adjust, direction = dir_lab[[dr]])]))
  if (nrow(ob)) {
    ob <- cog_axis(ob)
    ob[, direction := factor(direction, levels = dir_lab)]
    ob[, sig := factor(p.adjust <= PADJ, levels = c(TRUE, FALSE))]
    p <- ggplot(ob, aes(FoldEnrichment, lab, fill = -log10(pvalue), size = Count, colour = sig)) +
      geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
      bubble_layers("DEPs in\ncategory") +
      facet_grid(group ~ direction, scales = "free_y", space = "free_y") +
      labs(title = sprintf("COG over-representation - %s", cn),
           subtitle = sprintf("categories with >=1 DEP; universe = detected proteome with a COG letter (n=%d)\nfill = -log10 nominal p; black outline = BH FDR<=0.05",
                              length(annot)),
           x = "fold enrichment", y = NULL)
    save_fig(p, "enrich", cn, "cog_bubble_ora", CFG, width_mm = CFG$width_double_mm, height_mm = fig_h)
  }

  # ---- (2) log2FC of every detected protein per category (beeswarm) ---------
  st <- merge(C2G[gene %chin% universe & term %chin% present, .(protein_id = gene, letter = term)],
              d[, .(protein_id, log2FC, padj)], by = "protein_id")
  dl <- dir_labels(tg)
  st[, status := fifelse(!is.na(padj) & padj <= PADJ & log2FC > 0, dl[["Up"]],
                 fifelse(!is.na(padj) & padj <= PADJ & log2FC < 0, dl[["Down"]], "n.s."))]
  st[, status := factor(status, levels = c(unname(dl), "n.s."))]
  st <- cog_axis(st)
  nlab <- st[, .N, by = lab]
  st[, lab := factor(sprintf("%s (n=%d)", lab, nlab$N[match(lab, nlab$lab)]),
                     levels = sprintf("%s (n=%d)", levels(lab), nlab$N[match(levels(lab), nlab$lab)]))]
  med <- st[, .(med = median(log2FC)), by = .(lab, group)]
  st_cols <- setNames(c(REG_COLORS[["Up"]], REG_COLORS[["Down"]], "#bdbdbd"), levels(st$status))
  p <- ggplot(st, aes(log2FC, lab)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    geom_quasirandom(data = st[status == "n.s."], aes(colour = status), orientation = "y",
                     size = 0.45, alpha = 0.55, width = 0.38) +
    geom_quasirandom(data = st[status != "n.s."], aes(colour = status), orientation = "y",
                     size = 1.1, alpha = 0.95, width = 0.38) +
    geom_point(data = med, aes(med, lab), shape = 124, size = 5, colour = "grey10") +
    scale_colour_manual(values = st_cols, name = NULL, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    labs(title = sprintf("log2FC by COG category - %s", cn),
         subtitle = sprintf("each point = one detected protein; | = category median; log2FC > 0 = %s, < 0 = %s",
                            dir_labels(tg, lower = TRUE)[["Up"]], dir_labels(tg, lower = TRUE)[["Down"]]),
         x = sprintf("log2FC (%s / %s)", tg, cg), y = NULL) +
    cog_theme() + theme(legend.position = "top")
  save_fig(p, "enrich", cn, "cog_lfc_strip", CFG, width_mm = CFG$width_double_mm, height_mm = fig_h)
}

writeLines(c("step\tcog_enrichment",
  sprintf("time\t%s", Sys.time()),
  sprintf("proteome\t%s\txref_cache\t%s", CFG$go_proteome, xref_cache),
  sprintf("eggnog\tv5.0 level %d\tcache\t%s", LVL, og_cache),
  sprintf("map\tproteins_with_cog\t%d\tletters\t%d", uniqueN(C2G$gene), uniqueN(C2G$term)),
  sprintf("padj\t%.3g\tminGSSize\t%d\tmaxGSSize\tnone", PADJ, CFG$gsea_minGS),
  "de_set\tlimma_significant up/down\tuniverse\tdetected proteome"),
  file.path(CFG$dir_logs, "05c_cog_enrichment.log"))
message("\nDONE: COG ORA + figures")
