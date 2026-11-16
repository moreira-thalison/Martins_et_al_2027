#!/usr/bin/env Rscript
# =============================================================================
# enrich_utils.R
# Shared helpers for the enrichment steps (02 KEGG, 05 GO, 05c COG) and the
# bubble plots built from their tables (02b, 05b, 09b). Sourced after
# config.R and plot_utils.R.
#
# Sections
#   1. DE sets and rankings    enrichment_sets(), ranked_log2fc()
#   2. Annotation downloads    fetch_cached()
#   3. Running enrichment      run_ora(), run_gsea(), has_terms()
#   4. Writing results         write_enrichment_tables(), save_enrichment()
#   5. Bubble plots            read_result_table(), tidy_label(), order_ykey(),
#                              bubble_layers(), enrich_theme(), go_gsea_top()
#   6. COG categories          COG_CATEGORIES, cog_axis(), plot_cog_ora_bubble(),
#                              plot_cog_lfc_strip()
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })


# 1. DE sets and rankings ------------------------------------------------------
# universe = every tested protein (the detected proteome);
# up / down = limma-significant proteins higher / lower in the test group.
enrichment_sets <- function(d, padj_cut) {
  sig <- !is.na(d$padj) & d$padj <= padj_cut
  list(universe = unique(d$protein_id),
       up       = d$protein_id[which(sig & d$log2FC > 0)],
       down     = d$protein_id[which(sig & d$log2FC < 0)])
}

# GSEA input: every tested protein ranked by limma log2FC (decreasing).
ranked_log2fc <- function(d) {
  dd <- d[is.finite(log2FC) & !is.na(protein_id)]
  ranks <- setNames(dd$log2FC, dd$protein_id)
  sort(ranks[!duplicated(names(ranks))], decreasing = TRUE)
}


# 2. Annotation downloads ------------------------------------------------------
# Download url to cache once; later runs reuse the cache unless force = TRUE.
fetch_cached <- function(url, cache, force = FALSE) {
  if (!force && file.exists(cache) && file.info(cache)$size > 0) {
    message("  cache: ", cache)
    return(invisible(cache))
  }
  message("  downloading ", url)
  utils::download.file(url, cache, quiet = TRUE, mode = "wb")
  invisible(cache)
}


# 3. Running enrichment --------------------------------------------------------
# Every enrichment keeps the FULL table (pvalueCutoff = 1 in the callers); FDR
# filtering and reporting happen here and in section 4.

has_terms <- function(r) !is.null(r) && nrow(as.data.frame(r)) > 0

# Over-representation: enrich_fn(genes) runs the test. Fewer than 3 genes is
# skipped; an error is reported and returns NULL, so one failing set does not
# stop the step.
run_ora <- function(kind, label, genes, enrich_fn, padj_cut) {
  genes <- unique(genes[!is.na(genes)])
  if (length(genes) < 3) {
    message(sprintf("  [skip] %s ORA %s: <3 genes", kind, label))
    return(NULL)
  }
  r <- tryCatch(enrich_fn(genes), error = function(e) {
    message("  [error] ", kind, " ORA ", label, ": ", conditionMessage(e))
    NULL
  })
  tab <- if (is.null(r)) NULL else as.data.frame(r)
  top <- if (NROW(tab)) sprintf(" (top nominal: %s, p=%.3g)",
                                tab$Description[which.min(tab$pvalue)], min(tab$pvalue)) else ""
  message(sprintf("  %s ORA %-28s : %d terms tested, %d pass FDR%s",
                  kind, label, NROW(tab), sum(tab$p.adjust <= padj_cut), top))
  r
}

# Gene-set enrichment: gsea_fn() runs the test (set the seed just before).
run_gsea <- function(kind, label, gsea_fn, padj_cut) {
  r <- tryCatch(gsea_fn(), error = function(e) {
    message("  [error] ", kind, " GSEA ", label, ": ", conditionMessage(e))
    NULL
  })
  if (!has_terms(r)) {
    message("  [empty] ", kind, " GSEA ", label)
    return(NULL)
  }
  g <- as.data.frame(r)
  top <- which.min(g$pvalue)
  message(sprintf("  %s GSEA %s : %d sets tested, %d pass FDR (top nominal: %s, NES=%.2f, p=%.3g)",
                  kind, label, nrow(g), sum(g$p.adjust <= padj_cut),
                  g$Description[top], g$NES[top], g$pvalue[top]))
  r
}


# 4. Writing results -----------------------------------------------------------
# Writes results/prot__enrich__<key>__<name>_all.tsv (full table, read by the
# bubble-plot scripts) and, when any term passes FDR, <...>__<name>.tsv with
# only those terms. Returns the FDR-passing rows.
write_enrichment_tables <- function(tab, key, name, CFG) {
  tab  <- as.data.frame(tab)
  stem <- file.path(CFG$dir_results, sprintf("prot__enrich__%s__%s", key, name))
  fwrite(tab, paste0(stem, "_all.tsv"), sep = "\t")
  sig <- tab[tab$p.adjust <= CFG$padj, , drop = FALSE]
  if (nrow(sig)) fwrite(sig, paste0(stem, ".tsv"), sep = "\t")
  invisible(sig)
}

# Tables as above, plus a clusterProfiler dotplot when any term passes FDR
# (nominal trends are shown by the bubble plots instead; no bar charts).
save_enrichment <- function(r, key, name, title, CFG, gsea = FALSE) {
  sig <- write_enrichment_tables(r, key, name, CFG)
  if (!nrow(sig)) return(invisible(NULL))
  p <- if (gsea) {
    enrichplot::dotplot(r, showCategory = 20, split = ".sign") + facet_grid(. ~ .sign)
  } else {
    enrichplot::dotplot(r, showCategory = 20)
  }
  save_fig(p + ggtitle(title) + theme_mo(8), "enrich", key, name, CFG,
           width_mm = CFG$width_double_mm, height_mm = CFG$height_mm * 1.3)
}


# 5. Bubble plots --------------------------------------------------------------
# Shared encoding: x = fold enrichment or NES, size = overlap / set size,
# fill = -log10(nominal p), black outline = BH FDR <= 0.05.

GO_ASPECT_LABELS <- c(BP = "Biological process", MF = "Molecular function",
                      CC = "Cellular component")

# An enrichment table written to results/, or NULL if that step produced none.
read_result_table <- function(file, CFG) {
  p <- file.path(CFG$dir_results, file)
  if (file.exists(p)) fread(p) else NULL
}

# Cap label length, then wrap to lines of <= width characters.
tidy_label <- function(x, max = 80, width = 42) {
  x <- ifelse(nchar(x) > max, paste0(substr(x, 1, max - 3), "..."), x)
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
}

# Sort by ord_cols and give each row its own y position, so the order holds
# even when a term appears in several facets.
order_ykey <- function(dt, ord_cols) {
  setorderv(dt, ord_cols)
  dt[, ykey := factor(.I, levels = .I, labels = Description)]
  dt
}

bubble_layers <- function(size_name, size_range = c(2, 8), key_size = 3.5) {
  list(
    # show.legend = TRUE draws the FDR key even when no term passes
    geom_point(shape = 21, stroke = 0.5, alpha = 0.9, show.legend = TRUE),
    scale_fill_gradient(low = "grey85", high = "#d73027",
                        name = expression(-log[10]~italic(p))),
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey70"),
                        labels = c(`TRUE` = "FDR<=0.05", `FALSE` = "n.s."),
                        name = "BH", drop = FALSE),
    scale_size_continuous(range = size_range, name = size_name),
    # display = "rectangles": vector colourbar, so SVG / PDF editors get no bitmap
    guides(colour = guide_legend(override.aes = list(size = key_size, fill = "white", alpha = 1)),
           fill = guide_colourbar(display = "rectangles"))
  )
}

# Theme for KEGG / COG bubbles and the COG strip (horizontal row strips).
enrich_theme <- function() {
  list(theme_mo(8),
       theme(plot.title.position = "plot",
             panel.grid.major.y = element_line(colour = "grey92"),
             strip.text = element_text(face = "bold"),
             strip.text.y = element_text(angle = 0, hjust = 0)))
}

# Top n_per_aspect GO GSEA sets (by nominal p) per aspect, ready to plot.
go_gsea_top <- function(contrast, n_per_aspect, CFG) {
  gse <- rbindlist(lapply(names(GO_ASPECT_LABELS), function(asp) {
    x <- read_result_table(sprintf("prot__enrich__%s__go_%s_gsea_all.tsv", contrast, tolower(asp)), CFG)
    if (is.null(x) || !nrow(x)) return(NULL)
    x[, aspect := GO_ASPECT_LABELS[[asp]]]
    x[order(pvalue)][seq_len(min(n_per_aspect, .N))]
  }), fill = TRUE)
  if (!nrow(gse)) return(gse)
  gse[, aspect := factor(aspect, levels = GO_ASPECT_LABELS)]
  gse[, sig := p.adjust <= CFG$padj]
  gse[, neglogp := -log10(pvalue)]
  gse[, Description := tidy_label(Description)]
  order_ykey(gse, c("aspect", "NES"))
}


# 6. COG categories ------------------------------------------------------------
## 6.1 Category table (letter, short name, supergroup) ----
COG_GROUPS <- c("Information storage\n& processing", "Cellular processes\n& signaling",
                "Metabolism", "Poorly\ncharacterized")
COG_CATEGORIES <- data.table(
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
  group = factor(rep(COG_GROUPS, c(5, 10, 8, 2)), levels = COG_GROUPS)
)
COG_CATEGORIES[, lab := sprintf("%s  %s", letter, name)]

# Add the category label and supergroup; labels run in COG order top -> bottom.
cog_axis <- function(dt) {
  dt <- merge(dt, COG_CATEGORIES[, .(letter, lab, group)], by = "letter")
  dt[, lab := factor(lab, levels = rev(COG_CATEGORIES$lab))]
  dt
}

## 6.2 ORA bubble: category x direction ----
# ob: one row per category and direction (letter, FoldEnrichment, Count, pvalue,
# p.adjust, direction = "Higher in <test>" / "Lower in <test>").
plot_cog_ora_bubble <- function(ob, test, padj_cut, title, subtitle, size_range = c(2, 8)) {
  ob <- cog_axis(ob)
  ob[, direction := factor(direction, levels = dir_labels(test))]
  ob[, sig := factor(p.adjust <= padj_cut, levels = c(TRUE, FALSE))]
  ggplot(ob, aes(FoldEnrichment, lab, fill = -log10(pvalue), size = Count, colour = sig)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
    bubble_layers("DEPs in\ncategory", size_range = size_range) +
    facet_grid(group ~ direction, scales = "free_y", space = "free_y") +
    labs(title = title, subtitle = subtitle, x = "fold enrichment", y = NULL) +
    enrich_theme()
}

## 6.3 log2FC strip: one point per detected protein and category ----
# st: one row per protein and category letter (protein_id, letter, log2FC, padj,
# regulation). Rows keep their input order inside each category.
plot_cog_lfc_strip <- function(st, test, control, title, subtitle) {
  st <- copy(st)
  st[, status := direction_status(regulation, "n.s.", test)]
  st <- cog_axis(st)
  n_per_cat <- st[, .N, by = lab]
  with_n <- function(x) sprintf("%s (n=%d)", x, n_per_cat$N[match(x, n_per_cat$lab)])
  st[, lab := factor(with_n(lab), levels = with_n(levels(lab)))]
  med  <- st[, .(med = median(log2FC)), by = .(lab, group)]
  cols <- setNames(unname(REG_COLORS[c("Up", "Down", "NS")]), levels(st$status))

  ggplot(st, aes(log2FC, lab)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    ggbeeswarm::geom_quasirandom(data = st[status == "n.s."], aes(colour = status),
                                 orientation = "y", size = 0.45, alpha = 0.55, width = 0.38) +
    ggbeeswarm::geom_quasirandom(data = st[status != "n.s."], aes(colour = status),
                                 orientation = "y", size = 1.1, alpha = 0.95, width = 0.38) +
    geom_point(data = med, aes(med, lab), shape = 124, size = 5, colour = "grey10") +
    scale_colour_manual(values = cols, name = NULL, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    labs(title = title, subtitle = subtitle,
         x = sprintf("log2FC (%s / %s)", test, control), y = NULL) +
    enrich_theme() +
    theme(legend.position = "top")
}
