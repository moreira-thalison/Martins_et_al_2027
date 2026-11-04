#!/usr/bin/env Rscript
# ============================================================================
# 06_ggkegg_pathways.R  -  KEGG pathway maps with limma log2FC overlaid, using
# ggkegg (modern Pathview alternative: parses KGML -> ggplot/tidygraph, vector
# output). Colours each enzyme/gene box by the mean limma log2FC of the mapped
# protein(s). Red = higher in test (Lactate), blue = lower in test.
# Only boxes with |log2FC| > THRESH are coloured (default 1); the rest stay as
# the plain reference map, so the figure highlights the substantial changes.
#
# ID mapping: UniProt `protein_id` -> KEGG `afm:AFUA_*` gene id via
# KEGGREST::keggConv("afm","uniprot"). Multiple UniProt -> same afm gene are
# averaged; multiple genes per KGML box are combined by ggkegg (mean).
#
# Pathways = ALL the shortlisted cell-wall/glycan and mitochondria/energy maps.
# Edit PATHWAYS to add/remove. Requires internet (KGML + reference image each).
# Output: figures/KEGG Pathways/prot__pathway__<contrast>_<pid>__ggkegg.{png,pdf}
# (maps with no detected protein are skipped; coverage is reported in the log).
#
# Run:  Rscript scripts/06_ggkegg_pathways.R
# ============================================================================
suppressPackageStartupMessages({
  library(ggkegg); library(tidygraph); library(dplyr); library(ggraph)
  library(data.table); library(KEGGREST); library(ggplot2); library(scales)
})
options(timeout = 600)
args <- commandArgs(FALSE); self <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
SDIR <- dirname(self); source(file.path(SDIR, "config.R"))
PADJ  <- CFG$padj
THRESH <- 0            # colour only boxes with |log2FC| > THRESH
# only these KGML node shapes are drawable boxes; genes rendered as glycan-
# structure "line" nodes (e.g. in N-glycan maps) have NA coords and cannot be
# coloured, so they must NOT be counted as covered.
DRAWABLE <- c("rectangle", "roundrectangle")
# KEGG default node fills to make transparent so our log2FC rects show through
# (white reference boxes + green "gene present in afm" boxes + grey map links).
TRANSP <- c("#FFFFFF", "#BFFFBF", "#90EE90", "#BFBFFF", "#B3B3B3", "#CCCCCC", "#999999")

# ---- pathways to draw (all shortlisted cell-wall + energetic maps) ----------
PATHWAYS <- rbindlist(list(
  # cell wall / glycan / envelope
  data.table(pid="afm00520", cat="cell_wall", title="Amino sugar & nucleotide sugar metabolism (chitin)"),
  data.table(pid="afm00051", cat="cell_wall", title="Fructose and mannose metabolism (mannan)"),
  data.table(pid="afm00500", cat="cell_wall", title="Starch and sucrose metabolism (glucan/trehalose)"),
  data.table(pid="afm00040", cat="cell_wall", title="Pentose and glucuronate interconversions"),
  data.table(pid="afm00100", cat="cell_wall", title="Steroid biosynthesis (ergosterol)"),
  data.table(pid="afm00510", cat="cell_wall", title="N-Glycan biosynthesis"),
  data.table(pid="afm00513", cat="cell_wall", title="Various types of N-glycan biosynthesis"),
  data.table(pid="afm00515", cat="cell_wall", title="Mannose type O-glycan biosynthesis"),
  data.table(pid="afm00514", cat="cell_wall", title="Other types of O-glycan biosynthesis"),
  data.table(pid="afm00563", cat="cell_wall", title="GPI-anchor biosynthesis"),
  # mitochondria / energy / carbon
  data.table(pid="afm00190", cat="energetic", title="Oxidative phosphorylation"),
  data.table(pid="afm00020", cat="energetic", title="Citrate cycle (TCA cycle)"),
  data.table(pid="afm00620", cat="energetic", title="Pyruvate metabolism"),
  data.table(pid="afm00630", cat="energetic", title="Glyoxylate and dicarboxylate metabolism"),
  data.table(pid="afm00010", cat="energetic", title="Glycolysis / Gluconeogenesis"),
  data.table(pid="afm00071", cat="energetic", title="Fatty acid degradation"),
  data.table(pid="afm00061", cat="energetic", title="Fatty acid biosynthesis"),
  data.table(pid="afm01212", cat="energetic", title="Fatty acid metabolism (overview)"),
  data.table(pid="afm00280", cat="energetic", title="Valine, leucine and isoleucine degradation"),
  data.table(pid="afm00640", cat="energetic", title="Propanoate metabolism"),
  data.table(pid="afm00650", cat="energetic", title="Butanoate metabolism"),
  data.table(pid="afm01200", cat="energetic", title="Carbon metabolism (overview)"),
  # cellular processes
  data.table(pid="afm04138", cat="cellular", title="Autophagy - yeast"),
  data.table(pid="afm04139", cat="cellular", title="Mitophagy - yeast"),
  data.table(pid="afm04111", cat="cellular", title="Cell cycle - yeast"),
  data.table(pid="afm04149", cat="cellular", title="Peroxisome - yeast")
))

# ---- log2FC keyed by afm gene id -------------------------------------------
ct <- CFG$contrasts[[1]]; cn <- ct$name; tg <- ct$test; cg <- ct$control
d  <- fread(file.path(CFG$dir_results, sprintf("prot__de__%s__results_limma.tsv", cn)))
cv <- keggConv("afm", "uniprot")                        # up:UNIPROT -> afm:GENE
d[, afm := cv[paste0("up:", protein_id)]]
message(sprintf("mapped %d/%d proteins to afm gene ids", sum(!is.na(d$afm)), nrow(d)))
gl <- d[!is.na(afm), .(log2FC = mean(log2FC, na.rm = TRUE),
                       padj   = min(padj,   na.rm = TRUE)), by = afm]
vec <- setNames(gl$log2FC, gl$afm)

# symmetric colour limit shared across maps (robust: 99th pct, capped 2..4)
L <- max(2, min(4, ceiling(quantile(abs(vec), 0.99, na.rm = TRUE))))

# ---- draw each pathway ------------------------------------------------------
cov <- data.table()
for (i in seq_len(nrow(PATHWAYS))) {
  pid <- PATHWAYS$pid[i]; ttl <- PATHWAYS$title[i]; cate <- PATHWAYS$cat[i]
  message("\n=== ", pid, " [", cate, "] : ", ttl, " ===")
  g <- tryCatch(
    pathway(pid) |>
      mutate(lfc = node_numeric(vec)) |>
      mutate(drawable = graphics_type %in% DRAWABLE) |>
      mutate(lfc_thr = ifelse(drawable & !is.na(lfc) & abs(lfc) > THRESH, lfc, NA_real_)),
    error = function(e) { message("  [error] ", conditionMessage(e)); NULL })
  if (is.null(g)) { cov <- rbind(cov, data.table(pathway=pid, category=cate, title=ttl,
                        gene_nodes=NA, with_data=NA, coloured=NA, max_abs_lfc=NA_real_,
                        status="error")); next }
  nd <- as_tibble(as_tbl_graph(g), active = "nodes")
  drw    <- nd$drawable & !is.na(nd$lfc)                  # drawable boxes with data
  n_gene <- sum(nd$type == "gene" & nd$drawable, na.rm = TRUE)
  n_dat  <- sum(drw)
  n_col  <- sum(!is.na(nd$lfc_thr))
  mx     <- if (n_dat > 0) max(abs(nd$lfc[drw]), na.rm = TRUE) else NA_real_
  message(sprintf("  gene boxes: %d total, %d with data, %d coloured (|log2FC|>%g); max|log2FC|=%.2f",
                  n_gene, n_dat, n_col, THRESH, ifelse(is.na(mx), 0, mx)))

  # nothing detected in this map (and some, e.g. afm00514, lack a KEGG image) -> skip
  if (n_dat == 0) {
    cov <- rbind(cov, data.table(pathway=pid, category=cate, title=ttl, gene_nodes=n_gene,
                  with_data=0L, coloured=0L, max_abs_lfc=NA_real_, status="no_data_skipped"))
    message("  [skip] no detected proteins in this map"); next
  }
  cov <- rbind(cov, data.table(pathway=pid, category=cate, title=ttl, gene_nodes=n_gene,
                with_data=n_dat, coloured=n_col, max_abs_lfc=round(mx,2),
                status=ifelse(n_col>0,"ok","no_box_over_threshold")))
  xr <- range(nd$x, na.rm = TRUE); yr <- range(nd$y, na.rm = TRUE)

  p <- ggraph(g, layout = "manual", x = x, y = y) +
    geom_node_rect(aes(fill = lfc_thr, filter = !is.na(lfc_thr))) +
    overlay_raw_map(pid, transparent_colors = TRANSP) +
    scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                         midpoint = 0, limits = c(-0.588, 0.588), oob = scales::squish,
                         name = "log2FC") +
    labs(title = sprintf("%s  -  %s", pid, ttl),
         subtitle = sprintf("%s: red = %s, blue = %s | |log2FC|>%g: %d/%d boxes coloured | max|log2FC|=%.2f",
                            cn, dir_labels(tg, lower = TRUE)[["Up"]], dir_labels(tg, lower = TRUE)[["Down"]],
                            THRESH, n_col, n_dat, mx)) +
    theme_void() +
    theme(plot.background  = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 10, colour = "grey30"),
          legend.position = c(0.98, 0.5))

  stem <- sprintf("prot__pathway__%s_%s__ggkegg", cn, pid)
  base <- file.path(CFG$dir_figures_kegg, stem)
  base_pdf <- file.path(CFG$dir_figures_kegg_pdf, stem)
  w <- diff(xr) / 100; h <- diff(yr) / 100
  ok <- tryCatch({
    ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = 600, bg = "white", limitsize = FALSE)
    ggsave(paste0(base_pdf, ".pdf"), p, width = w, height = h, bg = "white", limitsize = FALSE)
    TRUE
  }, error = function(e) { message("  [error] save ", pid, ": ", conditionMessage(e)); FALSE })
  if (!ok) {
    for (f in c(paste0(base, ".png"), paste0(base_pdf, ".pdf"))) if (file.exists(f)) file.remove(f)
    cov[pathway == pid, status := "render_failed"]
  }
}

setorder(cov, category, -coloured, na.last = TRUE)
message("\nPathway coverage:"); print(cov)
writeLines(c("step\tggkegg_pathways", sprintf("time\t%s", Sys.time()),
  sprintf("contrast\t%s", cn),
  sprintf("mapping\tkeggConv afm<-uniprot\tcolour_limit\t+/-%g\tthreshold\t|log2FC|>%g", L, THRESH),
  sprintf("pathways\t%d (%d drawn ok)", nrow(PATHWAYS), sum(cov$status=="ok", na.rm=TRUE)),
  "colour\tmean limma log2FC per gene box (red=higher in test, blue=lower in test), only |log2FC|>THRESH"),
  file.path(CFG$dir_logs, "06_ggkegg_pathways.log"))
message("\nDONE: ggkegg pathway maps (", sum(cov$status=="ok", na.rm=TRUE), "/", nrow(PATHWAYS), " with coloured boxes)")
