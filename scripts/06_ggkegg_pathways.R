#!/usr/bin/env Rscript
# =============================================================================
# 06_ggkegg_pathways.R
# KEGG pathway maps with the limma log2FC overlaid, drawn with ggkegg (KGML ->
# ggplot / tidygraph, vector output). Each enzyme / gene box is filled with the
# mean log2FC of its mapped protein(s): red = higher in test (Lactate), blue =
# lower in test. Boxes with |log2FC| <= THRESH stay as the plain reference map.
#
# ID mapping: UniProt protein_id -> KEGG afm:AFUA_* via
# KEGGREST::keggConv("afm", "uniprot"). Several UniProt ids on one afm gene are
# averaged; several genes in one KGML box are combined by ggkegg (mean).
#
# Pathways: the shortlisted cell-wall / glycan, energy and cellular-process
# maps in PATHWAYS (edit to add or remove). Needs internet (KGML + reference
# image per map); run from inside cache/kegg_kgml (run_all.sh does this).
#
# Input  : results/prot__de__<contrast>__results_limma.tsv
# Outputs: figures/KEGG Pathways/prot__pathway__<contrast>_<pid>__ggkegg.png
#          figures/KEGG Pathways/KEGG Figures PDF/prot__pathway__<contrast>_<pid>__ggkegg.pdf
#          logs/06_ggkegg_pathways.log (maps without detected proteins are skipped;
#          coverage per map is printed)
# Run    : Rscript scripts/06_ggkegg_pathways.R
#
# Steps
#   0. Setup and drawing constants
#   1. Pathway list
#   2. log2FC per afm gene
#   3. Draw each pathway
#      3.1 Parse KGML and attach log2FC
#      3.2 Coverage
#      3.3 Plot
#      3.4 Save
#   4. Coverage report and run log
# =============================================================================


# 0. Setup and drawing constants -----------------------------------------------
suppressPackageStartupMessages({
  library(ggkegg); library(tidygraph); library(dplyr); library(ggraph)
  library(data.table); library(KEGGREST); library(ggplot2); library(scales)
})
options(timeout = 600)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
SDIR <- if (length(script_file)) dirname(normalizePath(script_file)) else "scripts"
source(file.path(SDIR, "config.R"))
source(file.path(SDIR, "io.R"))

THRESH     <- 0       # colour only boxes with |log2FC| > THRESH
FILL_LIMIT <- 0.588   # colour scale saturates at |log2FC| = FILL_LIMIT (shared by all maps)
# Only these KGML node shapes are drawable boxes. Genes drawn as glycan-structure
# "line" nodes (e.g. in N-glycan maps) have NA coordinates and cannot be
# coloured, so they are not counted as covered.
DRAWABLE <- c("rectangle", "roundrectangle")
# KEGG default node fills made transparent so the log2FC boxes show through
# (white reference boxes, green "gene present in afm" boxes, grey map links).
TRANSP <- c("#FFFFFF", "#BFFFBF", "#90EE90", "#BFBFFF", "#B3B3B3", "#CCCCCC", "#999999")


# 1. Pathway list --------------------------------------------------------------
PATHWAYS <- rbindlist(list(
  # cell wall / glycan / envelope
  data.table(pid = "afm00520", cat = "cell_wall", title = "Amino sugar & nucleotide sugar metabolism (chitin)"),
  data.table(pid = "afm00051", cat = "cell_wall", title = "Fructose and mannose metabolism (mannan)"),
  data.table(pid = "afm00500", cat = "cell_wall", title = "Starch and sucrose metabolism (glucan/trehalose)"),
  data.table(pid = "afm00040", cat = "cell_wall", title = "Pentose and glucuronate interconversions"),
  data.table(pid = "afm00100", cat = "cell_wall", title = "Steroid biosynthesis (ergosterol)"),
  data.table(pid = "afm00510", cat = "cell_wall", title = "N-Glycan biosynthesis"),
  data.table(pid = "afm00513", cat = "cell_wall", title = "Various types of N-glycan biosynthesis"),
  data.table(pid = "afm00515", cat = "cell_wall", title = "Mannose type O-glycan biosynthesis"),
  data.table(pid = "afm00514", cat = "cell_wall", title = "Other types of O-glycan biosynthesis"),
  data.table(pid = "afm00563", cat = "cell_wall", title = "GPI-anchor biosynthesis"),
  # mitochondria / energy / carbon
  data.table(pid = "afm00190", cat = "energetic", title = "Oxidative phosphorylation"),
  data.table(pid = "afm00020", cat = "energetic", title = "Citrate cycle (TCA cycle)"),
  data.table(pid = "afm00620", cat = "energetic", title = "Pyruvate metabolism"),
  data.table(pid = "afm00630", cat = "energetic", title = "Glyoxylate and dicarboxylate metabolism"),
  data.table(pid = "afm00010", cat = "energetic", title = "Glycolysis / Gluconeogenesis"),
  data.table(pid = "afm00071", cat = "energetic", title = "Fatty acid degradation"),
  data.table(pid = "afm00061", cat = "energetic", title = "Fatty acid biosynthesis"),
  data.table(pid = "afm01212", cat = "energetic", title = "Fatty acid metabolism (overview)"),
  data.table(pid = "afm00280", cat = "energetic", title = "Valine, leucine and isoleucine degradation"),
  data.table(pid = "afm00640", cat = "energetic", title = "Propanoate metabolism"),
  data.table(pid = "afm00650", cat = "energetic", title = "Butanoate metabolism"),
  data.table(pid = "afm01200", cat = "energetic", title = "Carbon metabolism (overview)"),
  # cellular processes
  data.table(pid = "afm04138", cat = "cellular", title = "Autophagy - yeast"),
  data.table(pid = "afm04139", cat = "cellular", title = "Mitophagy - yeast"),
  data.table(pid = "afm04111", cat = "cellular", title = "Cell cycle - yeast"),
  data.table(pid = "afm04149", cat = "cellular", title = "Peroxisome - yeast")
))


# 2. log2FC per afm gene -------------------------------------------------------
ct <- CFG$contrasts[[1]]
cn <- ct$name
d  <- read_de_results(CFG, cn)
uniprot_to_afm <- keggConv("afm", "uniprot")             # up:UNIPROT -> afm:GENE
d[, afm := uniprot_to_afm[paste0("up:", protein_id)]]
message(sprintf("mapped %d/%d proteins to afm gene ids", sum(!is.na(d$afm)), nrow(d)))
gene_lfc <- d[!is.na(afm), .(log2FC = mean(log2FC, na.rm = TRUE)), by = afm]
lfc_by_gene <- setNames(gene_lfc$log2FC, gene_lfc$afm)


# 3. Draw each pathway ---------------------------------------------------------
coverage_row <- function(pid, cate, ttl, status, gene_nodes = NA, with_data = NA,
                         coloured = NA, max_abs_lfc = NA_real_)
  data.table(pathway = pid, category = cate, title = ttl, gene_nodes = gene_nodes,
             with_data = with_data, coloured = coloured, max_abs_lfc = max_abs_lfc, status = status)

dl_low <- dir_labels(ct$test, lower = TRUE)
cov <- data.table()
for (i in seq_len(nrow(PATHWAYS))) {
  pid <- PATHWAYS$pid[i]; ttl <- PATHWAYS$title[i]; cate <- PATHWAYS$cat[i]
  message("\n=== ", pid, " [", cate, "] : ", ttl, " ===")

  ## 3.1 Parse KGML and attach log2FC ----
  g <- tryCatch(
    pathway(pid) |>
      mutate(lfc = node_numeric(lfc_by_gene)) |>
      mutate(drawable = graphics_type %in% DRAWABLE) |>
      mutate(lfc_thr = ifelse(drawable & !is.na(lfc) & abs(lfc) > THRESH, lfc, NA_real_)),
    error = function(e) { message("  [error] ", conditionMessage(e)); NULL })
  if (is.null(g)) {
    cov <- rbind(cov, coverage_row(pid, cate, ttl, "error"))
    next
  }

  ## 3.2 Coverage ----
  nd     <- as_tibble(as_tbl_graph(g), active = "nodes")
  drw    <- nd$drawable & !is.na(nd$lfc)                   # drawable boxes with data
  n_gene <- sum(nd$type == "gene" & nd$drawable, na.rm = TRUE)
  n_dat  <- sum(drw)
  n_col  <- sum(!is.na(nd$lfc_thr))
  mx     <- if (n_dat > 0) max(abs(nd$lfc[drw]), na.rm = TRUE) else NA_real_
  message(sprintf("  gene boxes: %d total, %d with data, %d coloured (|log2FC|>%g); max|log2FC|=%.2f",
                  n_gene, n_dat, n_col, THRESH, ifelse(is.na(mx), 0, mx)))
  # nothing detected in this map (some, e.g. afm00514, also lack a KEGG image)
  if (n_dat == 0) {
    cov <- rbind(cov, coverage_row(pid, cate, ttl, "no_data_skipped", n_gene, 0L, 0L))
    message("  [skip] no detected proteins in this map")
    next
  }
  cov <- rbind(cov, coverage_row(pid, cate, ttl, ifelse(n_col > 0, "ok", "no_box_over_threshold"),
                                 n_gene, n_dat, n_col, round(mx, 2)))

  ## 3.3 Plot ----
  p <- ggraph(g, layout = "manual", x = x, y = y) +
    geom_node_rect(aes(fill = lfc_thr, filter = !is.na(lfc_thr))) +
    overlay_raw_map(pid, transparent_colors = TRANSP) +
    scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0,
                         limits = c(-FILL_LIMIT, FILL_LIMIT), oob = scales::squish, name = "log2FC") +
    labs(title = sprintf("%s  -  %s", pid, ttl),
         subtitle = sprintf("%s: red = %s, blue = %s | |log2FC|>%g: %d/%d boxes coloured | max|log2FC|=%.2f",
                            cn, dl_low[["Up"]], dl_low[["Down"]], THRESH, n_col, n_dat, mx)) +
    theme_void() +
    theme(plot.background  = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 10, colour = "grey30"),
          legend.position = c(0.98, 0.5))

  ## 3.4 Save (size follows the KGML map; partial files removed on failure) ----
  stem <- sprintf("prot__pathway__%s_%s__ggkegg", cn, pid)
  png_file <- file.path(CFG$dir_figures_kegg, paste0(stem, ".png"))
  pdf_file <- file.path(CFG$dir_figures_kegg_pdf, paste0(stem, ".pdf"))
  w <- diff(range(nd$x, na.rm = TRUE)) / 100
  h <- diff(range(nd$y, na.rm = TRUE)) / 100
  saved <- tryCatch({
    ggsave(png_file, p, width = w, height = h, dpi = 600, bg = "white", limitsize = FALSE)
    ggsave(pdf_file, p, width = w, height = h, bg = "white", limitsize = FALSE)
    TRUE
  }, error = function(e) { message("  [error] save ", pid, ": ", conditionMessage(e)); FALSE })
  if (!saved) {
    for (f in c(png_file, pdf_file)) if (file.exists(f)) file.remove(f)
    cov[pathway == pid, status := "render_failed"]
  }
}


# 4. Coverage report and run log -----------------------------------------------
setorder(cov, category, -coloured, na.last = TRUE)
message("\nPathway coverage:")
print(cov)
n_ok <- sum(cov$status == "ok", na.rm = TRUE)
writeLines(c(
  "step\tggkegg_pathways",
  sprintf("time\t%s", Sys.time()),
  sprintf("contrast\t%s", cn),
  sprintf("mapping\tkeggConv afm<-uniprot\tcolour_limit\t+/-%g\tthreshold\t|log2FC|>%g", FILL_LIMIT, THRESH),
  sprintf("pathways\t%d (%d drawn ok)", nrow(PATHWAYS), n_ok),
  "colour\tmean limma log2FC per gene box (red=higher in test, blue=lower in test), only |log2FC|>THRESH"),
  file.path(CFG$dir_logs, "06_ggkegg_pathways.log"))
message("\nDONE: ggkegg pathway maps (", n_ok, "/", nrow(PATHWAYS), " with coloured boxes)")
