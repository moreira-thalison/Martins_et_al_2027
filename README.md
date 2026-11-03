# *Aspergillus fumigatus* Af293 conidia DIA proteomics: Lactate vs Saline

Analysis code for Martins et al. (2027). Differential proteomics of
*A. fumigatus* Af293 conidia exposed to lactate versus saline, from
FragPipe / DIA-NN `combined_protein.tsv` (MaxLFQ), analysed with limma.

## Design

| Group   | N | Samples                                              | Acquisition |
|---------|---|------------------------------------------------------|-------------|
| Saline  | 3 | Saline_R1_rep2_1, Saline_R2_rep2_2, Saline_R3_rep2_3 | 2026-08-19  |
| Lactate | 3 | Lactate_R2_2, Lactate_R3_3, Lactate_R4_4             | 2026-08-07  |

The samples are chosen with `CFG$keep_samples` in `scripts/config.R` from
the full 14-sample matrix in `data/combined_protein.tsv`.

**Caveat:** each group comes from a different acquisition day, so group and
batch cannot be separated. The model is `~0 + group` without batch
correction.

## Methods (summary)

- 0 -> NA, log2; contaminants and reverse hits removed; a protein is kept if
  it has at least 2 valid values in at least 1 group.
- limma: quantile normalisation, Perseus-style imputation (downshift 1.8 SD,
  width 0.3, seed 123), `arrayWeights()`; significance = BH padj <= 0.05,
  no log2FC cutoff.
- Enrichment uses the detected proteome as the universe:
  KEGG (`afm`, ORA + GSEA), GO (curated UniProt, proteome UP000002530),
  COG (UniProt eggNOG xrefs -> eggNOG 5.0 level 2759).
- KEGG pathway maps (ggkegg), DEP heatmaps, detected-proteome Venn and
  composite panel figures.

## Pipeline

| Step | Script                            | Output                               |
|------|-----------------------------------|--------------------------------------|
| 00   | `00_make_samplesheet.R`           | `config/samples.tsv`                 |
| 01   | `01_DE_proteome_limma.R`          | DE tables, PCA, volcano, MA          |
| 03   | `03_qc_rle.R`                     | RLE QC figure                        |
| 02   | `02_kegg_enrichment.R`, `02b_...` | KEGG ORA + GSEA tables, bubble plots |
| 05   | `05_go_enrichment.R`, `05b_...`   | GO ORA + GSEA tables, bubble plots   |
| 05c  | `05c_cog_enrichment.R`            | COG ORA, profile and strip plots     |
| 06   | `06_ggkegg_pathways.R`            | coloured KEGG pathway maps           |
| 07   | `07_heatmap_lactate_cellwall.R`   | cell-wall protein heatmap            |
| 08   | `08_global_heatmap_venn.R`        | DEP heatmaps, detection Venn         |
| 09b  | `09b_panel_variants.R`            | composite panels                     |

`config.R` holds every path and threshold; `io.R` and `plot_utils.R` are
shared helpers.

## Layout

```
scripts/      pipeline
config/       samples.tsv (sample sheet that was run), experiment_annotation.tsv
data/         combined_protein.tsv, experiment_annotation.tsv (raw files not tracked)
annotation/   cached UniProt GO / eggNOG reference tables
env/          dependencies.R installer + sessionInfo of a known-good run
results/ figures/ objects/ logs/   generated, git-ignored
cache/        KEGG KGML downloads, git-ignored
```


## Raw data

The Thermo `.raw` files are kept in the ProteomeXchange data base under the ID XXXXXXX.
<!-- TODO: add PRIDE / MassIVE accession once deposited. -->
