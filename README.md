# *Aspergillus fumigatus* Af293 conidia DIA proteomics: Lactate vs Saline

Analysis code for Martins et al. (2027). Differential proteomics of
*A. fumigatus* Af293 conidia exposed to lactate versus saline, from
FragPipe / DIA-NN `combined_protein.tsv` (MaxLFQ), analysed with limma.

## Design

| Group   | N | Samples                                              | 
|---------|---|------------------------------------------------------|
| Saline  | 3 | Saline_R1, Saline_R2, Saline_R3                      |
| Lactate | 3 | Lactate_R1, Lactate_R2, Lactate_R3                   | 

The Uniprot proteome of the strain Af293 was used as database for the 
spectral serach instead of the CEA10's one due to be more complete
and with more confident proteins. The similarity bewtween the strains makes the 
results fully applicable to CEA10 though.

## Methods (summary)

- 0 -> NA, log2; contaminants and reverse hits removed; a protein is kept if
  it has at least 2 valid values in at least 1 group.
- limma: quantile normalisation, Perseus-style imputation (downshift 1.8 SD,
  width 0.3, seed 123), `arrayWeights()`; significance = BH padj <= 0.05,
  no log2FC cutoff.
- Enrichment uses the detected proteome as the universe:
  KEGG (`afm`, ORA + GSEA), GO (curated UniProt, proteome Af293 UP000002530),
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

`config.R` holds every path and threshold; `io.R` (input / DE tables),
`plot_utils.R` (figures) and `enrich_utils.R` (KEGG / GO / COG enrichment and
bubble plots) are shared helpers. Every script opens with a header listing its
inputs, outputs and numbered steps, and each step is marked in the code by a
`# N. Step name ----` line (shown in the RStudio document outline).

## Layout

```
scripts/      pipeline
config/       samples.tsv, experiment_annotation.tsv (output from FragPipe)
data/         combined_protein.tsv (output from FragPipe), experiment_annotation.tsv (output from FragPipe)
annotation/   cached UniProt GO / eggNOG reference tables
env/          dependencies.R installer + sessionInfo
results/ figures/ objects/ logs/   generated, git-ignored
cache/        KEGG KGML downloads, git-ignored
```

## Reproduce

```bash
Rscript env/dependencies.R
bash scripts/run_all.sh #ties all the scripts together
```

Developed on R 4.5.2. Core steps (00, 01, 03) run offline. The enrichment
and pathway steps need internet access (KEGG; UniProt/eggNOG only if
`annotation/` is missing). To run from another location, set
`DIA_PROJECT_ROOT=/path/to/project`.

## Raw data

The raw files are available at PRIDE under the project with the 
ID XXXXXXXX.

## Publication
This proteome analysis is part of the a paper coordinated by Flávia Martins (1st author)
and Frederico Soriani (last author) from Universidade Federal de Minas Gerais; the whole proteomics
(from conidia to data analysis) were coordinated by Thalison Moreira (PhD Candidate) and Clayton Borges
from Universidade Federal de Goiás at the Molecular Biology Laboratory and Multiuser Center for Omics 
Analysis of Biomolecules (CMAS).
