# Human submandibular-gland DESeq2 analysis

This folder contains the human submandibular-gland (SM) sex-differential-expression analysis used for Table S11 and the human component of the Figure 2b volcano plot.

## Files

- `deseq2_hum.R` fits a sex-only DESeq2 model to six human SM samples from `miscelaneous_sheets/gene_expression_matrix_human_RAW_COUNTS.csv`. Female is the reference level, so positive log2 fold changes indicate higher expression in males.
- `DESeq2_sex_SM.csv` contains the complete human SM DESeq2 results.
- `Volcano_plot_sex_SM_hum_mouse.R` compares the human and mouse SM sex-DESeq2 results and produces the Figure 2b source plot.

## Usage

Run `deseq2_hum.R` from the repository root after creating the human raw-count matrix. The volcano-plot script additionally requires the mouse SM result from `deseq2/deseq2_results/`. Both scripts contain project-specific paths that should be reviewed before execution.

Main R packages: `DESeq2`, `tidyverse`, `ggplot2`, `ggrepel`, and `patchwork`.
