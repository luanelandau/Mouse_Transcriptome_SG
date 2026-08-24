# Figure 3: chromosome 7 and kallikrein expression

This folder contains the scripts and source plots for Figure 3, focused on the chromosome 7 kallikrein (`Klk`) locus.

## Panels

- `Fig3a_chr7_Klk_density_prot+RNA.R`: plots chromosome 7 positions and spatial density for significant transcriptome and proteome genes, highlighting genes shared by both datasets.
- `Fig3b_Bar_plots_Klk_expression.R`: compares the proportion of submandibular-gland expression from Klk genes in male and female samples and prints a Wilcoxon test.
- `Fig_3c_Plot_KLK_genes.R`: compares mean female and male expression for significant Klk genes on a log10 scale with standard-deviation bars.

Required R packages: `tidyverse`, `ggrepel`, and `patchwork`.

## Usage

Run the scripts from the repository root after generating the mouse DESeq2 tables and annotated expression master sheets. The scripts contain project-specific working-directory paths that should be reviewed first. PNG and SVG source plots are written under `Figure_3/`; final panel assembly is not performed by these scripts.
