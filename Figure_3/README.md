# Figure 3

Scripts and output files for Figure 3. Run the scripts from R; each script sets the working directory to the project folder and writes its output to `Figure_3/`.

- `Fig3a_chr7_Klk_density_prot+RNA.R`: plots chromosome 7 positions and spatial density for significant transcriptome and proteome genes, highlighting genes shared by both datasets.
- `Fig3b_Bar_plots_Klk_expression.R`: compares the proportion of submandibular-gland expression from Klk genes in male and female samples and prints a Wilcoxon test.
- `Fig_3c_Plot_KLK_genes.R`: compares mean female and male expression for significant Klk genes on a log10 scale with standard-deviation bars.

Required R packages: `tidyverse`, `ggrepel`, and `patchwork`.
