# Figure 1 analysis and source plots

This folder contains the scripts and computational source plots for Figure 1 panels 1a–c and 1e. Panel 1d and final figure assembly are not included.

| Panel | Analysis | Script |
|---|---|---|
| 1a–b | Top secreted genes in human and mouse salivary glands | `Figure1a-b_human_and_mouse_bubbles.R` |
| 1c | Contribution of orthology categories to secreted-gene expression | `Fig1c_barplots_orthology_secreted_genes_human_mouse.R` |
| 1e | Mouse saliva proteome versus gland transcript abundance | `Figure1e_proteome_correlation.R` |

The analyses use annotated human and mouse expression master sheets under `miscelaneous_sheets/`. Panel 1e also uses `Saliva_mouse_summary.csv`; its tables and plots are written under `outputs_proteome_x_expression/` and also contribute to Figure S6.

Run the scripts from the repository root after reviewing their project-specific paths:

```bash
Rscript Figure_1/Figure1a-b_human_and_mouse_bubbles.R
Rscript Figure_1/Fig1c_barplots_orthology_secreted_genes_human_mouse.R
Rscript Figure_1/Figure1e_proteome_correlation.R
```

PNG and SVG files are quantitative source plots. Selected outputs were edited and assembled in Inkscape for the final publication figure.
