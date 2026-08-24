# Figure S2: Ensembl orthology verification

## Overview

This directory contains an independent verification of the mouse–human orthology analysis used in the paper. The original analysis used orthology assignments from The Jackson Laboratory; the verification presented here instead uses the Ensembl mouse–human orthology table downloaded from BioMart.

The workflow generates Figure S2, which reproduces the analyses and layout of Figure 1 using Ensembl orthology annotations. Comparing the two results provides a check that the conclusions are robust to the choice of orthology source.

## Input orthology data

`mouse_human_orthologs_ENS_v115.csv` is the mouse–human orthology table downloaded from Ensembl BioMart (Ensembl release 115).

`orthologs_with_classification_ensembl.csv` contains the Ensembl ortholog records formatted and classified for use by the downstream annotation scripts.

The scripts also use the original human and mouse expression mastersheets stored elsewhere in the parent project.

## Workflow and scripts

Run the scripts in the following order:

1. **Summarize and compare orthology assignments**

   `summary_orthology_comparison.R` compares the Ensembl orthology classifications with the Jackson Laboratory classifications used in the original analysis. It produces summary tables and a visualization of the agreement and differences between the two orthology sources.

2. **Add Ensembl orthology annotations to the expression mastersheets**

   - `adding_orthology_to_mouse_sheets_update_with_ensembl.R` annotates the mouse mastersheets with human orthologs and Ensembl orthology classifications.
   - `adding_orthology_to_human_sheets_update_with_ensembl.R` annotates the human mastersheets with mouse orthologs and Ensembl orthology classifications.

   These scripts follow the same annotation procedure used previously for the Jackson Laboratory orthology data, but substitute the Ensembl-derived orthology table. Genes without an identified ortholog are classified as mouse-specific or human-specific, as appropriate.

3. **Generate the Figure S2 bubble plots**

   `Figure1_human_and_mouse_bubbles_ensembl.R` reproduces the Figure 1 bubble-plot analysis using the Ensembl-annotated mastersheets. It shows the top secreted genes in the three mouse and human salivary glands, colored by orthology class.

4. **Generate the Figure S2 bar plots**

   `Fig1C_barplots_orthology_secreted_genes_human_mouse_ENSEMBLorthology.R` reproduces the Figure 1C bar plots using Ensembl orthology. The plots summarize the percentage of total expression contributed by one-to-one, lineage-specific, and other orthology categories.

## Outputs

### `mastersheets_ensembl/`

Contains the human and mouse expression mastersheets annotated with Ensembl orthology assignments. These are the Ensembl-based counterparts of the mastersheets used for the original Figure 1 analysis.

### `figures/`

Contains the PNG and SVG files used to assemble Figure S2:

- `bubbles_mouse_top_human_bottom_3glands_ensembl.*` — mouse and human salivary-gland bubble plots.
- `percent_total_expression_orthology_secreted_combined_horizontal_ensembl.*` — bar plots showing the contribution of each orthology category to total expression.

### Orthology-comparison outputs

- `orthology_summary_ensembl_vs_jackson.csv` — summarized comparison of Ensembl and Jackson Laboratory orthology classifications.
- `orthology_type_counts_ensembl_vs_jackson.csv` — counts by orthology category and source.
- `orthology_type_counts_ensembl_vs_jackson.png` — visualization of the orthology-category counts.

## R dependencies

The scripts use the following R packages:

- `dplyr`
- `forcats`
- `ggforce`
- `ggplot2`
- `packcircles`
- `patchwork`
- `purrr`
- `readr`
- `stringr`
- `tidyr`

## Reproducibility note

Several scripts contain absolute paths to the project directory. If the project is moved to another location, update `project_dir`, `base_dir`, `in_dir`, `human_dir`, or `setwd()` in the relevant scripts before running the workflow.
