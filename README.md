# Mouse salivary-gland transcriptome and secretome analyses

This repository contains the analysis scripts and supporting data used for:

> **Gene Expansion and Regulatory Rewiring Shape Sex-Biased Evolution of the Mouse Submandibular Gland Secretome**
>
> Luane Jandira Bueno Landau, Shikha Jain, Nathan Griffin, Achisha Saikia, Jill M. Kramer, Sarah Knox, Stefan Ruhl, and Omer Gokcumen
>
> [https://doi.org/10.64898/2026.03.18.712472](https://doi.org/10.64898/2026.03.18.712472)

The repository documents the workflows used to process mouse and human gene-expression data, annotate orthology and secretome status, compare expression across species, and generate source plots for the manuscript. Some final publication figures were assembled or visually edited after the quantitative source plots were generated; those cases are noted in the relevant directory README.

## Analysis overview

The main workflow is:

1. Generate gene-level counts from the raw sequencing reads using the cluster RNA-seq workflow.
2. Calculate TPM values and combine individual samples into mouse and human expression matrices (`TPM/`).
3. Assign human–mouse orthology classes (`orthology_assignment/`).
4. Add secretome, functional, and orthology annotations to the expression master sheets (`Mouse_secretome_annotaton/` and `human_secretome_annotation/`).
5. Use the annotated master sheets for differential-expression and cross-species correlation analyses (`deseq2/`, `deseq2_hum/`, and `expression_correlations/`).
6. Generate the quantitative source plots for Figures 1–4 (`Figure_1/` through `Figure_4/`).
7. Analyze genomic clustering and promoter motifs (`cluster_analysis/`).

Run scripts from the repository root unless a directory README says otherwise. Several scripts contain project-specific paths or parameters that should be reviewed before execution.

## Manuscript figures and supplementary materials

Figures 1–4 have dedicated directories (`Figure_1/` through `Figure_4/`); brief panel-level descriptions are provided below and details are in each directory README.

### Supplementary figures

| Manuscript item | Description and provenance | Main scripts or source files |
|---|---|---|
| **Figure S1** | Principal-component analysis of the mouse expression data. | `deseq2/PCA.R` |
| **Figure S2** | Ensembl-based verification of the mouse–human orthology analyses in Figure 1. | Scripts, annotated master sheets, and plots under `Figure_S2/` |
| **Figure S3** | Mouse–human expression correlations across three salivary glands, pancreas, and liver. | `expression_correlations/analyze_mouse_human_correlations_by_tissue.R`; outputs under `expression_correlations/` |
| **Figure S4** | Bootstrap comparisons of the tissue correlations shown in Figure S3. | `expression_correlations/compare_correlations_between_tissues_bootstrap.R`; `expression_correlations/figures/ALL_5_tissues_spearman_bootstrap_comparisons.png` |
| **Figure S5** | Gel image. | Not generated or stored in this repository. |
| **Figure S6** | Correlations between the mouse salivary proteome and gene expression in each salivary gland. | `Figure_1/Figure1e_proteome_correlation.R`; `Figure_1/outputs_proteome_x_expression/` |
| **Figure S7** | Sex-stratified expression of sialyltransferase genes. | `Figure_2/Fig2e_sialic_acids_plot.R`; `Figure_2/boxplot_sialic_acids.png` |
| **Figure S8** | Immunofluorescence figure. | Not generated or stored in this repository. |
| **Figure S9** | CAFE alignment figure. | Not generated or stored in this repository. |
| **Figure S10** | Methods schematic. | Created in BioRender; not stored in this repository. |

### Supplementary tables

| Manuscript item | Description and provenance | Main scripts or source files |
|---|---|---|
| **Table S1** | Raw gene-level expression counts for the CD1 and C57 mouse samples, produced from FASTQ files and combined into a master sheet. | `processing_fastq/`; `readCounts_CD1/`; `readCounts_C57/`; `TPM/generate_mastersheets_RAW_and_TPM.R` |
| **Table S2** | Mouse saliva proteomics data. | `Figure_1/Saliva_mouse_summary.csv` (input dataset; not generated here) |
| **Table S3** | Mouse gland expression master sheets annotated for secretion, function, and human–mouse orthology. | `orthology_assignment/`; `Mouse_secretome_annotaton/`; outputs under `miscelaneous_sheets/mouse_expression/` |
| **Table S4** | Salmon-based mouse TPM master sheets and summary of secreted-gene expression. | `normalization_checks/salmon_reanalysis/Table_S4_annotate_salmon_tpm_and_summarize_secreted_genes.R`; outputs under `normalization_checks/salmon_reanalysis/salmon_annotated_secretion_orthology/` |
| **Table S5** | Human salivary-gland, pancreas, and liver expression master sheets with secretion and orthology annotations. Salivary-gland counts were processed with the mouse FASTQ workflow; GTEx pancreas and liver TPMs were averaged across individuals and annotated separately. | `readCounts_human/`; `TPM/generate_mastersheets_counts_human.R`; `human_secretome_annotation/`; outputs under `miscelaneous_sheets/human expression/` |
| **Table S6** | Pairwise bootstrap comparisons of mouse–human tissue correlations. | `expression_correlations/ALL_5_tissues_spearman_bootstrap_comparisons.csv` |
| **Table S7** | Additional proteomics results. | Not generated or stored in this repository. |
| **Table S8** | Sex differences in mouse saliva protein abundance. | `Figure_2/Saliva_mouse_sex_diff.csv` (input dataset; not generated here) |
| **Table S9** | Mouse DESeq2 results for sex and strain in PAR, SM, SL, pancreas, and liver. | `deseq2/deseq2_sex_strain_final.R`; outputs under `deseq2/deseq2_results/` |
| **Table S10** | Chi-square comparison of the number/proportion of sex-associated DEGs among tissues. | `deseq2/Table_S10_chi_square_DEGs.R`; `deseq2/Table_S10_chi_square_DEGs.csv` |
| **Table S11** | Human submandibular-gland DESeq2 sex-differential-expression results. | `deseq2_hum/deseq2_hum.R`; `deseq2_hum/DESeq2_sex_SM.csv` |
| **Table S12** | Cluster analysis results. | Scripts and results under `deseq2_hum/` |
| **Table S13** | Mouse sex-differential-expression results from the Salmon normalization check. | `normalization_checks/salmon_reanalysis/deseq2_sex.R`; outputs under `normalization_checks/salmon_reanalysis/deseq2_sex/` |
| **Table S14** | HOMER promoter-motif results for the `Klk` cluster. | `cluster_analysis/cluster_klk/` |
| **Table S15** | BLAST results for `Klk` genes in other species. | BLAST command is described in the paper; scripts and results are not stored here. |

### Figure 1

The scripts for Figure 1 panels **1a, 1b, 1c, and 1e** are in `Figure_1/`:

| Panel | Analysis | Script |
|---|---|---|
| **1a–b** | Human and mouse salivary-gland expression bubble plots | `Figure_1/Figure1a-b_human_and_mouse_bubbles.R` |
| **1c** | Contribution of orthology categories to secreted-gene expression | `Figure_1/Fig1c_barplots_orthology_secreted_genes_human_mouse.R` |
| **1e** | Mouse saliva proteome versus gland transcript abundance | `Figure_1/Figure1e_proteome_correlation.R` |

These analyses use the annotated human and mouse master sheets described above. See [`Figure_1/README.md`](Figure_1/README.md) for inputs, outputs, dependencies, and notes about final figure assembly.

### Figure 2

The `Figure_2/` directory contains the scripts, inputs, and source plots for panels **2a–c and 2e**:

| Panel | Analysis | Script |
|---|---|---|
| **2a** | Overlap of mouse sex-biased genes across PAR, SM, SL, pancreas, and liver | `Figure_2/Fig2a_upset_plot.R` |
| **2b** | Mouse and human SM sex-differential-expression volcano plots | `Figure_2/Fig2b_volcano_plot_sex_SM_hum_mouse.R` |
| **2c** | Integration of sex-dimorphic saliva proteins with gland transcriptomes | `Figure_2/Fig2c_Sex_analysis_proteome_2.R` |
| **2e** | Sex-stratified expression of sialyltransferase genes in mouse SM | `Figure_2/Fig2e_sialic_acids_plot.R` |

The folder also contains the external saliva differential-proteomics input, PNG/SVG plot outputs, and retained earlier plot versions. Experimental images and final panel assembly are not generated by these scripts. See [`Figure_2/README.md`](Figure_2/README.md) for thresholds, inputs, outputs, dependencies, and figure-assembly notes.

### Figure 3

The `Figure_3/` directory contains three analyses focused on the chromosome 7 kallikrein (`Klk`) locus:

| Panel | Analysis | Script |
|---|---|---|
| **3a** | Chromosome 7 positions and spatial density of significant transcriptomic and proteomic genes | `Figure_3/Fig3a_chr7_Klk_density_prot+RNA.R` |
| **3b** | Proportion of SM expression contributed by `Klk` genes in male and female samples | `Figure_3/Fig3b_Bar_plots_Klk_expression.R` |
| **3c** | Mean female-versus-male expression of significant `Klk` genes | `Figure_3/Fig_3c_Plot_KLK_genes.R` |

The directory includes the resulting PNG and SVG source plots. See [`Figure_3/README.md`](Figure_3/README.md) for a concise script and dependency summary.

### Figure 4: TAD analyses

The TAD component of **Figure 4** is under `Figure_4/`. The main workflow uses `Figure_4/hic_25kb.sh` with `Figure_4/CH12_25kb_select.ini`. The nested `Figure_4/Fig_4c_TADs_sex_differences/` workflow compares male and female TADs using `hic_25kb_new.sh` and `ser_gra_svg_25kb_new.ini`.

The accompanying chromatin tracks were obtained as processed bigWig files from ENCODE and include ATAC-seq, CTCF, RAD21, and H3K4me3 datasets. See [`Figure_4/README.md`](Figure_4/README.md) and [`Figure_4/Fig_4c_TADs_sex_differences/README.md`](Figure_4/Fig_4c_TADs_sex_differences/README.md) for accessions and analysis notes.

## Repository guide

| Directory | Contents |
|---|---|
| [`processing_fastq/`](processing_fastq/) | FASTQ trimming, alignment, and gene-counting shell scripts; see [`processing_fastq/README.md`](processing_fastq/README.md) |
| [`processing_RNA/`](processing_RNA/) | Reserved directory for upstream RNA-processing materials; currently empty |
| [`readCounts_C57/`](readCounts_C57/) and [`readCounts_CD1/`](readCounts_CD1/) | Mouse featureCounts outputs used to build the expression matrices |
| [`readCounts_human/`](readCounts_human/) | Human count inputs |
| [`TPM/`](TPM/) | TPM calculation and expression-matrix generation; see [`TPM/README.md`](TPM/README.md) |
| [`orthology_assignment/`](orthology_assignment/) | Human–mouse orthology classification; see [`orthology_assignment/README.md`](orthology_assignment/README.md) |
| [`Mouse_secretome_annotaton/`](Mouse_secretome_annotaton/) | Mouse secretome and orthology annotation; see [`Mouse_secretome_annotaton/README.md`](Mouse_secretome_annotaton/README.md) |
| [`human_secretome_annotation/`](human_secretome_annotation/) | Human secretome, HPA, and orthology annotation; see [`human_secretome_annotation/README.md`](human_secretome_annotation/README.md) |
| [`miscelaneous_sheets/`](miscelaneous_sheets/) | Intermediate and final expression master sheets used by downstream analyses |
| [`deseq2/`](deseq2/) | Differential-expression and PCA analyses |
| [`deseq2_hum/`](deseq2_hum/) | Human SM sex differential-expression analysis and mouse–human volcano-plot source script |
| [`expression_correlations/`](expression_correlations/) | Cross-species expression-correlation analyses; see [`expression_correlations/README.md`](expression_correlations/README.md) |
| [`Figure_1/`](Figure_1/) | Figure 1 source-plot scripts, proteomics input, and outputs |
| [`Figure_2/`](Figure_2/) | Figure 2 sex-bias, proteome/transcriptome integration, and sialyltransferase source plots; see [`Figure_2/README.md`](Figure_2/README.md) |
| [`Figure_3/`](Figure_3/) | Figure 3 chromosome 7 and `Klk`-expression source plots; see [`Figure_3/README.md`](Figure_3/README.md) |
| [`Figure_4/`](Figure_4/) | Figure 4 TAD workflows, including the nested male–female comparison; see [`Figure_4/README.md`](Figure_4/README.md) |
| [`Figure_S2/`](Figure_S2/) | Figure S2 Ensembl-based orthology verification; see [`Figure_S2/README.md`](Figure_S2/README.md) |
| [`cluster_analysis/`](cluster_analysis/) | SM sex-DEG genomic-clustering permutations and HOMER motif analysis of the `Klk` cluster; see the [`permutation`](cluster_analysis/permutations/README.md) and [`Klk motif`](cluster_analysis/cluster_klk/README.md) READMEs |

The directory names `Mouse_secretome_annotaton` and `miscelaneous_sheets` retain their original spellings so that existing script paths continue to work.

## Data availability and large files

A path on a computing cluster cannot be linked directly from GitHub for public access: it is only meaningful to users who can access that cluster and filesystem. For a file needed to reproduce the analysis, use one of these approaches:

- Commit it directly when it is modest in size and redistribution is permitted.
- Use Git LFS for a versioned file that is too large for ordinary Git history.
- Deposit raw sequencing reads and large processed datasets in an appropriate public repository (for example GEO/SRA, Zenodo, or Figshare), then add the accession or DOI here.
- Use a GitHub Release for a stable downloadable snapshot that does not need to live in the normal source history.

For long-term reproducibility, a public data-repository accession is preferable to a cluster path. Raw FASTQ files generally belong in a domain repository such as SRA; GitHub should contain the scripts, small inputs, metadata, and links needed to reconstruct the analysis.

## Reproducing the analyses

Detailed package requirements and commands are documented in the directory-level READMEs. A typical execution order is:

```bash
Rscript TPM/Calculate_TPM.R <featureCounts-file>
Rscript TPM/generate_mastersheets_RAW_and_TPM.R
Rscript orthology_assignment/orthology_assignment.R
Rscript Mouse_secretome_annotaton/01_secretome_annotation_dictionary.R
Rscript Mouse_secretome_annotaton/02-annotate_each_gland_secreted_non_secreted_mouse.R
Rscript Mouse_secretome_annotaton/03-adding_orthology.R
Rscript human_secretome_annotation/annotate_human_expression_HPA_orthology.R
```

The figure and correlation scripts can then be run after their documented input tables have been created. This repository does not yet provide a single automated end-to-end workflow, and external input datasets must be obtained as described in the relevant README files.

## Citation

If you use this repository, please cite the paper above. A versioned archival DOI for the code repository should also be added here when one is available.
