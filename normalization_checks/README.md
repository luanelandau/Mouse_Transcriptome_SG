# Normalization and expression checks

This folder contains supporting checks of expression thresholds and an independent Salmon-based reanalysis of the mouse RNA-seq data.

## Contents

- `number_of_genes_expressed_per_gland/count_expressed_genes_by_tissue.R` counts genes with TPM greater than 2 in each sample and summarizes the counts by PAR, SM, SL, pancreas, and liver. Its CSV outputs are stored in the same directory.
- `salmon_reanalysis/` contains FASTQ processing, Salmon quantification, transcript-to-gene summarization, TPM comparisons, Table S4 generation, and the Salmon-based sex-DESeq2 results used for Table S13. See `salmon_reanalysis/README.md`.

The scripts use annotated expression master sheets and project-specific paths. Review those paths and the nested workflow instructions before running the analyses.
