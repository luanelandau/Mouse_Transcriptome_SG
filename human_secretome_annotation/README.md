# Human salivary-gland expression annotation

This directory contains the R script used to add Human Protein Atlas (HPA)
annotations and human–mouse orthology to the human salivary-gland expression
matrix.

The workflow produces separate annotated expression tables for the parotid
(PAR), sublingual (SL), and submandibular (SM) glands.

## Main script

`annotate_human_expression_HPA_orthology.R`

The script:

1. reads the gene-level human TPM matrix;
2. adds `Secretome location` and `Secretome function` from the HPA Protein
   Atlas table;
3. selects salivary-gland immunohistochemistry records and adds `Level`,
   `Reliability`, and `Cell type`;
4. adds human–mouse orthology classifications;
5. combines multiple mouse orthologs for the same human gene using ` | `;
6. labels genes without a mouse ortholog as `human-specific`;
7. calculates mean TPM separately for PAR, SL, and SM; and
8. writes one annotated table per gland.

## R requirements

The analysis requires R and the following packages:

```r
install.packages(c("dplyr", "readr"))
```

## Input files generated in this study

The script should be run from the root of the `Mouse_Transcriptome_SG`
repository.

### Human expression matrix

```text
miscelaneous_sheets/gene_expression_matrix_human_TPMs.csv
```

This is the gene-level human TPM matrix. The `Geneid` column contains HGNC
gene symbols. Sample columns follow these naming conventions:

```text
adult_PAR_*
adult_SL_*
adult_SM_*
```

### Human–mouse orthology

```text
orthology_assignment/orthologs_with_classification.csv
```

The columns used by the script are:

```text
human_gene
mouse_gene
ortholog_type
```

## Human Protein Atlas inputs

The full `proteinatlas.tsv` file is not included in this repository because of
its size. Download and decompress it locally. If `normal_ihc_data.tsv` is also
absent from a local copy of the repository, download and decompress that file
as described below. Set `hpa_data_dir` near the beginning of the R script to
the directory containing both files.

### Protein Atlas table

Required filename:

```text
proteinatlas.tsv
```

The HPA Protein Atlas table is available from the
[Human Protein Atlas downloadable-data page](https://www.proteinatlas.org/about/download).
The compressed download is named `proteinatlas.tsv.zip`.

The workflow uses:

```text
Gene
Secretome location
Secretome function
```

All non-missing HPA Secretome annotations are retained, including entries
classified as `Intracellular and membrane`.

### Normal-tissue immunohistochemistry

Required filename:

```text
normal_ihc_data.tsv
```

The immunohistochemistry table used for the original analysis was downloaded
as `normal_ihc_data.tsv.zip` from the
[HPA v24.1 tissue-data page](https://v24.proteinatlas.org/humanproteome/tissue/data).

The workflow selects:

```r
Tissue == "Salivary gland"
```

and retains:

```text
Gene name
Level
Reliability
Cell type
```

If a gene has multiple salivary-gland IHC records, distinct annotations are
combined using `; `.

## Directory configuration

Edit this line in `annotate_human_expression_HPA_orthology.R`:

```r
hpa_data_dir <- "/path/to/directory/containing/HPA/files"
```

The selected directory should contain:

```text
proteinatlas.tsv
normal_ihc_data.tsv
```

## Running the analysis

From the root of the `Mouse_Transcriptome_SG` repository, run:

```bash
Rscript human_secretome_annotation/annotate_human_expression_HPA_orthology.R
```

## Outputs

The script writes the following tables to
`miscelaneous_sheets/human expression/`:

```text
PAR_human_mastersheet_TPMs_annotated_with_orthology_HPA_v25.csv
SL_human_mastersheet_TPMs_annotated_with_orthology_HPA_v25.csv
SM_human_mastersheet_TPMs_annotated_with_orthology_HPA_v25.csv
```

Each output contains:

```text
Geneid
mouse_gene
ortholog_type
Secretome location
Secretome function
Level
Reliability
Cell type
Mean_TPM
gland-specific TPM columns
```
