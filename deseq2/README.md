# DESeq2 Analyses

This folder contains the R scripts and outputs used for differential 
gene-expression and principal component analyses of the mouse RNA-seq 
dataset.

## Files

* `deseq2_sex_strain_final.R`
  Performs separate DESeq2 analyses for each tissue using the design `~ 
strain + sex`. It generates differential-expression results for sex and 
strain, as well as combined master tables.

* `PCA.R`
  Performs PCA using variance-stabilized gene-expression counts. It 
generates PCA plots for all tissues and for the salivary glands only.

## Output folders

* `deseq2_results/`
  Contains DESeq2 result tables, filtered differentially expressed gene 
tables, and master sheets.

* `figures/`
  Contains PCA figures in PNG and SVG formats.

## Usage

Run the scripts from the main `Mouse_Transcriptome_SG` project directory:

```r
source("deseq2/deseq2_sex_strain_final.R")
source("deseq2/PCA.R")
```

