# Salmon reanalysis workflow

Run the scripts in the following order:

1. **Process NCBI FASTQ files with fastp and Salmon**
   ```bash
   bash fastp+salmon.sh
   ```

2. **Create the transcript-to-gene index**
   ```bash
   Rscript tx2gene.R
   ```

3. **Import Salmon outputs and summarize transcript abundance at the gene level**
   ```bash
   Rscript transcript_to_gene_abundance.R
   ```

4. **Create species- and tissue-specific master sheets**
   ```bash
   Rscript create_mastersheet_species_tissue.R
   ```

5. **Run the sex-based DESeq2 analysis**
   ```bash
   Rscript deseq2_sex.R
   ```

6. **Annotate the Salmon TPM master sheets and generate Table S4**
   ```bash
   Rscript Table_S4_annotate_salmon_tpm_and_summarize_secreted_genes.R
   ```

   Annotated tissue master sheets and the secreted-expression summary are written under `salmon_annotated_secretion_orthology/`.

## Optional normalization check

Check the distribution of Salmon TPM values:

```bash
Rscript TPM_distribution_by_species_comparison_salmon_hisat.R
```

The workflow uses project-specific paths and a `quant_paths.txt` file listing Salmon outputs. Review these inputs before execution. The DESeq2 result tables are written under `deseq2_sex/`, and associated volcano plots under `figures_sex/`.
