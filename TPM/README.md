## Calculate TPM and Generate Gene-Expression Mastersheets

This workflow uses two R scripts:

1. `Calculate_TPM.R` calculates TPM values for each individual featureCounts file.
2. `generate_mastersheets_RAW_and_TPM.R` combines all individual C57 and CD1 files into raw-count and TPM mastersheets.

### Requirements

Install the required R packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("edgeR")
install.packages(c("readr", "dplyr", "purrr", "stringr"))
```

### Step 1: Calculate TPM for Each Sample

The input files should be featureCounts output files ending in:

```text
.counts.txt
```

Run the script for all C57 and CD1 samples:

```bash
for directory in readCounts_C57 readCounts_CD1; do
    for i in "$directory"/*.counts.txt; do
        Rscript Calculate_TPM.R "$i"
    done
done
```

For each input file, the script creates a new file ending in:

```text
.counts_with_tpm.txt
```

For example:

```text
Mous-1B-PAR-Mal-L.counts.txt
Mous-1B-PAR-Mal-L.counts_with_tpm.txt
```

### Step 2: Generate the Mastersheets

After all TPM files have been created, run:

```bash
Rscript generate_mastersheets_RAW_and_TPM.R
```

Before running the script, update its working directory so that it points to the main project folder.

The script reads the `.counts.txt` and `.counts_with_tpm.txt` files from:

```text
readCounts_C57/
readCounts_CD1/
```

It creates two combined matrices:

```text
miscelaneous_sheets/gene_expression_matrix_C57_CD1_RAW_COUNTS.csv
miscelaneous_sheets/gene_expression_matrix_C57_CD1_TPMs.csv
```

The raw-count matrix contains the original featureCounts values, while the TPM matrix contains the calculated TPM values. Each sample is stored in a separate column.
