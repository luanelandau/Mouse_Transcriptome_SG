# Mouse–Human Expression Correlations

This folder contains scripts and output files used to compare the expression of secreted one-to-one orthologs between mouse and human across five tissues:

- PAR: parotid gland
- SM: submandibular gland
- SL: sublingual gland
- PANC: pancreas
- LIV: liver

## Scripts

### `correlations_secreted_genes_with_CI_TPM>2_stack.R`

This script builds the unified mouse–human expression tables for each tissue and calculates the mouse–human Spearman correlation.

For each tissue, it:

1. Reads the human and mouse annotated expression mastersheets.
2. Retains secreted one-to-one orthologs.
3. Joins the human and mouse tables using the ortholog assignments.
4. Filters genes according to `MIN_TPM` and `FILTER_MODE`.
5. Calculates the Spearman correlation and a bootstrap 95% confidence interval.
6. Creates a combined correlation figure with all five tissues using shared axis limits.

Important parameters near the beginning of the script are:

```r
MIN_TPM <- 0
FILTER_MODE <- "both"
N_BOOTSTRAPS <- 5000
```

- `MIN_TPM` sets the minimum mean TPM.
- `FILTER_MODE = "both"` requires the gene to meet the TPM threshold in both species.
- `FILTER_MODE = "either"` requires the gene to meet the threshold in at least one species.

Although the script filename contains `TPM>2`, the current script uses `MIN_TPM <- 0`. Change this value to `2` when a TPM threshold of 2 is required.

Main outputs include:

```text
<TISSUE>_secreted_one_to_one_TPM<MIN_TPM>_<FILTER_MODE>_unified.csv
figures/ALL_secreted_one_to_one_TPM<MIN_TPM>_<FILTER_MODE>_correlations_two_per_row_shared_axes.png
```

The unified files contain:

- `human_gene`
- `mouse_gene`
- `TPM_human`
- `TPM_mouse`

### `correlation_glands_confidence_interval.R`

This script compares the mouse–human Spearman correlations among the five tissues.

For every pair of tissues, it:

1. Keeps only ortholog pairs present in both tissue tables.
2. Calculates the Spearman correlation for each tissue using the same genes.
3. Calculates the difference between the correlations.
4. Resamples the shared ortholog pairs 10,000 times.
5. Calculates a bootstrap 95% confidence interval for the correlation difference.
6. Calculates bootstrap p-values and applies Benjamini–Hochberg correction across the 10 tissue comparisons.
7. Creates a plot showing the observed differences and confidence intervals.

The difference is always calculated as:

```text
correlation of tissue 1 - correlation of tissue 2
```

A positive value means tissue 1 has the higher correlation. A negative value means tissue 2 has the higher correlation.

Main outputs are:

```text
ALL_5_tissues_spearman_bootstrap_comparisons.csv
ALL_5_tissues_spearman_bootstrap_values.csv
figures/ALL_5_tissues_spearman_bootstrap_comparisons.png
```

The comparisons file contains one row per tissue comparison. The bootstrap-values file contains the results from every bootstrap iteration and is much larger.

## Recommended order

Run the scripts from the project directory:

```bash
Rscript "expression_correlations/correlations_secreted_genes_with_CI_TPM>2_stack.R"
Rscript "expression_correlations/correlation_glands_confidence_interval.R"
```

Run the correlation script first because it creates the unified tissue tables used by the tissue-comparison script.

## Current folder contents

The files named like:

```text
PAR_secreted_one_to_one_TPM0_both_unified.csv
SM_secreted_one_to_one_TPM0_both_unified.csv
SL_secreted_one_to_one_TPM0_both_unified.csv
PANC_secreted_one_to_one_TPM0_both_unified.csv
LIV_secreted_one_to_one_TPM0_both_unified.csv
```

were generated using a minimum TPM of 0 in both species.
