# Count the number of genes with TPM > 2 in every sample.
# Columns 1 through 9 contain annotations and summary values.
# Therefore, the individual sample TPM columns begin at column 10.

# Read the five tissue files --------------------------------------------------

PAR <- read.csv(
  "../../miscelaneous_sheets/mouse_expression/PAR_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv",
  check.names = FALSE
)

SM <- read.csv(
  "../../miscelaneous_sheets/mouse_expression/SM_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv",
  check.names = FALSE
)

SL <- read.csv(
  "../../miscelaneous_sheets/mouse_expression/SL_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv",
  check.names = FALSE
)

PANC <- read.csv(
  "../../miscelaneous_sheets/mouse_expression/PANC_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv",
  check.names = FALSE
)

LIV <- read.csv(
  "../../miscelaneous_sheets/mouse_expression/LIV_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv",
  check.names = FALSE
)


# Count genes with TPM > 2 in each sample ------------------------------------

PAR_counts <- colSums(PAR[, 10:ncol(PAR)] > 2, na.rm = TRUE)
SM_counts <- colSums(SM[, 10:ncol(SM)] > 2, na.rm = TRUE)
SL_counts <- colSums(SL[, 10:ncol(SL)] > 2, na.rm = TRUE)
PANC_counts <- colSums(PANC[, 10:ncol(PANC)] > 2, na.rm = TRUE)
LIV_counts <- colSums(LIV[, 10:ncol(LIV)] > 2, na.rm = TRUE)


# Put the counts into one table ----------------------------------------------

PAR_results <- data.frame(
  Tissue = "PAR",
  Sample = names(PAR_counts),
  Expressed_genes_TPM_greater_than_2 = as.integer(PAR_counts)
)

SM_results <- data.frame(
  Tissue = "SM",
  Sample = names(SM_counts),
  Expressed_genes_TPM_greater_than_2 = as.integer(SM_counts)
)

SL_results <- data.frame(
  Tissue = "SL",
  Sample = names(SL_counts),
  Expressed_genes_TPM_greater_than_2 = as.integer(SL_counts)
)

PANC_results <- data.frame(
  Tissue = "PANC",
  Sample = names(PANC_counts),
  Expressed_genes_TPM_greater_than_2 = as.integer(PANC_counts)
)

LIV_results <- data.frame(
  Tissue = "LIV",
  Sample = names(LIV_counts),
  Expressed_genes_TPM_greater_than_2 = as.integer(LIV_counts)
)

all_results <- rbind(
  PAR_results,
  SM_results,
  SL_results,
  PANC_results,
  LIV_results
)


# Display and save the results -----------------------------------------------

print(all_results, row.names = FALSE)

write.csv(
  all_results,
  "normalization_checks/number_of_genes_expressed_per_gland/expressed_gene_counts_per_sample_TPM_greater_than_2.csv",
  row.names = FALSE
)
