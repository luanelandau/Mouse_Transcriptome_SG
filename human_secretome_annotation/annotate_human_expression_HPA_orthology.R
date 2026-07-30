#!/usr/bin/env Rscript

# Human salivary-gland expression annotation
#
# This script adds Human Protein Atlas (HPA) annotations and human–mouse
# orthology to a gene-level TPM matrix, calculates mean TPM within each gland,
# and writes separate tables for parotid (PAR), sublingual (SL), and
# submandibular (SM) glands.
#
# Required R packages: dplyr and readr
#
# INPUTS GENERATED IN THIS STUDY
#
# 1. miscelaneous_sheets/gene_expression_matrix_human_TPMs.csv
#    Gene-level human TPM matrix. Geneid contains HGNC gene symbols and the
#    sample columns are named adult_PAR_*, adult_SL_*, and adult_SM_*.
#
# 2. orthology_assignment/orthologs_with_classification.csv
#    Human–mouse orthology assignments. Required columns are human_gene,
#    mouse_gene, and ortholog_type.
#
# EXTERNAL INPUTS
#
# 3. proteinatlas.tsv
#    Full gene-level table downloaded from the Human Protein Atlas (HPA)
#    version 25. The columns "Secretome location" and "Secretome function"
#    are used here.
#
# 4. normal_ihc_data.tsv
#    Human Protein Atlas v24.1 normal-tissue immunohistochemistry data,
#    downloaded as normal_ihc_data.tsv.zip from:
#    https://v24.proteinatlas.org/humanproteome/tissue/data
#
#    After decompression, rows with Tissue == "Salivary gland" were used
#    to annotate protein-detection Level, Reliability, and Cell type.
#
# HPA data are not distributed with this repository because of file size.
# Download the Protein Atlas and normal-tissue IHC tables from:
# https://www.proteinatlas.org/about/download
#
# Place both decompressed TSV files in a local directory and update
# hpa_data_dir below. Run this script from the Mouse_Transcriptome_SG
# repository root.

library(dplyr)
library(readr)


# -------------------------------------------------------------------------
# File locations
# -------------------------------------------------------------------------

# Directory containing the two downloaded HPA files.
hpa_data_dir <- path.expand(
  "~/Library/CloudStorage/Box-Box/SalivaryGlands_LL"
)

# Files included in this repository.
expression_file <- file.path(
  "miscelaneous_sheets",
  "gene_expression_matrix_human_TPMs.csv"
)

orthology_file <- file.path(
  "orthology_assignment",
  "orthologs_with_classification.csv"
)

# Large external files downloaded from the HPA.
protein_atlas_file <- file.path(
  hpa_data_dir,
  "proteinatlas.tsv"
)

normal_ihc_file <- file.path(
  hpa_data_dir,
  "normal_ihc_data.tsv"
)

# Output directory and filenames.
output_folder <- file.path(
  "miscelaneous_sheets",
  "human expression"
)

PAR_output_file <- file.path(
  output_folder,
  "PAR_human_mastersheet_TPMs_annotated_with_orthology.csv"
)

SL_output_file <- file.path(
  output_folder,
  "SL_human_mastersheet_TPMs_annotated_with_orthology.csv"
)

SM_output_file <- file.path(
  output_folder,
  "SM_human_mastersheet_TPMs_annotated_with_orthology.csv"
)

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)


# -------------------------------------------------------------------------
# Read gene expression
# -------------------------------------------------------------------------

human_expression <- read_csv(
  expression_file,
  show_col_types = FALSE
)


# -------------------------------------------------------------------------
# Add HPA Secretome location and function
#
# All non-missing HPA Secretome annotations are retained. This includes
# entries classified by HPA as "Intracellular and membrane"; no additional
# predicted-secretome filter is applied at this stage.
# -------------------------------------------------------------------------

protein_atlas <- read_tsv(
  protein_atlas_file,
  show_col_types = FALSE
)

protein_atlas_annotation <- protein_atlas %>%
  distinct(Gene, .keep_all = TRUE) %>%
  select(
    Geneid = Gene,
    `Secretome location`,
    `Secretome function`
  )


# -------------------------------------------------------------------------
# Add salivary-gland IHC annotations
# -------------------------------------------------------------------------

normal_ihc <- read_tsv(
  normal_ihc_file,
  show_col_types = FALSE
)

salivary_gland_ihc <- normal_ihc %>%
  filter(Tissue == "Salivary gland") %>%
  group_by(`Gene name`) %>%
  summarise(
    Level = paste(
      unique(Level[!is.na(Level)]),
      collapse = "; "
    ),
    Reliability = paste(
      unique(Reliability[!is.na(Reliability)]),
      collapse = "; "
    ),
    `Cell type` = paste(
      unique(`Cell type`[!is.na(`Cell type`)]),
      collapse = "; "
    ),
    .groups = "drop"
  ) %>%
  rename(Geneid = `Gene name`)

# Convert empty strings produced for missing IHC annotations to NA.
salivary_gland_ihc$Level[
  salivary_gland_ihc$Level == ""
] <- NA

salivary_gland_ihc$Reliability[
  salivary_gland_ihc$Reliability == ""
] <- NA

salivary_gland_ihc$`Cell type`[
  salivary_gland_ihc$`Cell type` == ""
] <- NA

human_annotated <- human_expression %>%
  left_join(
    protein_atlas_annotation,
    by = "Geneid"
  ) %>%
  left_join(
    salivary_gland_ihc,
    by = "Geneid"
  )


# -------------------------------------------------------------------------
# Add human–mouse orthology
#
# When a human gene has multiple mouse orthologs, their symbols are retained
# in one row and separated by " | ". Genes without a mouse ortholog are
# classified as human-specific.
# -------------------------------------------------------------------------

orthology <- read_csv(
  orthology_file,
  show_col_types = FALSE
)

orthology <- orthology %>%
  group_by(human_gene) %>%
  summarise(
    mouse_gene = paste(
      unique(mouse_gene[!is.na(mouse_gene)]),
      collapse = " | "
    ),
    ortholog_type = paste(
      unique(ortholog_type[!is.na(ortholog_type)]),
      collapse = " | "
    ),
    .groups = "drop"
  )

orthology$mouse_gene[
  orthology$mouse_gene == ""
] <- NA

orthology$ortholog_type[
  orthology$ortholog_type == ""
] <- NA

human_annotated <- human_annotated %>%
  left_join(
    orthology,
    by = c("Geneid" = "human_gene")
  ) %>%
  mutate(
    ortholog_type = if_else(
      is.na(mouse_gene),
      "human-specific",
      ortholog_type
    )
  )


# -------------------------------------------------------------------------
# Parotid gland
# -------------------------------------------------------------------------

PAR_annotated <- human_annotated %>%
  mutate(
    Mean_TPM = rowMeans(
      across(
        c(
          adult_PAR_01,
          adult_PAR_02,
          adult_PAR_03,
          adult_PAR_04
        )
      ),
      na.rm = TRUE
    )
  ) %>%
  select(
    Geneid,
    mouse_gene,
    ortholog_type,
    `Secretome location`,
    `Secretome function`,
    Level,
    Reliability,
    `Cell type`,
    Mean_TPM,
    adult_PAR_01,
    adult_PAR_02,
    adult_PAR_03,
    adult_PAR_04
  )

write_csv(
  PAR_annotated,
  PAR_output_file,
  na = "NA"
)


# -------------------------------------------------------------------------
# Sublingual gland
# -------------------------------------------------------------------------

SL_annotated <- human_annotated %>%
  mutate(
    Mean_TPM = rowMeans(
      across(
        c(
          adult_SL_01,
          adult_SL_02,
          adult_SL_03
        )
      ),
      na.rm = TRUE
    )
  ) %>%
  select(
    Geneid,
    mouse_gene,
    ortholog_type,
    `Secretome location`,
    `Secretome function`,
    Level,
    Reliability,
    `Cell type`,
    Mean_TPM,
    adult_SL_01,
    adult_SL_02,
    adult_SL_03
  )

write_csv(
  SL_annotated,
  SL_output_file,
  na = "NA"
)


# -------------------------------------------------------------------------
# Submandibular gland
# -------------------------------------------------------------------------

SM_annotated <- human_annotated %>%
  mutate(
    Mean_TPM = rowMeans(
      across(
        c(
          adult_SM_01,
          adult_SM_02,
          adult_SM_03,
          adult_SM_04,
          adult_SM_05,
          adult_SM_06
        )
      ),
      na.rm = TRUE
    )
  ) %>%
  select(
    Geneid,
    mouse_gene,
    ortholog_type,
    `Secretome location`,
    `Secretome function`,
    Level,
    Reliability,
    `Cell type`,
    Mean_TPM,
    adult_SM_01,
    adult_SM_02,
    adult_SM_03,
    adult_SM_04,
    adult_SM_05,
    adult_SM_06
  )

write_csv(
  SM_annotated,
  SM_output_file,
  na = "NA"
)


# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------

cat("\nHuman expression annotation completed.\n")
cat("Genes in expression matrix:", nrow(human_expression), "\n")
cat(
  "Genes with HPA Secretome location:",
  sum(!is.na(human_annotated$`Secretome location`)),
  "\n"
)
cat(
  "Genes with salivary-gland IHC data:",
  sum(!is.na(human_annotated$Level)),
  "\n"
)
cat("PAR output:", PAR_output_file, "\n")
cat("SL output:", SL_output_file, "\n")
cat("SM output:", SM_output_file, "\n")
