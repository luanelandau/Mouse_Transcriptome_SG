#!/usr/bin/env Rscript

## ============================================
## Build per-gland mastersheets + secreted flag
## (uses NEW annotated dictionary) + mean TPMs
## ============================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
})

## ----------------
## File paths
## ----------------
tpm_fp  <- "~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/miscelaneous_sheets/mouse_expression/gene_expression_matrix_C57_CD1_TPMs.csv"

# Use the new dictionary you generated (with manual annotations + functions)
dict_fp <- "~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/miscelaneous_sheets/mouse_proteome/Geneid_secretome_dictionary.csv"

out_dir <- "~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/miscelaneous_sheets/mouse_expression"

## ----------------
## Load data
## ----------------
tpm <- readr::read_csv(tpm_fp, show_col_types = FALSE)

## Drop contaminated PAR columns (provided as dot-style; convert to hyphenated)
drops_dot <- c(
  "Mous.1A.PAR.Mal.L","Mous.2A.PAR.Mal.L",
  "Mous.6A.PAR.Fem.L","Mous.4B.PAR.Fem.L",
  "Mous.5B.PAR.Fem.L","Mous.6B.PAR.Fem.L"
)
drops_hyphen <- gsub("\\.", "-", drops_dot)

to_drop <- intersect(names(tpm), drops_hyphen)
if (length(to_drop) > 0) {
  tpm <- tpm %>% select(-all_of(to_drop))
  message("🧹 Dropped contaminated columns: ", paste(to_drop, collapse = ", "))
} else {
  message("ℹ️ No contaminated columns found to drop (after hyphen conversion).")
}

# Read dictionary as character to avoid parsing warnings
dict <- readr::read_csv(
  dict_fp,
  col_types = readr::cols(.default = readr::col_character())
)

# Keep required + optional columns
req <- c("Geneid", "secreted")
if (!all(req %in% names(dict))) {
  stop("Dictionary must have columns 'Geneid' and 'secreted'.")
}
opt <- intersect(c("Reference", "Function"), names(dict))
dict_min <- dict %>%
  select(all_of(c(req, opt))) %>%
  distinct()

## ----------------
## Gland patterns
## ----------------
glands   <- c("PAR", "SM", "SL", "LIV", "PANC")
patterns <- setNames(paste0("-", glands, "-"), glands)  # matches e.g. "-PAR-"

## ----------------
## Helpers
## ----------------
# Safe row mean that returns NA if no columns were provided
safe_row_mean <- function(df, cols) {
  if (length(cols) == 0) return(rep(NA_real_, nrow(df)))
  out <- suppressWarnings(rowMeans(as.matrix(df[, cols, drop = FALSE]), na.rm = TRUE))
  out[is.nan(out)] <- NA_real_
  out
}

## ----------------
## Build one mastersheet
## ----------------
make_mastersheet <- function(gland_code, pattern_token) {
  # Identify columns: Geneid + any with the token (e.g., "-PAR-")
  col_keep <- names(tpm) == "Geneid" | str_detect(names(tpm), fixed(pattern_token))
  # Explicitly exclude coordinates/length if present
  drop_cols <- c("Chr", "Start", "End", "Strand", "Length")
  col_keep[names(tpm) %in% drop_cols] <- FALSE
  
  # Subset TPMs
  sub <- tpm[, col_keep, drop = FALSE]
  
  # Join dictionary
  sub_annot <- sub %>%
    left_join(dict_min, by = "Geneid")
  
  # Identify sample columns (TPM columns for this gland)
  front_cols <- c("Geneid", "secreted", opt)
  sample_cols <- setdiff(names(sub_annot), front_cols)
  
  # Split by sex using the column names.
  # Match "-Mal-" or terminal "Mal" (same for Fem).
  sample_cols_male <- sample_cols[str_detect(sample_cols, "(^|-)Mal(-|$)")]
  sample_cols_fem  <- sample_cols[str_detect(sample_cols, "(^|-)Fem(-|$)")]
  
  # Compute per-row means
  sub_annot <- sub_annot %>%
    mutate(
      Mean_TPM      = safe_row_mean(., sample_cols),
      Mean_TPM_Male = safe_row_mean(., sample_cols_male),
      Mean_TPM_Fem  = safe_row_mean(., sample_cols_fem)
    )
  
  # Order columns: Geneid, secreted, Reference, Function, mean columns, then sample columns
  sub_annot <- sub_annot %>%
    select(all_of(front_cols), Mean_TPM, Mean_TPM_Male, Mean_TPM_Fem, all_of(sample_cols))
  
  # Write
  out_fp <- file.path(out_dir, paste0(gland_code, "_mastersheet_TPMs_annotated_for_secretion.csv"))
  readr::write_csv(sub_annot, out_fp)
  message("✅ Wrote: ", out_fp,
          " | rows: ", nrow(sub_annot),
          " | samples: ", length(sample_cols),
          " | male cols: ", length(sample_cols_male),
          " | female cols: ", length(sample_cols_fem))
}

## ----------------
## Build all mastersheets
## ----------------
walk2(names(patterns), unname(patterns), make_mastersheet)
