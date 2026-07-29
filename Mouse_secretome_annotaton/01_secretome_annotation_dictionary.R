#!/usr/bin/env Rscript

## ===========================================
## Build a Geneid → secretome annotation dictionary
## (with manual overrides + function-only adds)
## ===========================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

## -------------
## File paths
## -------------
#raw counts for all samples
counts_fp    <- "~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/miscelaneous_sheets/mouse_expression/gene_expression_matrix_C57_CD1_RAW_COUNTS.csv"

#secreted_any is taken from secreted.tsv, which includes all annotations for all proteins from uniprot, but filtered for all secreted genes
secreted_fp  <- "~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/miscelaneous_sheets/mouse_proteome/secreted_any.tsv"

# Manual curation (secreted + refs/functions)
#For manual curations and function add additions, I have added manually some genes that I found to be secreted, but are not annotated by as secreted in this dataset. 
#I have included the annotation and the reference from where I took this into the sheet. 
manual_fp    <- "~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/miscelaneous_sheets/mouse_proteome/mannual_curation.csv"
# Function-only additions
function_add_fp <- "~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/miscelaneous_sheets/mouse_proteome/function_add.csv"

# Output dictionary (CSV + RDS)
out_dir      <- "~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/miscelaneous_sheets/mouse_proteome"
out_csv      <- file.path(out_dir, "Geneid_secretome_dictionary.csv")
out_rds      <- file.path(out_dir, "Geneid_secretome_dictionary.rds")

## -------------------------
## Load counts and secretome
## -------------------------
counts <- readr::read_csv(counts_fp, show_col_types = FALSE) %>%
  select(Geneid) %>%
  distinct() %>%
  mutate(Geneid_norm = toupper(Geneid))

secretome <- readr::read_tsv(secreted_fp, col_types = cols(.default = "c"))

if (!"Gene Names" %in% names(secretome)) {
  stop("Expected column 'Gene Names' not found in secretome TSV.")
}
if (!"Entry" %in% names(secretome)) {
  stop("Expected column 'Entry' not found in secretome TSV.")
}

## -------------------------------------------------
## Create alias map: one row per gene alias per entry
## -------------------------------------------------
alias_map <- secretome %>%
  dplyr::mutate(
    `Gene Names` = dplyr::if_else(
      is.na(`Gene Names`),
      "",
      `Gene Names`
    ),
    alias_vec = stringr::str_split(`Gene Names`, "\\s+")
  ) %>%
  tidyr::unnest(
    cols = alias_vec,
    keep_empty = TRUE
  ) %>%
  dplyr::rename(
    gene_alias = alias_vec
  ) %>%
  dplyr::filter(
    !is.na(gene_alias),
    gene_alias != ""
  ) %>%
  dplyr::mutate(
    gene_alias_norm = toupper(gene_alias)
  )
## -------------------------------------------------------------
## Collapse to per-alias summary (aggregating ALL secretome cols)
## -------------------------------------------------------------
collapse_unique <- function(x) {
  x <- unique(na.omit(x))
  x <- x[nzchar(x)]
  if (length(x) == 0) NA_character_ else paste(x, collapse = " | ")
}

secretome_cols <- colnames(secretome)

per_alias_summary <- alias_map %>%
  group_by(gene_alias_norm) %>%
  summarise(
    across(all_of(secretome_cols), collapse_unique),
    .groups = "drop"
  ) %>%
  mutate(secreted = "yes")

## ---------------------------------------------------
## Join onto the Geneid universe and finalize dictionary
## ---------------------------------------------------
dictionary <- counts %>%
  left_join(per_alias_summary, by = c("Geneid_norm" = "gene_alias_norm")) %>%
  select(-Geneid_norm)

# Mark non-matches
secretome_cols_present <- intersect(secretome_cols, names(dictionary))

dictionary <- dictionary %>%
  mutate(
    secreted = if_else(is.na(secreted), "No_annotation", secreted)
  )

# For genes with NOT "yes", blank the secretome columns to NA
dictionary <- dictionary %>%
  mutate(across(all_of(secretome_cols_present),
                ~ if_else(secreted == "yes", .x, NA_character_)))

## -----------------------
## Apply manual overrides
## -----------------------
manual_raw <- readr::read_csv(
  manual_fp,
  col_types = readr::cols(.default = readr::col_character())
)

manual <- manual_raw %>%
  transmute(
    Geneid       = trimws(Geneid),
    Geneid_norm  = toupper(Geneid),
    secreted     = trimws(coalesce(secreted, NA_character_)),
    Reference    = trimws(coalesce(Reference, NA_character_)),
    Function     = trimws(coalesce(Function, NA_character_))
  ) %>%
  filter(!is.na(Geneid_norm), Geneid_norm != "") %>%
  group_by(Geneid_norm) %>%
  summarise(
    secreted  = collapse_unique(secreted),
    Reference = collapse_unique(Reference),
    Function  = collapse_unique(Function),
    .groups = "drop"
  )

# Ensure columns exist
if (!"Reference" %in% names(dictionary)) dictionary$Reference <- NA_character_
if (!"Function"  %in% names(dictionary)) dictionary$Function  <- NA_character_

# Join and override from manual
dictionary <- dictionary %>%
  mutate(Geneid_norm = toupper(Geneid)) %>%
  left_join(manual, by = "Geneid_norm", suffix = c("", ".manual")) %>%
  mutate(
    secreted  = if_else(!is.na(secreted.manual) & nzchar(secreted.manual),
                        secreted.manual, secreted),
    Reference = coalesce(Reference.manual, Reference),
    Function  = coalesce(Function.manual, Function)
  ) %>%
  select(-ends_with(".manual"))

# Re-blank secretome columns for anything not "yes"
dictionary <- dictionary %>%
  mutate(across(all_of(secretome_cols_present),
                ~ if_else(secreted == "yes", .x, NA_character_)))

## --------------------------------------------
## Apply function-only additions (no other edits)
## --------------------------------------------
func_add_raw <- readr::read_csv(
  function_add_fp,
  col_types = readr::cols(.default = readr::col_character())
)

func_add <- func_add_raw %>%
  transmute(
    Geneid      = trimws(Geneid),
    Geneid_norm = toupper(Geneid),
    Function_add = trimws(coalesce(Function, NA_character_))
  ) %>%
  filter(!is.na(Geneid_norm), Geneid_norm != "") %>%
  group_by(Geneid_norm) %>%
  summarise(Function_add = collapse_unique(Function_add), .groups = "drop")

dictionary <- dictionary %>%
  left_join(func_add, by = "Geneid_norm") %>%
  mutate(
    # only fill Function where it is currently NA
    Function = if_else(is.na(Function) | !nzchar(Function),
                       Function_add, Function)
  ) %>%
  select(-Function_add, -Geneid_norm)

## ----------------
## Write outputs
## ----------------
readr::write_csv(dictionary, out_csv)
saveRDS(dictionary, out_rds)

## -------------
## Quick summary
## -------------
n_total  <- nrow(dictionary)
n_yes    <- sum(dictionary$secreted == "yes", na.rm = TRUE)
n_na     <- sum(dictionary$secreted == "No_annotation", na.rm = TRUE)
n_manual <- sum(dictionary$secreted == "manual_annotation", na.rm = TRUE)
n_func_filled <- sum(!is.na(dictionary$Function) & nzchar(dictionary$Function))

message("✅ Dictionary written:")
message("  CSV: ", out_csv)
message("  RDS: ", out_rds)
message("Counts: total Geneid = ", n_total,
        " | yes = ", n_yes,
        " | manual_annotation = ", n_manual,
        " | no annotation = ", n_na)
message("Function filled (non-empty): ", n_func_filled)
