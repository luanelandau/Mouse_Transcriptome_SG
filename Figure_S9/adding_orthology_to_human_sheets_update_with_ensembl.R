#!/usr/bin/env Rscript

## ============================================================
## Update the three human gland mastersheets with Ensembl orthology
## - Replace the existing mouse_gene and ortholog_type columns
## - Collapse multiple mouse paralogs into one pipe-separated value
## - Mark genes without a mouse ortholog as "human-specific"
## ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

## ----------------
## File paths
## ----------------
base_dir <- "~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/Figure_S9/"
human_dir <- "~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/miscelaneous_sheets/human expression/"

in_fp_PAR <- file.path(human_dir, "PAR_human_mastersheet_TPMs_annotated_with_orthology.csv")
in_fp_SL  <- file.path(human_dir, "SL_human_mastersheet_TPMs_annotated_with_orthology.csv")
in_fp_SM  <- file.path(human_dir, "SM_human_mastersheet_TPMs_annotated_with_orthology.csv")

in_fp_orth <- file.path(base_dir, "orthologs_with_classification_ensembl.csv")
out_dir_human <- file.path(base_dir, "mastersheets_ensembl")

## ----------------
## Load the three gland-specific human mastersheets
## ----------------
human <- bind_rows(
  PAR = readr::read_csv(in_fp_PAR, show_col_types = FALSE),
  SL  = readr::read_csv(in_fp_SL, show_col_types = FALSE),
  SM  = readr::read_csv(in_fp_SM, show_col_types = FALSE),
  .id = "gland"
) %>%
  select(-mouse_gene, -ortholog_type)

## ----------------
## Load Ensembl orthology, join, and collapse per gland and human gene
## ----------------
orth <- readr::read_csv(
  in_fp_orth,
  col_types = readr::cols(.default = readr::col_character())
) %>%
  transmute(
    human_gene    = str_trim(human_gene),
    mouse_gene    = str_trim(mouse_gene),
    ortholog_type = str_trim(ortholog_type)
  )

human_anno <- human %>%
  left_join(
    orth,
    by = c("Geneid" = "human_gene"),
    relationship = "many-to-many"
  ) %>%
  group_by(gland, Geneid) %>%
  summarise(
    across(-c(mouse_gene, ortholog_type), first),
    mouse_gene = str_c(sort(unique(na.omit(mouse_gene))), collapse = " | "),
    ortholog_type = str_c(sort(unique(na.omit(ortholog_type))), collapse = " | "),
    .groups = "drop"
  ) %>%
  mutate(
    mouse_gene = na_if(mouse_gene, ""),
    ortholog_type = coalesce(na_if(ortholog_type, ""), "human-specific")
  ) %>%
  relocate(mouse_gene, ortholog_type, .after = Geneid)

## ----------------
## Write the three updated gland-specific mastersheets
## ----------------
PAR_human_anno <- human_anno %>%
  filter(gland == "PAR") %>%
  select(Geneid:Mean_TPM, starts_with("adult_PAR_"))

SL_human_anno <- human_anno %>%
  filter(gland == "SL") %>%
  select(Geneid:Mean_TPM, starts_with("adult_SL_"))

SM_human_anno <- human_anno %>%
  filter(gland == "SM") %>%
  select(Geneid:Mean_TPM, starts_with("adult_SM_"))

out_fp_PAR <- file.path(out_dir_human, "PAR_human_mastersheet_TPMs_annotated_with_orthology_ensembl.csv")
out_fp_SL  <- file.path(out_dir_human, "SL_human_mastersheet_TPMs_annotated_with_orthology_ensembl.csv")
out_fp_SM  <- file.path(out_dir_human, "SM_human_mastersheet_TPMs_annotated_with_orthology_ensembl.csv")

readr::write_csv(PAR_human_anno, out_fp_PAR)
readr::write_csv(SL_human_anno, out_fp_SL)
readr::write_csv(SM_human_anno, out_fp_SM)

message("Wrote: ", out_fp_PAR, " | rows: ", nrow(PAR_human_anno))
message("Wrote: ", out_fp_SL,  " | rows: ", nrow(SL_human_anno))
message("Wrote: ", out_fp_SM,  " | rows: ", nrow(SM_human_anno))
message("All done.")
