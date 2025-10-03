#!/usr/bin/env Rscript

## ============================================================
## Add human orthology to existing per-gland mastersheets only
## (preserves one-to-many mappings; fills missing as mouse-specific)
## ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
})

## ----------------
## File paths
## ----------------
in_dir  <- "/Users/llandau/Library/CloudStorage/Box-Box/SalivaryGlands_LL/miscelaneous_sheets"
out_dir <- in_dir

# Target glands (we'll skip ones that don't exist)
glands <- c("PAR", "SM", "SL", "LIV", "PANC")

# Expected input/output names
in_files  <- setNames(file.path(in_dir,  paste0(glands, "_mastersheet_TPMs_annotated_for_secretion.csv")), glands)
out_files <- setNames(file.path(out_dir, paste0(glands, "_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv")), glands)

# Orthology dictionary
orth_fp <- file.path(in_dir, "orthologs_with_classification.csv")

## ----------------
## Load orthologs
## ----------------
orth <- readr::read_csv(
  orth_fp,
  col_types = readr::cols(.default = readr::col_character())
) %>%
  transmute(
    mouse_gene    = str_trim(mouse_gene),
    human_gene    = str_trim(human_gene),
    ortholog_type = str_trim(ortholog_type)
  )

if (!all(c("mouse_gene", "human_gene", "ortholog_type") %in% names(orth))) {
  stop("Ortholog file must contain columns: mouse_gene, human_gene, ortholog_type")
}

## ----------------
## Helpers
## ----------------
annotate_one <- function(g) {
  in_fp  <- in_files[[g]]
  out_fp <- out_files[[g]]
  
  if (!file.exists(in_fp)) {
    warning("Skipping ", g, " — file not found: ", in_fp)
    return(invisible(NULL))
  }
  
  df <- readr::read_csv(in_fp, show_col_types = FALSE)
  
  if (!"Geneid" %in% names(df)) {
    warning("Skipping ", g, " — 'Geneid' column missing in: ", in_fp)
    return(invisible(NULL))
  }
  
  # Left-join expands rows for one-to-many mappings automatically
  df2 <- df %>%
    left_join(orth, by = c("Geneid" = "mouse_gene")) %>%
    # If no human match, mark as mouse-specific
    mutate(
      ortholog_type = if_else(is.na(human_gene), "mouse-specific", ortholog_type)
    ) %>%
    relocate(human_gene, ortholog_type, .after = Geneid)
  
  readr::write_csv(df2, out_fp)
  message("✅ Wrote: ", out_fp,
          " | rows: ", nrow(df2),
          " | cols: ", ncol(df2))
}

## ----------------
## Run
## ----------------
existing <- glands[file.exists(in_files)]
missing  <- setdiff(glands, existing)

if (length(missing) > 0) {
  message("⚠️ Missing mastersheets (will skip): ", paste(missing, collapse = ", "))
  found <- list.files(in_dir, pattern = "_mastersheet_TPMs_annotated\\.csv$", full.names = TRUE)
  if (length(found)) {
    message("📄 Found in directory:\n  - ", paste(found, collapse = "\n  - "))
  }
}

walk(existing, annotate_one)

message("All done.")
