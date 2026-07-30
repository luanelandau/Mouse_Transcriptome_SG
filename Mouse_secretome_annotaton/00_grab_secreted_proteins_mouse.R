#!/usr/bin/env Rscript

## =======================================================
## Extract rows that mention "secreted" anywhere (inclusive)
## =======================================================
## Input : uniprotkb_proteome_UP000000589_2025_10_01.tsv -- This is the whole mouse proteome from UNIPROT
## Output: secreted_any.tsv (filtered to any row mentioning "secreted")
##         Plus a summary printed to console
##
## Notes:
## - Case-insensitive search for the literal substring "secreted".
## - Scans ALL columns; you’ll also get a column listing which fields matched.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
})

# ----------------
# File paths
# ----------------
#This file is no longer in the folder because of size. Please redownload it from UNIPROT
in_fp  <- "~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/miscelaneous_sheets/mouse_proteome/uniprotkb_proteome_UP000000589_2025_10_01.tsv"
out_fp <- "~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/miscelaneous_sheets/mouse_proteome/secreted_any.tsv"

# ----------------
# Load file (all char)
# ----------------
u <- readr::read_tsv(in_fp, col_types = cols(.default = "c"))

# ----------------
# Patterns to match
# ----------------
# You asked for anything that says "secreted".
# If you later want to broaden this, add terms like "secretion", "secretory", etc.
patterns <- c("secreted")  # <- add more terms here if desired

# precompile a single regex: secreted|secretory|... (case-insensitive)
rx <- str_c("(", str_c(patterns, collapse = "|"), ")", collapse = "")
rx <- regex(rx, ignore_case = TRUE)

# ----------------
# Scan all columns
# ----------------
# Make a logical matrix of matches per column
match_matrix <- map_lgl(u, ~ anyNA(.x)) # dummy init to get length
match_matrix <- NULL

col_hits <- lapply(u, function(col) {
  if (!is.character(col)) col <- as.character(col)
  str_detect(col, rx)
})

# Which rows have any "secreted" mention?
any_hit <- reduce(col_hits, `|`)

# Build a column listing which fields matched "secreted" for each row
hit_sources <- pmap_chr(
  as.list(data.frame(col_hits, check.names = FALSE)),
  function(...) {
    v <- c(...)
    if (all(is.na(v))) return(NA_character_)
    idx <- which(!is.na(v) & v)
    if (length(idx) == 0) return(NA_character_)
    paste(names(col_hits)[idx], collapse = " | ")
  }
)

# ----------------
# Assemble output
# ----------------
u_secreted <- u %>%
  mutate(.secreted_match_cols = hit_sources) %>%
  filter(any_hit)

# ----------------
# Write and summarize
# ----------------
readr::write_tsv(u_secreted, out_fp)

message("✅ Wrote: ", out_fp)
message("Rows total: ", nrow(u))
message("Rows with 'secreted' anywhere: ", nrow(u_secreted))
message("Examples of columns that matched (first 5 rows):")
print(head(u_secreted$.secreted_match_cols, 5))
