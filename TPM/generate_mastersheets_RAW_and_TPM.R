## ================================================================
## Build RAW COUNTS and TPM matrices, cleaning metadata to first part
## ================================================================
setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
})

dir_C57 <- "readCounts_C57"
dir_CD1 <- "readCounts_CD1"

files_counts <- c(
  list.files(dir_C57, pattern = "\\.counts\\.txt$", full.names = TRUE),
  list.files(dir_CD1, pattern = "\\.counts\\.txt$", full.names = TRUE)
)
files_counts <- files_counts[!grepl("counts_with_tpm\\.txt$", files_counts)]  # extra guard

files_counts_with_tpm <- c(
  list.files(dir_C57, pattern = "\\.counts_with_tpm\\.txt$", full.names = TRUE),
  list.files(dir_CD1, pattern = "\\.counts_with_tpm\\.txt$", full.names = TRUE)
)

stopifnot(length(files_counts) > 0, length(files_counts_with_tpm) > 0)

meta_cols <- c("Geneid","Chr","Start","End","Strand","Length")

# ---- helpers ----
clean_sample <- function(path) {
  b <- basename(path)
  b <- sub("\\.counts_with_tpm\\.txt$", "", b)
  b <- sub("\\.counts\\.txt$", "", b)
  b
}
read_tsv_safe <- function(f) {
  df <- read_tsv(
    f,
    col_types = cols(.default = col_character()),
    comment   = "#",          # <-- ignore the Program: line
    trim_ws   = TRUE,
    progress  = FALSE
  )
  keep <- !map_lgl(df, ~ all(is.na(.x) | .x == ""))
  df[, keep, drop = FALSE]
}
first_part <- function(x) sub(";.*$", "", as.character(x))

# Clean meta to the FIRST entry before any ';'
clean_meta <- function(df) {
  df %>%
    mutate(
      Chr    = first_part(Chr),
      Start  = suppressWarnings(as.integer(first_part(Start))),
      End    = suppressWarnings(as.integer(first_part(End))),
      Strand = first_part(Strand),
      Length = suppressWarnings(as.integer(Length))
    )
}

# Canonical metadata from one file (after cleaning)
canonical_meta <- function(f) {
  read_tsv_safe(f) %>% clean_meta() %>% select(all_of(meta_cols)) %>% distinct(Geneid, .keep_all = TRUE)
}

# (Geneid + sample) for RAW counts (from *.counts.txt)
gene_sample_counts <- function(f) {
  df <- read_tsv_safe(f) %>% clean_meta()
  samp <- clean_sample(f)
  counts_col <- setdiff(names(df), meta_cols)
  counts_col <- counts_col[length(counts_col)]  # last non-meta col = counts
  tibble(Geneid = df$Geneid, !!samp := suppressWarnings(as.numeric(df[[counts_col]])))
}

# (Geneid + sample) for TPM (from *.counts_with_tpm.txt)
gene_sample_tpm <- function(f) {
  df <- read_tsv_safe(f) %>% clean_meta()
  samp <- clean_sample(f)
  non_meta <- setdiff(names(df), meta_cols)
  tpm_col <- non_meta[str_detect(non_meta, "(?i)\\bTPM\\b|_TPM$|TPM$")]
  if (length(tpm_col) == 0) tpm_col <- non_meta[length(non_meta)]  # assume last if unnamed
  tpm_col <- tpm_col[1]
  tibble(Geneid = df$Geneid, !!samp := suppressWarnings(as.numeric(df[[tpm_col]])))
}

# Reduce (by Geneid) and attach canonical metadata up front
wide_by_gene <- function(lst_gene_sample, meta_tbl) {
  w <- purrr::reduce(lst_gene_sample, ~ dplyr::full_join(.x, .y, by = "Geneid"))
  meta_tbl %>%
    dplyr::left_join(w, by = "Geneid") %>%
    dplyr::select(dplyr::all_of(meta_cols), dplyr::everything()) %>%
    dplyr::arrange(Geneid)
}

## ---------------- RAW COUNTS ----------------
meta_counts <- canonical_meta(files_counts[1])
lst_counts  <- map(files_counts, gene_sample_counts)
wide_counts <- wide_by_gene(lst_counts, meta_counts)

## ---------------- TPM ----------------
meta_tpm <- canonical_meta(files_counts_with_tpm[1])
lst_tpm  <- map(files_counts_with_tpm, gene_sample_tpm)
wide_tpm <- wide_by_gene(lst_tpm, meta_tpm)

## ---------------- Write ----------------
out_counts <- "miscelaneous_sheets/gene_expression_matrix_C57_CD1_RAW_COUNTS.csv"
out_tpm    <- "miscelaneous_sheets/gene_expression_matrix_C57_CD1_TPMs.csv"
write_csv(wide_counts, out_counts)
write_csv(wide_tpm,    out_tpm)

cat("Done.\nRAW -> ", out_counts, "\nTPM -> ", out_tpm, "\n", sep = "")




