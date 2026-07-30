#!/usr/bin/env Rscript

# Convert GTEx GCT v1.2 files to CSV, retain PAR genes, calculate mean TPM
# across samples, and attach the human/mouse orthology annotations.
#
# Usage:
#   Rscript HUM_pancreas_liver_annotating_TPMs_github.R [data_directory]
#
# If data_directory is omitted, the script uses the Box directory below.
# Pass a different directory when running the script on another computer.
# The directory must contain:
#   - gene_tpm_v10_pancreas.gct
#   - gene_tpm_v10_liver.gct
#   - PAR_human_mastersheet_TPMs_annotated_with_orthology.csv

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required. Install it with install.packages('data.table').")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required. Install it with install.packages('dplyr').")
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) {
  stop("Usage: Rscript HUM_pancreas_liver_annotating_TPMs_github.R [data_directory]")
}

default_data_dir <- paste0(
  path.expand("~/Library/CloudStorage/Box-Box/"),
  "SalivaryGlands_LL/Mouse_Transcriptome_SG/",
  "miscelaneous_sheets/human expression"
)
data_dir <- if (length(args) == 1L) args[[1L]] else default_data_dir
data_dir <- normalizePath(data_dir, mustWork = TRUE)

required_files <- c(
  "gene_tpm_v10_pancreas.gct",
  "gene_tpm_v10_liver.gct",
  "PAR_human_mastersheet_TPMs_annotated_with_orthology.csv"
)
missing_files <- required_files[!file.exists(file.path(data_dir, required_files))]
if (length(missing_files) > 0L) {
  stop("Missing required file(s) in ", data_dir, ": ", paste(missing_files, collapse = ", "))
}

convert_gct_v12_to_csv <- function(path) {
  message("Converting: ", basename(path))
  header <- readLines(path, n = 2L, warn = FALSE)
  if (length(header) < 2L || !grepl("^#?\\s*1\\.2\\s*$", header[[1L]])) {
    stop("Expected a GCT v1.2 file: ", path)
  }

  dimensions <- strsplit(header[[2L]], "\t", fixed = TRUE)[[1L]]
  if (length(dimensions) != 2L || any(!grepl("^\\d+$", dimensions))) {
    stop("Invalid GCT dimension line in: ", path)
  }

  table <- fread(path, skip = 2L, sep = "\t", header = TRUE, data.table = TRUE)
  expected_rows <- as.integer(dimensions[[1L]])
  expected_samples <- as.integer(dimensions[[2L]])
  if (nrow(table) != expected_rows || ncol(table) != expected_samples + 2L) {
    stop(
      "GCT dimensions do not match the data in ", basename(path),
      ": expected ", expected_rows, " rows and ", expected_samples,
      " samples; found ", nrow(table), " rows and ", ncol(table) - 2L, " samples."
    )
  }
  if (!all(c("Name", "Description") %in% names(table))) {
    stop("GCT must contain Name and Description columns: ", path)
  }

  output_path <- sub("\\.gct$", ".csv", path, ignore.case = TRUE)
  fwrite(table, output_path)
  message("  -> wrote: ", output_path)
  output_path
}

annotation_path <- file.path(
  data_dir,
  "PAR_human_mastersheet_TPMs_annotated_with_orthology.csv"
)
annotation_full <- fread(annotation_path)
annotation_columns <- c(
  "Geneid", "mouse_gene", "ortholog_type",
  "Secretome location", "Secretome function",
  "Level", "Reliability", "Cell type"
)
missing_annotation_columns <- setdiff(annotation_columns, names(annotation_full))
if (length(missing_annotation_columns) > 0L) {
  stop(
    "Annotation file is missing column(s): ",
    paste(missing_annotation_columns, collapse = ", ")
  )
}

annotation <- annotation_full |>
  select(all_of(annotation_columns)) |>
  distinct()
par_genes <- unique(annotation$Geneid)

process_organ <- function(gct_filename, organ_tag) {
  csv_path <- convert_gct_v12_to_csv(file.path(data_dir, gct_filename))
  expression <- fread(csv_path)
  sample_columns <- setdiff(names(expression), c("Name", "Description"))
  nonnumeric_samples <- sample_columns[
    !vapply(expression[, ..sample_columns], is.numeric, logical(1L))
  ]
  if (length(nonnumeric_samples) > 0L) {
    stop(
      "Non-numeric GTEx sample column(s) in ", basename(csv_path), ": ",
      paste(nonnumeric_samples, collapse = ", ")
    )
  }

  filtered <- expression[Description %in% par_genes]
  raw_output <- sub(
    "\\.csv$", "_filteredBy_PARgenes_raw.csv", csv_path, ignore.case = TRUE
  )
  fwrite(filtered, raw_output)
  message("  -> wrote raw: ", raw_output, " (", nrow(filtered), " rows)")

  filtered[, Mean_TPM := rowMeans(.SD, na.rm = TRUE), .SDcols = sample_columns]
  means <- filtered[, .(Mean_TPM = mean(Mean_TPM, na.rm = TRUE)), by = Description]

  annotated_means <- means |>
    inner_join(annotation, by = c("Description" = "Geneid")) |>
    transmute(
      Geneid = Description,
      mouse_gene,
      ortholog_type,
      `Secretome location`,
      `Secretome function`,
      Level,
      Reliability,
      `Cell type`,
      Mean_TPM
    ) |>
    mutate(
      across(
        c(
          mouse_gene, `Secretome location`, `Secretome function`,
          Level, Reliability, `Cell type`
        ),
        ~ ifelse(ortholog_type == "human-specific", NA, .x)
      )
    )

  means_output <- sub(
    "\\.csv$",
    "_filteredBy_PARgenes_annotated_means.csv",
    csv_path,
    ignore.case = TRUE
  )
  fwrite(annotated_means, means_output)
  message(
    "  -> wrote means: ", means_output,
    " (", nrow(annotated_means), " annotated rows)"
  )

  data.frame(
    organ = organ_tag,
    total_symbols = dplyr::n_distinct(expression$Description),
    matched_symbols = dplyr::n_distinct(annotated_means$Geneid),
    pct_matched = 100 * dplyr::n_distinct(annotated_means$Geneid) /
      dplyr::n_distinct(expression$Description)
  )
}

stats <- rbind(
  process_organ("gene_tpm_v10_pancreas.gct", "pancreas"),
  process_organ("gene_tpm_v10_liver.gct", "liver")
)
print(stats, row.names = FALSE)
