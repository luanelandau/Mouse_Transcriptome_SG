# =============================== #
#   FULL PIPELINE: Per-gland DESeq2 + Master Sheet Generation
# =============================== #

setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/")

# --- Libraries ---
library(DESeq2)
library(tidyverse)
library(readr)
library(tidyr)
library(stringr)
library(purrr)

# --- Load and clean count data (raw integer counts required) ---
dat <- read.csv("miscelaneous_sheets/gene_expression_matrix_C57_CD1.csv")

# Drop metadata and contaminated samples (fix syntax)
drops <- c(
  "Chr","Start","End","Strand","Length",
  "Mous.1A.PAR.Mal.L","Mous.2A.PAR.Mal.L",
  "Mous.6A.PAR.Fem.L","Mous.4B.PAR.Fem.L",
  "Mous.5B.PAR.Fem.L","Mous.6B.PAR.Fem.L"
)
dat <- dat[, !(names(dat) %in% drops)]

# Create count matrix
count_data <- dat[, -1, drop = FALSE]
rownames(count_data) <- dat[, 1, drop = TRUE]

# --- Build metadata (info) ---
info <- tibble(
  ID = colnames(count_data)
) %>%
  mutate(
    sex = case_when(
      str_detect(ID, "Mal") ~ "male",
      str_detect(ID, "Fem") ~ "female",
      TRUE ~ NA_character_
    ),
    strain = case_when(
      str_detect(ID, "B") ~ "C57",
      str_detect(ID, "A") ~ "CD1",
      TRUE ~ NA_character_
    ),
    gland = str_extract(ID, "LIV|PANC|PAR|SL|SM")
  ) %>%
  mutate(
    # set references explicitly
    sex    = factor(sex,    levels = c("female","male")),
    strain = factor(strain, levels = c("CD1","C57")),
    gland  = factor(gland,  levels = c("PAR","SL","SM","LIV","PANC"))
  )

# Check alignment
stopifnot(all(info$ID == colnames(count_data)))

# --- Output folder ---
dir.create("deseq2/deseq2_results", showWarnings = FALSE)

# --- Per-gland DESeq2 with per-gland prefiltering ---
gland_list <- unique(info$gland)

for (g in gland_list) {
  message("\n========== Gland: ", g, " ==========")
  info_sub   <- dplyr::filter(info, gland == g)
  counts_sub <- count_data[, info_sub$ID, drop = FALSE]
  
  # Per-gland low-count filter: keep genes with >= 10 total counts within this gland
  keep <- rowSums(counts_sub) >= 10
  counts_sub <- counts_sub[keep, , drop = FALSE]
  
  # Safety: ensure integers (DESeq2 expects integers)
  counts_sub[] <- round(as.matrix(counts_sub))
  
  # Build DESeq2 object
  dds_sub <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData   = droplevels(as.data.frame(info_sub)),
    design    = ~ strain + sex
  )
  
  # Fit (poscounts helps with many zeros)
  dds_sub <- DESeq(dds_sub, sfType = "poscounts")
  
  # Main-effect contrasts
  res_strain <- results(dds_sub, contrast = c("strain","C57","CD1"))
  res_sex    <- results(dds_sub, contrast = c("sex","male","female"))
  
  # (Optional but recommended) LFC shrinkage for more stable LFCs
  # If apeglm is installed; otherwise comment these two lines out.
  if (requireNamespace("apeglm", quietly = TRUE)) {
    suppressPackageStartupMessages(library(apeglm))
    res_strain <- lfcShrink(dds_sub, contrast = c("strain","C57","CD1"), type = "apeglm")
    res_sex    <- lfcShrink(dds_sub, contrast = c("sex","male","female"), type = "apeglm")
  }
  
  # Save results
  write.csv(as.data.frame(res_strain),
            file = file.path("deseq2/deseq2_results", paste0("DESeq2_strain_", g, ".csv")),
            row.names = TRUE)
  write.csv(as.data.frame(res_sex),
            file = file.path("deseq2/deseq2_results", paste0("DESeq2_sex_", g, ".csv")),
            row.names = TRUE)
  
  # Quick summary
  cat("Strain results summary:\n"); print(summary(res_strain))
  cat("Sex results summary:\n");    print(summary(res_sex))
}

# =============================== #
#   Combine into master sheet
# =============================== #

result_folder <- "deseq2/deseq2_results"
glands       <- c("PAR","SL","SM","LIV","PANC")
comparisons  <- c("strain","sex")

all_results <- list()
for (g in glands) {
  for (cmp in comparisons) {
    file_path <- file.path(result_folder, paste0("DESeq2_", cmp, "_", g, ".csv"))
    stopifnot(file.exists(file_path))
    
    res <- read.csv(file_path, row.names = 1, stringsAsFactors = FALSE) %>%
      tibble::rownames_to_column(var = "gene") %>%
      dplyr::select(gene, baseMean, log2FoldChange, padj) %>%
      dplyr::mutate(comparison = cmp, gland = g)
    
    all_results[[paste0(cmp, "_", g)]] <- res
  }
}

master_df <- bind_rows(all_results)

master_wide <- master_df %>%
  pivot_wider(
    id_cols = gene,
    names_from  = c("comparison","gland"),
    values_from = c("log2FoldChange","padj","baseMean"),
    names_glue  = "{.value}_{comparison}_{gland}"
  )

# Order columns by blocks (metric x comparison x gland)
glands      <- c("PAR","SL","SM","LIV","PANC")
comparisons <- c("strain","sex")
metrics     <- c("log2FoldChange","padj","baseMean")

ordered_cols <- tidyr::expand_grid(metric = metrics, cmp = comparisons, gland = glands) %>%
  mutate(col = paste0(metric, "_", cmp, "_", gland)) %>%
  pull(col)

ordered_cols <- c("gene", intersect(ordered_cols, colnames(master_wide)))
master_wide  <- master_wide[, ordered_cols, drop = FALSE]

write_csv(master_wide, file.path(result_folder, "DESeq2_master_sheet_by_gene.csv"))

# =============================== #
#   Optional: save filtered DEG-only CSVs
# =============================== #

setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/deseq2_results")
if (!dir.exists("filtered")) dir.create("filtered")

all_files <- list.files(pattern = "DESeq2_(strain|sex)_.+\\.csv$")
filter_and_save <- function(file) {
  df <- read_csv(file, show_col_types = FALSE)
  logfc_col <- names(df)[str_detect(names(df), "log2FoldChange")]
  padj_col  <- names(df)[str_detect(names(df), "^padj$")]
  
  if (length(logfc_col) == 1 && length(padj_col) == 1) {
    df_filtered <- df %>%
      filter(!is.na(.data[[logfc_col]]), !is.na(.data[[padj_col]])) %>%
      filter(.data[[padj_col]] < 0.05, abs(.data[[logfc_col]]) > 1)
    
    out_name <- str_replace(file, "\\.csv$", "_filtered.csv")
    write_csv(df_filtered, file.path("filtered", out_name))
    message("Saved: ", file.path("filtered", out_name))
  } else {
    warning("Missing logFC or padj column in ", file)
  }
}
invisible(lapply(all_files, filter_and_save))

# =============================== #
#   Master sheet BY GENE (ANNOTATED)
#   - reads from deseq2_results/annotated
#   - keeps human_gene, mgi_id, ortholog_type
# =============================== #

library(tidyverse)
library(readr)
library(stringr)

result_folder <- "deseq2_results/annotated"   # <- annotated files live here
glands        <- c("PAR","SL","SM","LIV","PANC")
comparisons   <- c("strain","sex")

# collect all (cmp,gland) tables
all_results <- list()
for (g in glands) {
  for (cmp in comparisons) {
    file_path <- file.path(result_folder, paste0("DESeq2_", cmp, "_", g, ".csv"))
    stopifnot(file.exists(file_path))
    
    res <- readr::read_csv(file_path, show_col_types = FALSE) %>%
      # expected columns already include: gene, baseMean, log2FoldChange, padj, human_gene, mgi_id, ortholog_type
      dplyr::select(gene, human_gene, mgi_id, ortholog_type,
                    baseMean, log2FoldChange, padj) %>%
      dplyr::mutate(comparison = cmp, gland = g)
    
    all_results[[paste0(cmp, "_", g)]] <- res
  }
}

master_df <- dplyr::bind_rows(all_results)

# Wide master with annotation columns preserved as IDs
master_wide <- master_df %>%
  tidyr::pivot_wider(
    id_cols    = c(gene, human_gene, mgi_id, ortholog_type),
    names_from = c("comparison","gland"),
    values_from = c("log2FoldChange","padj","baseMean"),
    names_glue = "{.value}_{comparison}_{gland}"
  )

# Order columns by blocks (metric x comparison x gland)
glands      <- c("PAR","SL","SM","LIV","PANC")
comparisons <- c("strain","sex")
metrics     <- c("log2FoldChange","padj","baseMean")

ordered_cols <- tidyr::expand_grid(metric = metrics, cmp = comparisons, gland = glands) %>%
  dplyr::mutate(col = paste0(metric, "_", cmp, "_", gland)) %>%
  dplyr::pull(col)

# keep the 4 ID columns first, then the ordered metric blocks
id_cols <- c("gene","human_gene","mgi_id","ortholog_type")
master_wide <- master_wide[, c(id_cols, intersect(ordered_cols, colnames(master_wide))), drop = FALSE]

# write output next to inputs
readr::write_csv(master_wide, file.path(result_folder, "DESeq2_annotated_master_sheet_by_gene.csv"))

message("✅ Wrote: ", file.path(result_folder, "DESeq2_annotated_master_sheet_by_gene.csv"))
