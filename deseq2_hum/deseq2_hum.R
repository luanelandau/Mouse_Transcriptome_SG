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
dat <- read.csv("miscelaneous_sheets/gene_expression_matrix_human_RAW_COUNTS.csv")

# Drop metadata and contaminated samples (fix syntax)
drops <- c(
  "Chr","Start","End","Strand","Length"
)
dat <- dat[, !(names(dat) %in% drops)]

# Create count matrix
count_data <- dat[, -1, drop = FALSE]
rownames(count_data) <- dat[, 1, drop = TRUE]

# --- Subset to SM-only samples ---
sm_cols <- grep("^adult_SM_", colnames(count_data), value = TRUE)
count_sm <- count_data[, sm_cols, drop = FALSE]

# --- Build metadata (info) only for SM samples ---
info <- tibble(
  ID = colnames(count_sm)
) %>%
  mutate(
    ## Sex: assign manually based on sample ID
    sex = case_when(
      ID %in% c("adult_SM_02", "adult_SM_04", "adult_SM_05") ~ "male",
      TRUE                                                   ~ "female"
    )
  ) %>%
  mutate(
    sex = factor(sex, levels = c("female","male"))
  )

# quick check
count_sm[1:5, 1:5]
info

# Check alignment
stopifnot(all(info$ID == colnames(count_data)))

# --- Output folder ---
#dir.create("deseq2_results_hum", showWarnings = FALSE)

# --- Per-gland DESeq2 with per-gland prefiltering ---
info_list <- unique(info$sex)

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
            file = file.path("deseq2_hum/", paste0("DESeq2_strain_", g, ".csv")),
            row.names = TRUE)
  write.csv(as.data.frame(res_sex),
            file = file.path("deseq2_hum/", paste0("DESeq2_sex_", g, ".csv")),
            row.names = TRUE)
  
  # Quick summary
  cat("Strain results summary:\n"); print(summary(res_strain))
  cat("Sex results summary:\n");    print(summary(res_sex))
}

# Assumes you already created:
#   count_sm  <- count_data[, grep("^adult_SM_", colnames(count_data)), drop = FALSE]
#   info      <- tibble(ID = colnames(count_sm)) %>%
#                 mutate(sex = if_else(ID %in% c("adult_SM_02","adult_SM_04","adult_SM_05"), "male", "female")) %>%
#                 mutate(sex = factor(sex, levels = c("female","male")))

library(DESeq2)

# Per-gland (SM) low-count filter: keep genes with >= 10 total counts in SM
keep <- rowSums(count_sm) >= 10
counts_sub <- count_sm[keep, , drop = FALSE]

# Ensure integers
counts_sub[] <- round(as.matrix(counts_sub))

# Build DESeq2 object (sex-only)
dds_sm <- DESeqDataSetFromMatrix(
  countData = counts_sub,
  colData   = as.data.frame(info) %>% tibble::column_to_rownames("ID"),
  design    = ~ sex
)

# Fit (poscounts is good with many zeros)
dds_sm <- DESeq(dds_sm, sfType = "poscounts")

# Sex contrast (male vs female)
res_sex <- results(dds_sm, contrast = c("sex","male","female"))

# Optional: LFC shrinkage (if apeglm available)
if (requireNamespace("apeglm", quietly = TRUE)) {
  suppressPackageStartupMessages(library(apeglm))
  res_sex <- lfcShrink(dds_sm, contrast = c("sex","male","female"), type = "apeglm")
}

# Save results
dir.create("deseq2_results", showWarnings = FALSE)
write.csv(as.data.frame(res_sex),
          file = file.path("deseq2_results_hum", "DESeq2_sex_SM.csv"),
          row.names = TRUE)

# Quick summary
cat("Sex results summary (SM):\n")
print(summary(res_sex))

