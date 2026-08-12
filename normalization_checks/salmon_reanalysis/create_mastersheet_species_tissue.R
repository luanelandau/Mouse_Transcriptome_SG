library(readr)
library(dplyr)
library(purrr)

setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/normalization_checks/salmon_reanalysis/")

gene_dir <- "genelevel"
samples <- read_csv("sample_info.csv", show_col_types = FALSE)

dir.create("tpm", showWarnings = FALSE)

combos <- samples %>%
  distinct(species, tissue)

for (i in 1:nrow(combos)) {
  
  sp  <- combos$species[i]
  tis <- combos$tissue[i]
  
  cat("Processing:", sp, tis, "\n")
  
  sub_samples <- samples %>%
    filter(species == sp, tissue == tis)
  
  files <- paste0(gene_dir, "/gene_level_", sub_samples$sample, ".csv")
  
  keep <- file.exists(files)
  files <- files[keep]
  sub_samples <- sub_samples[keep, ]
  
  if (length(files) == 0) next
  
  df_list <- map2(files, sub_samples$sample, function(f, s) {
    
    df <- read_csv(f, show_col_types = FALSE) %>%
      select(gene_id, gene_name, TPM)
    
    colnames(df)[colnames(df) == "TPM"] <- s
    
    df
  })
  
  merged <- purrr::reduce(df_list, full_join, by = c("gene_id", "gene_name"))
  
  sample_cols <- sub_samples$sample
  
  merged <- merged %>%
    mutate(
      mean_tpm = rowMeans(across(all_of(sample_cols)), na.rm = TRUE),
      sd_tpm = apply(across(all_of(sample_cols)), 1, sd, na.rm = TRUE)
    ) %>%
    select(gene_id, gene_name, mean_tpm, sd_tpm, everything())
  
  outname <- paste0("tpm/", sp, "_", tis, "_TPM_mastersheet.csv")
  write_csv(merged, outname)
}
