library(tximport)
library(readr)
library(dplyr)
library(tibble)

#BiocManager::install("tximport")


setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/normalization_checks/salmon_reanalysis/")

#mouse only
samples <- read_tsv(
  "quant_paths.txt",
  show_col_types = FALSE
)
dir.create("genelevel", showWarnings = FALSE)

# output folder
outdir <- "genelevel"

for (i in 1:nrow(samples)) {
  
  species <- samples$species[i]
  sample_name <- samples$sample[i]
  quant_file <- samples$path[i]
  
  cat("Processing:", species, sample_name, "\n")
  
  tx2gene = read_csv("mouse_tx2gene.csv", show_col_types = FALSE)
  ####
  
  tx2gene2 <- tx2gene %>%
    select(transcript_id, gene_id)
  
  tx2gene2 <- tx2gene %>%
    mutate(
      transcript_id = sub("\\.[0-9]+$", "", transcript_id)
    ) %>%
    select(transcript_id, gene_id) %>%
    distinct()
  ###############
  gene_annot <- tx2gene %>%
    select(gene_id, gene_name) %>%
    distinct()
  
  # ---- run tximport ----
  txi <- tximport(
    quant_file,
    type = "salmon",
    tx2gene = tx2gene2,
    ignoreTxVersion = TRUE
  )
  
  # ---- build dataframe ----
  df <- data.frame(
    gene_id = rownames(txi$abundance),
    TPM     = txi$abundance[,1],
    counts  = txi$counts[,1],
    length  = txi$length[,1]
  )
  
  df <- df %>%
    left_join(gene_annot, by = "gene_id") %>%
    select(gene_id, gene_name, TPM, counts, length)
  
  # ---- save ----
  write_csv(
    df,
    file.path(outdir, paste0("gene_level_", sample_name, ".csv"))
  )
}