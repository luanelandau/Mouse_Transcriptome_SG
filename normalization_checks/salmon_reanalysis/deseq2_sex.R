############################################################
# DESeq2 per tissue: Male vs Female
# Mouse salivary gland and other tissues
############################################################

library(DESeq2)
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(stringr)

############################################################
# Set working directory
############################################################

setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/normalization_checks/salmon_reanalysis/")

############################################################
# Create output folders
############################################################

dir.create("deseq2_sex", showWarnings = FALSE)
dir.create("figures_sex", showWarnings = FALSE)

############################################################
# Read metadata
############################################################

##IMPORTANT: I am not taking into consideration the strains here, as I want to 
#treat this as mainly sex differences across strains. I could, and therefore the
#variation would be CORRECTED for strain. But since i didnt do this in the initial
#pipeline, I do not want to complicate things. 

# The metadata file should contain columns similar to:
#
# sample_name          ind       tissue   sex
# Mous-1A-LIV-Mal      Mous-1A   LIV      Male
# Mous-1A-LIV-Fem      Mous-1A   LIV      Female
#
# Change the file name below if necessary.

info <- read_tsv("sample_info_complete.txt", show_col_types = FALSE)

# Rename sample_name to ID because the script uses ID below.
# Remove this line if your metadata already has a column called ID.

info <- info %>%
  dplyr::rename(ID = sample_name)

# Standardize sex labels to Male and Female.
# This also recognizes values such as male, female, Mal and Fem.

info <- info %>%
  mutate(
    sex = case_when(
      str_to_lower(sex) %in% c("male", "mal", "m") ~ "Male",
      str_to_lower(sex) %in% c("female", "fem", "f") ~ "Female",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(sex))

# Tissues present in the metadata

tissues <- sort(unique(info$tissue))

print(table(info$tissue, info$sex))

############################################################
# Empty table to collect DEGs from all tissues
############################################################

all_deg_list <- list()

############################################################
# Loop through tissues
############################################################

for (this_tissue in tissues) {
  
  message("========================================")
  message("Processing tissue: ", this_tissue)
  message("========================================")
  
  ##########################################################
  # Subset metadata for this tissue
  ##########################################################
  
  info_tissue <- info %>%
    filter(tissue == this_tissue)
  
  # Check whether both sexes are present
  
  group_counts <- table(info_tissue$sex)
  print(group_counts)
  
  if (!all(c("Female", "Male") %in% names(group_counts))) {
    message(
      "Skipping ", this_tissue,
      ": missing Male or Female group."
    )
    next
  }
  
  if (any(group_counts[c("Female", "Male")] < 2)) {
    message(
      "Skipping ", this_tissue,
      ": fewer than 2 samples in one sex."
    )
    next
  }
  
  ##########################################################
  # Read gene-level files for this tissue
  ##########################################################
  
  count_list <- list()
  gene_annot <- NULL
  
  for (i in seq_len(nrow(info_tissue))) {
    
    sample_id <- info_tissue$ID[i]
    
    file_path <- file.path(
      "genelevel",
      paste0("gene_level_", sample_id, ".csv")
    )
    
    if (!file.exists(file_path)) {
      warning("File not found: ", file_path)
      next
    }
    
    message("Reading: ", file_path)
    
    df <- read_csv(file_path, show_col_types = FALSE)
    
    # Keep gene annotation from the first file
    
    if (is.null(gene_annot)) {
      
      gene_annot <- df %>%
        select(gene_id, gene_name) %>%
        distinct(gene_id, .keep_all = TRUE)
    }
    
    # Extract counts
    
    sample_counts <- df %>%
      select(gene_id, counts)
    
    colnames(sample_counts)[
      colnames(sample_counts) == "counts"
    ] <- sample_id
    
    count_list[[sample_id]] <- sample_counts
  }
  
  ##########################################################
  # Skip if files were missing
  ##########################################################
  
  if (length(count_list) < 4) {
    message(
      "Skipping ", this_tissue,
      ": fewer than 4 count files found."
    )
    next
  }
  
  ##########################################################
  # Merge counts into one matrix
  ##########################################################
  
  counts_df <- Reduce(
    function(x, y) full_join(x, y, by = "gene_id"),
    count_list
  )
  
  # Replace missing counts with zero
  
  counts_df[is.na(counts_df)] <- 0
  
  # Convert to matrix
  
  counts_mat <- counts_df %>%
    column_to_rownames("gene_id") %>%
    as.matrix()
  
  # DESeq2 expects integer counts
  
  counts_mat <- round(counts_mat)
  
  storage.mode(counts_mat) <- "integer"
  
  ##########################################################
  # Prepare metadata
  ##########################################################
  
  coldata <- info_tissue %>%
    filter(ID %in% colnames(counts_mat)) %>%
    as.data.frame()
  
  rownames(coldata) <- coldata$ID
  
  # Make metadata order match count-matrix columns
  
  coldata <- coldata[colnames(counts_mat), , drop = FALSE]
  
  # Check that metadata and counts match
  
  if (!all(rownames(coldata) == colnames(counts_mat))) {
    stop(
      "Metadata and count matrix do not match for tissue: ",
      this_tissue
    )
  }
  
  ##########################################################
  # Recheck sex counts after removing missing files
  ##########################################################
  
  final_group_counts <- table(coldata$sex)
  print(final_group_counts)
  
  if (!all(c("Female", "Male") %in% names(final_group_counts))) {
    message(
      "Skipping ", this_tissue,
      ": missing Male or Female after checking count files."
    )
    next
  }
  
  if (any(final_group_counts[c("Female", "Male")] < 2)) {
    message(
      "Skipping ", this_tissue,
      ": fewer than 2 samples in one sex after checking count files."
    )
    next
  }
  
  ##########################################################
  # Set factor levels
  ##########################################################
  
  # Female is the reference level.
  # Therefore:
  #
  # Positive log2FoldChange = higher expression in Male
  # Negative log2FoldChange = higher expression in Female
  
  coldata$sex <- factor(
    coldata$sex,
    levels = c("Female", "Male")
  )
  
  ##########################################################
  # Run DESeq2
  ##########################################################
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData = coldata,
    design = ~ sex
  )
  
  # Remove very low-count genes
  
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  
  dds <- DESeq(dds)
  
  res <- results(
    dds,
    contrast = c("sex", "Male", "Female"),
    alpha = 0.05
  )
  
  res_df <- as.data.frame(res)
  
  ##########################################################
  # Add gene_id and gene_name
  ##########################################################
  
  res_df$gene_id <- rownames(res_df)
  
  res_df <- res_df %>%
    left_join(gene_annot, by = "gene_id") %>%
    select(gene_id, gene_name, everything()) %>%
    arrange(padj)
  
  ##########################################################
  # Collect significant DEGs for final combined table
  ##########################################################
  
  deg_this_tissue <- res_df %>%
    filter(
      !is.na(padj),
      padj < 0.05,
      !is.na(log2FoldChange),
      abs(log2FoldChange) >= 1
    ) %>%
    mutate(
      higher_in = case_when(
        log2FoldChange >= 1 ~ "Male",
        log2FoldChange <= -1 ~ "Female"
      )
    ) %>%
    transmute(
      ens_id = gene_id,
      gene = gene_name,
      tissue = this_tissue,
      higher_in = higher_in,
      log2fold = log2FoldChange,
      padj = padj,
      basemean = baseMean
    )
  
  all_deg_list[[this_tissue]] <- deg_this_tissue
  
  ##########################################################
  # Save DESeq2 table
  ##########################################################
  
  out_csv <- file.path(
    "deseq2_sex",
    paste0(
      "DESeq2_",
      this_tissue,
      "_Male_vs_Female.csv"
    )
  )
  
  write.csv(
    res_df,
    out_csv,
    row.names = FALSE
  )
  
  message("Saved DESeq2 table: ", out_csv)
  
  ##########################################################
  # Volcano plot
  ##########################################################
  
  volcano_df <- res_df %>%
    filter(
      !is.na(log2FoldChange),
      !is.na(padj)
    ) %>%
    mutate(
      padj_plot = ifelse(
        padj == 0,
        1e-300,
        padj
      ),
      neglog10_padj = -log10(padj_plot),
      significance = case_when(
        padj < 0.05 & log2FoldChange >= 1 ~
          "Higher in Male",
        
        padj < 0.05 & log2FoldChange <= -1 ~
          "Higher in Female",
        
        TRUE ~
          "Not significant"
      )
    )
  
  # Label the 15 most significant DEGs
  
  label_df <- volcano_df %>%
    filter(
      significance != "Not significant",
      !is.na(gene_name),
      gene_name != ""
    ) %>%
    arrange(padj) %>%
    slice_head(n = 15)
  
  p <- ggplot(
    volcano_df,
    aes(
      x = log2FoldChange,
      y = neglog10_padj
    )
  ) +
    geom_point(
      aes(color = significance),
      alpha = 0.7,
      size = 2
    ) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    geom_text_repel(
      data = label_df,
      aes(label = gene_name),
      size = 3,
      max.overlaps = 20
    ) +
    scale_color_manual(
      values = c(
        "Higher in Male" = "#043a5c",
        "Higher in Female" = "#d87171",
        "Not significant" = "grey70"
      )
    ) +
    labs(
      title = paste0(
        "Volcano plot: ",
        this_tissue,
        " Male vs Female"
      ),
      subtitle = paste0(
        "Male: n = ",
        final_group_counts["Male"],
        "; Female: n = ",
        final_group_counts["Female"]
      ),
      x = "log2 fold change (Male / Female)",
      y = "-log10 adjusted p-value",
      color = NULL
    ) +
    theme_classic()
  
  out_png <- file.path(
    "figures_sex",
    paste0(
      "volcano_",
      this_tissue,
      "_Male_vs_Female.png"
    )
  )
  
  ggsave(
    out_png,
    plot = p,
    width = 7,
    height = 6,
    dpi = 300
  )
  
  message("Saved volcano plot: ", out_png)
}
