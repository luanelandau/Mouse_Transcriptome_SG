# Human submandibular gland (SM): male versus female DESeq2

setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/")

library(DESeq2)
library(tidyverse)

# Load the human raw-count matrix.
dat <- read.csv("miscelaneous_sheets/gene_expression_matrix_human_RAW_COUNTS.csv")

# Remove featureCounts annotation columns, leaving Geneid and sample counts.
drops <- c("Chr", "Start", "End", "Strand", "Length")
dat <- dat[, !(names(dat) %in% drops)]

# Make a count matrix with gene names as row names.
count_data <- dat[, -1, drop = FALSE]
rownames(count_data) <- dat[, 1]

# Keep only the six human submandibular gland samples.
sm_cols <- grep("^adult_SM_", colnames(count_data), value = TRUE)
count_sm <- count_data[, sm_cols, drop = FALSE]

# Describe the sex of each SM sample.
# Female is the reference level, so positive fold changes mean higher in males.
info <- tibble(
  ID = colnames(count_sm)
) %>%
  mutate(
    sex = case_when(
      ID %in% c("adult_SM_02", "adult_SM_04", "adult_SM_05") ~ "male",
      TRUE ~ "female"
    ),
    sex = factor(sex, levels = c("female", "male"))
  )

# Confirm that metadata rows and count-matrix columns are in the same order.
stopifnot(all(info$ID == colnames(count_sm)))

# Keep genes with at least 10 total reads across the six SM samples.
keep <- rowSums(count_sm) >= 10
counts_sub <- count_sm[keep, , drop = FALSE]

# DESeq2 requires integer raw counts.
counts_sub[] <- round(as.matrix(counts_sub))

# Move sample IDs to row names because DESeq2 matches these to count columns.
info_deseq <- as.data.frame(info) %>%
  column_to_rownames("ID")

# Construct and fit a sex-only DESeq2 model.
dds_sm <- DESeqDataSetFromMatrix(
  countData = counts_sub,
  colData = info_deseq,
  design = ~ sex
)

dds_sm <- DESeq(dds_sm, sfType = "poscounts")

# Compare males with females.
res_sex <- results(
  dds_sm,
  contrast = c("sex", "male", "female")
)

# Save the complete DESeq2 results table.
write.csv(
  as.data.frame(res_sex),
  file = "deseq2_hum/DESeq2_sex_SM.csv",
  row.names = TRUE
)

summary(res_sex)
