#!/usr/bin/env Rscript

# GO Biological Process enrichment of all significant DESeq2 genes in each
# mouse salivary gland, without separating genes by orthology category.
#
# Significant foreground:
#   abs(log2FoldChange) > 1 and padj < 0.05
#
# Gland-specific background:
#   all genes with a non-missing DESeq2 adjusted p-value. These are the genes
#   that had an opportunity to satisfy the foreground selection criteria.
#
# The enrichGO pvalueCutoff and qvalueCutoff below apply to GO terms. They are
# separate from the DESeq2 padj threshold used to select significant genes.

library(clusterProfiler)
library(org.Mm.eg.db)

input_folder <- "~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/deseq2/deseq2_results"
output_folder <- "~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/go_enrichment/all_significant_genes"

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# PAROTID GLAND (PAR)
# -----------------------------------------------------------------------------

PAR <- read.csv(
  file.path(input_folder, "DESeq2_sex_PAR.csv"),
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
PAR$Geneid <- rownames(PAR)

# All PAR genes eligible to meet the padj threshold.
PAR_background <- PAR[!is.na(PAR$padj), ]

# All significant PAR genes, regardless of orthology category.
PAR_significant <- PAR_background[
  abs(PAR_background$log2FoldChange) > 1 &
    PAR_background$padj < 0.05,
]

PAR_background_ids <- bitr(
  unique(PAR_background$Geneid),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

PAR_significant_ids <- bitr(
  unique(PAR_significant$Geneid),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

PAR_GO_BP <- enrichGO(
  gene = unique(PAR_significant_ids$ENTREZID),
  universe = unique(PAR_background_ids$ENTREZID),
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

#write.csv(
#  PAR_background,
#  file.path(output_folder, "PAR_sex_background_nonmissing_padj.csv"),
#  row.names = FALSE
#)
#write.csv(
#  PAR_significant,
#  file.path(output_folder, "PAR_sex_all_significant_absLFC1_padj0.05_genes.csv"),
#  row.names = FALSE
#)
#write.csv(
#  PAR_background_ids,
#  file.path(output_folder, "PAR_sex_background_nonmissing_padj_mapped_ids.csv"),
#  row.names = FALSE
#)
#write.csv(
#  PAR_significant_ids,
#  file.path(output_folder, "PAR_sex_all_significant_absLFC1_padj0.05_mapped_ids.csv"),
#  row.names = FALSE
#)
write.csv(
  as.data.frame(PAR_GO_BP),
  file.path(output_folder, "PAR_sex_all_significant_absLFC1_padj0.05_GO_BP.csv"),
  row.names = FALSE
)

cat(
  "PAR background:", nrow(PAR_background), "genes;",
  nrow(PAR_background_ids), "mapped. Significant:",
  nrow(PAR_significant), "genes;", nrow(PAR_significant_ids), "mapped. GO BP terms:",
  nrow(as.data.frame(PAR_GO_BP)), "\n"
)


# -----------------------------------------------------------------------------
# SUBLINGUAL GLAND (SL)
# -----------------------------------------------------------------------------

SL <- read.csv(
  file.path(input_folder, "DESeq2_sex_SL.csv"),
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
SL$Geneid <- rownames(SL)

# All SL genes eligible to meet the padj threshold.
SL_background <- SL[!is.na(SL$padj), ]

# All significant SL genes, regardless of orthology category.
SL_significant <- SL_background[
  abs(SL_background$log2FoldChange) > 1 &
    SL_background$padj < 0.05,
]

SL_background_ids <- bitr(
  unique(SL_background$Geneid),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

SL_significant_ids <- bitr(
  unique(SL_significant$Geneid),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

SL_GO_BP <- enrichGO(
  gene = unique(SL_significant_ids$ENTREZID),
  universe = unique(SL_background_ids$ENTREZID),
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

#write.csv(
#  SL_background,
#  file.path(output_folder, "SL_sex_background_nonmissing_padj.csv"),
#  row.names = FALSE
#)
#write.csv(
#  SL_significant,
#  file.path(output_folder, "SL_sex_all_significant_absLFC1_padj0.05_genes.csv"),
#  row.names = FALSE
#)
#write.csv(
#  SL_background_ids,
#  file.path(output_folder, "SL_sex_background_nonmissing_padj_mapped_ids.csv"),
#  row.names = FALSE
#)
#write.csv(
#  SL_significant_ids,
#  file.path(output_folder, "SL_sex_all_significant_absLFC1_padj0.05_mapped_ids.csv"),
#  row.names = FALSE
#)
write.csv(
  as.data.frame(SL_GO_BP),
  file.path(output_folder, "SL_sex_all_significant_absLFC1_padj0.05_GO_BP.csv"),
  row.names = FALSE
)

cat(
  "SL background:", nrow(SL_background), "genes;",
  nrow(SL_background_ids), "mapped. Significant:",
  nrow(SL_significant), "genes;", nrow(SL_significant_ids), "mapped. GO BP terms:",
  nrow(as.data.frame(SL_GO_BP)), "\n"
)


# -----------------------------------------------------------------------------
# SUBMANDIBULAR GLAND (SM)
# -----------------------------------------------------------------------------

SM <- read.csv(
  file.path(input_folder, "DESeq2_sex_SM.csv"),
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
SM$Geneid <- rownames(SM)

# All SM genes eligible to meet the padj threshold.
SM_background <- SM[!is.na(SM$padj), ]

# All significant SM genes, regardless of orthology category.
SM_significant <- SM_background[
  abs(SM_background$log2FoldChange) > 1 &
    SM_background$padj < 0.05,
]

SM_background_ids <- bitr(
  unique(SM_background$Geneid),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

SM_significant_ids <- bitr(
  unique(SM_significant$Geneid),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

SM_GO_BP <- enrichGO(
  gene = unique(SM_significant_ids$ENTREZID),
  universe = unique(SM_background_ids$ENTREZID),
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

#write.csv(
#  SM_background,
#  file.path(output_folder, "SM_sex_background_nonmissing_padj.csv"),
#  row.names = FALSE
#)
#write.csv(
#  SM_significant,
#  file.path(output_folder, "SM_sex_all_significant_absLFC1_padj0.05_genes.csv"),
#  row.names = FALSE
#)
#write.csv(
#  SM_background_ids,
#  file.path(output_folder, "SM_sex_background_nonmissing_padj_mapped_ids.csv"),
#  row.names = FALSE
#)
#write.csv(
#  SM_significant_ids,
#  file.path(output_folder, "SM_sex_all_significant_absLFC1_padj0.05_mapped_ids.csv"),
#  row.names = FALSE
#)
write.csv(
  as.data.frame(SM_GO_BP),
  file.path(output_folder, "SM_sex_all_significant_absLFC1_padj0.05_GO_BP.csv"),
  row.names = FALSE
)

cat(
  "SM background:", nrow(SM_background), "genes;",
  nrow(SM_background_ids), "mapped. Significant:",
  nrow(SM_significant), "genes;", nrow(SM_significant_ids), "mapped. GO BP terms:",
  nrow(as.data.frame(SM_GO_BP)), "\n"
)


# Summary of the three enrichment analyses.
analysis_summary <- data.frame(
  gland = c("PAR", "SL", "SM"),
  background_genes = c(
    nrow(PAR_background),
    nrow(SL_background),
    nrow(SM_background)
  ),
  background_mapped_ids = c(
    nrow(PAR_background_ids),
    nrow(SL_background_ids),
    nrow(SM_background_ids)
  ),
  significant_genes = c(
    nrow(PAR_significant),
    nrow(SL_significant),
    nrow(SM_significant)
  ),
  significant_mapped_ids = c(
    nrow(PAR_significant_ids),
    nrow(SL_significant_ids),
    nrow(SM_significant_ids)
  ),
  significant_GO_BP_terms = c(
    nrow(as.data.frame(PAR_GO_BP)),
    nrow(as.data.frame(SL_GO_BP)),
    nrow(as.data.frame(SM_GO_BP))
  )
)

#write.csv(
#  analysis_summary,
#  file.path(output_folder, "all_glands_sex_absLFC1_padj0.05_GO_BP_summary.csv"),
#  row.names = FALSE
#)
#