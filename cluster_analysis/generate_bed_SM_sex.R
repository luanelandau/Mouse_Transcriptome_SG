# Load libraries
library(dplyr)
library(readr)
library(stringr)

setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/")

# Load DEGs (significant genes)
deg <- read_csv("deseq2/deseq2_results/DESeq2_sex_SM.csv")

colnames(deg)[1] <- "Geneid"

#filter by significant genes
deg_ss = subset(deg, deg$padj<0.05 & abs(deg$log2FoldChange)>1,)

# Load the coords file
coords <- read_csv("miscelaneous_sheets/gene_expression_matrix_C57_CD1_TPMs.csv") %>%
  select(Geneid, Chr, Start, End)

# Clean the Chr, Start, End to use only the first value
coords_clean <- coords %>%
  mutate(
    Chr = str_split(Chr, ";", simplify = TRUE)[,1],
    Start = str_split(Start, ";", simplify = TRUE)[,1],
    End = str_split(End, ";", simplify = TRUE)[,1]
  )

deg_coords <- inner_join(deg_ss, coords_clean, by = "Geneid")

bed <- deg_coords %>%
  select(Chr, Start, End, Geneid) %>%
  arrange(Chr, as.numeric(Start))  # numeric sort is important

write.table(bed, "cluster_analysis/SM_sex_DEGs.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


#Here I need to generate a bed file for all genes in DEG too 
deg_coords_all <- inner_join(deg, coords_clean, by = "Geneid")

bed_all <- deg_coords_all %>%
  select(Chr, Start, End, Geneid) %>%
  arrange(Chr, as.numeric(Start))  # numeric sort is important

write.table(bed_all, "cluster_analysis/DEGs.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

