setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/")
library("DESeq2")
library("ggplot2")

#input the csv file
dat <- read.csv("miscelaneous_sheets/gene_expression_matrix_C57_CD1.csv")

#dropping the columns we dont want just for the purpose of this analysis
drops <- c("Chr","Start","End","Strand","Length")
head(drops)
dat <- dat[ , !(names(dat) %in% drops)] #excludes columns with the names defined on drops
colnames(dat)

#Note: B is for C57, whereas A is for CD1. 

#################Muscle contamination
#After analyzing the dataset, we found likely muscle contamination on some individuals, mainly based on the gene expression. Here are the
#individuals that should be dropped from the dataset: 

drops=c("Mous.1A.PAR.Mal.L", "Mous.2A.PAR.Mal.L", "Mous.6A.PAR.Fem.L", "Mous.4B.PAR.Fem.L", "Mous.5B.PAR.Fem.L", "Mous.6B.PAR.Fem.L")
dat <- dat[ , !(names(dat) %in% drops)] #excludes columns with the names defined on drops
colnames(dat)

##############################################################

count_data <- dat[, -1] #make sure this is relevant
rownames(count_data) <- dat[, 1] # Set the first column (geneID) as row names

library(dplyr)
library(stringr)

# Example data frame
info <- data.frame(
  ID = colnames(count_data),
  sex = NA,
  strain = NA,
  gland = NA,
  stringsAsFactors = FALSE
)

# Transform the dataframe
info <- info %>%
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
  )

# View result
print(info)

## =========================
## PCA (recommended design)
## =========================

colors_gland=c("SM"="#56B4E9",
               "SL"="#CC79A7", 
               "PANC"="#D55E00",
               "LIV"="#F5C710",
               "PAR"="#009E73")

# Rebuild DESeq object controlling for gland (and optionally sex)
dds_strain2 <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData   = info,
  design    = ~ gland + strain + sex   # add + sex if you want: ~ gland + sex + strain
)

# Filter low counts
dds_strain2 <- dds_strain2[rowSums(counts(dds_strain2)) >= 10, ]

#total number of genes kept
nrow(dds_strain2)

# Run DESeq
dds_strain2 <- DESeq(dds_strain2)

# VST and PCA
vsdata2 <- vst(dds_strain2, blind = FALSE)

pcaData2 <- plotPCA(vsdata2, intgroup = c("gland", "strain"), returnData = TRUE)
percentVar2 <- round(100 * attr(pcaData2, "percentVar"), 2)

# Plot
pca_strains2 <- ggplot(pcaData2, aes(PC1, PC2, color = gland, shape = strain)) +
  geom_point(size = 3) +
  theme_minimal() +
  scale_color_manual(values = colors_gland) +
  labs(
    title = "PCA all glands (design: gland + strain)",
    x = paste0("PC1: ", percentVar2[1], "% variance"),
    y = paste0("PC2: ", percentVar2[2], "% variance")
  )

pca_strains2

ggsave("deseq2/figures/PCA_allglands_new.svg", plot = pca_strains2)
ggsave("deseq2/figures/PCA_allglands_new.png", plot = pca_strains2, dpi=300, height = 5, width = 5 )

##############################################################
################### PCA SALIVARY GLANDS ONLY ##################
##############################################################

library(DESeq2)
library(dplyr)
library(ggplot2)

## Make sure metadata rows correspond exactly to count-data columns
rownames(info) <- info$ID

## Convert metadata variables to factors
info$sex <- factor(info$sex, levels = c("female", "male"))
info$strain <- factor(info$strain, levels = c("CD1", "C57"))
info$gland <- factor(info$gland)

## Confirm that the metadata and count matrix are in the same order
all(rownames(info) == colnames(count_data))
# This should return TRUE


##############################################################
## Select PAR, SL, and SM samples
##############################################################

salivary_samples <- info$gland %in% c("PAR", "SL", "SM")
count_data_salivary <- count_data[, salivary_samples]
info_salivary <- info[salivary_samples, , drop = FALSE]

## Remove unused LIV and PANC factor levels
info_salivary$gland <- droplevels(info_salivary$gland)
info_salivary$strain <- droplevels(info_salivary$strain)
info_salivary$sex <- droplevels(info_salivary$sex)

##############################################################
## Build a NEW DESeq2 object using only salivary glands
##############################################################

dds_salivary <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(count_data_salivary)),
  colData   = info_salivary,
  design    = ~ gland + strain + sex)

## Remove genes with fewer than 10 total counts across the
## salivary-gland samples
dds_salivary <- dds_salivary[rowSums(counts(dds_salivary)) >= 10,]

nrow(dds_salivary)

## Run DESeq2 using only the salivary-gland dataset
dds_salivary <- DESeq(dds_salivary)

## Variance-stabilizing transformation
vsdata_salivary <- vst(dds_salivary,blind = FALSE)

##############################################################
## PCA SGs: color by strain and shape by sex
##############################################################

pcaData_salivary <- plotPCA(
  vsdata_salivary,
  intgroup = c("gland", "strain", "sex"),
  returnData = TRUE
)

percentVar_salivary <- round(
  100 * attr(pcaData_salivary, "percentVar"),
  2
)

colors_strain <- c(
  "CD1" = "#f0debe",
  "C57" = "#700909"
)

pca_salivary_strain_sex <- ggplot(
  pcaData_salivary,
  aes(
    x = PC1,
    y = PC2,
    color = strain,
    shape = sex
  )
) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = colors_strain) +
  theme_minimal() +
  labs(
    title = "PCA of Salivary Gland Gene Expression",
    subtitle = "Color indicates strain; shape indicates sex",
    x = paste0(
      "PC1: ",
      percentVar_salivary[1],
      "% variance"
    ),
    y = paste0(
      "PC2: ",
      percentVar_salivary[2],
      "% variance"
    ),
    color = "Strain",
    shape = "Sex"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

pca_salivary_strain_sex

ggsave(
  filename = "deseq2/figures/PCA_salivary_strain_sex.svg",
  plot = pca_salivary_strain_sex,
  width = 6,
  height = 5
)

ggsave(
  filename = "deseq2/figures/PCA_salivary_strain_sex.png",
  plot = pca_salivary_strain_sex,
  width = 6,
  height = 5,
  dpi = 300
)