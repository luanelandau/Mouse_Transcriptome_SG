#Creating tx2gene (transcript to gene) annotation files

setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/normalization_checks/salmon_reanalysis/")

#Importing gtf files across species:
library(rtracklayer)
library(readr)
library(dplyr)
library(stringr)
library(tibble)

##################### ENSEMBL GTFs ######################################################

#mouse
#this is not in the folder but it is the gtf file for GRCm39 version 115 in ensembl.
mouse_gtf <- read_tsv("~/Library/CloudStorage/Box-Box/ancestral_sg/salmon/gtf/Mus_musculus.GRCm39.115.gtf.gz", 
                         comment = "#",col_names = FALSE, show_col_types = FALSE) 

colnames(mouse_gtf) <- c("seqname", "source", "feature", "start", "end","score", "strand", "frame", "attribute")
tx_mouse <- mouse_gtf %>% filter(feature == "transcript")

tx2gene_mouse <- tx_mouse %>%
  transmute(
    transcript_id = str_match(attribute, 'transcript_id "([^"]+)"')[,2],
    gene_id       = str_match(attribute, 'gene_id "([^"]+)"')[,2],
    gene_name     = str_match(attribute, 'gene_name "([^"]+)"')[,2]
  ) %>%
  distinct()


#Clean for NAs
tx2gene_mouse$gene_name[is.na(tx2gene_mouse$gene_name)] <- tx2gene_mouse$gene_id[is.na(tx2gene_mouse$gene_name)]
tx2gene_mouse <- tx2gene_mouse %>%
  filter(!is.na(transcript_id), !is.na(gene_id))

head(tx2gene_mouse)

write_csv(tx2gene_mouse, "mouse_tx2gene.csv")
