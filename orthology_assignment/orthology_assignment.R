
setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/orthology_assignment/")

#File with ortholog names
orthologs <- read.table("Human_to_Mouse_Orthologs_HMD_HumanPhenotype.rpt", 
                                 sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")

colnames(orthologs) <- c("human_gene", "human_id", "mouse_gene", "mouse_id", "phenotype", "extra")

library(dplyr)

# Count how many mouse genes each human gene maps to
human_to_mouse_counts <- orthologs %>%
  group_by(human_gene) %>%
  summarize(n_mouse = n_distinct(mouse_gene), .groups = "drop")

# Count how many human genes each mouse gene maps to
mouse_to_human_counts <- orthologs %>%
  group_by(mouse_gene) %>%
  summarize(n_human = n_distinct(human_gene), .groups = "drop")

# Merge counts back to the main dataframe
orthologs_classified <- orthologs %>%
  left_join(human_to_mouse_counts, by = "human_gene") %>%
  left_join(mouse_to_human_counts, by = "mouse_gene") %>%
  mutate(ortholog_type = case_when(
    n_mouse == 1 & n_human == 1 ~ "one-to-one",
    n_mouse > 1 & n_human == 1 ~ "one-to-many (human-to-mouse)",
    n_mouse == 1 & n_human > 1 ~ "one-to-many (mouse-to-human)",
    n_mouse > 1 & n_human > 1 ~ "many-to-many"
  ))

# View result
head(orthologs_classified)

# Save it if needed
write.csv(orthologs_classified, "orthologs_with_classification.csv", row.names = FALSE)
