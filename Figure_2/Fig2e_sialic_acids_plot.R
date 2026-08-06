library(ggplot2)
library(dplyr)
library(tidyr)

setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/")

genes <- c("St3gal1", "St3gal4", "St3gal5", "St6galnac2", "St6galnac1")

sm_sheet <- read.csv("miscelaneous_sheets/mouse_expression/SM_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv")

sial_acid_long <- sm_sheet %>%
  filter(Geneid %in% genes) %>%
  select(
    Geneid,
    starts_with("Mous.")
  ) %>%
  pivot_longer(
    cols = -Geneid,
    names_to = "Sample",
    values_to = "TPM"
  ) %>%
  mutate(
    Sex = case_when(
      grepl("Mal", Sample) ~ "Male",
      grepl("Fem", Sample) ~ "Female",
      TRUE ~ NA_character_
    ),
    Geneid = factor(Geneid, levels = genes),
    Sex = factor(Sex, levels = c("Male", "Female"))
  ) %>%
  filter(!is.na(Sex))

sial_acid_long <- sm_sheet %>%
  filter(Geneid %in% genes) %>%
  select(
    Geneid,
    starts_with("Mous.")
  ) %>%
  pivot_longer(
    cols = -Geneid,
    names_to = "Sample",
    values_to = "TPM"
  ) %>%
  mutate(
    Sex = case_when(
      grepl("Mal", Sample) ~ "Male",
      grepl("Fem", Sample) ~ "Female",
      TRUE ~ NA_character_
    ),
    Geneid = factor(Geneid, levels = genes),
    Sex = factor(Sex, levels = c("Male", "Female"))
  ) %>%
  filter(!is.na(Sex))

#this is figure S5
boxplot_sialic=ggplot(sial_acid_long, aes(x = Geneid, y = TPM, fill = Sex)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.shape = NA, alpha = 0.4
  ) +
  geom_jitter(
    aes(color = Sex),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
    size = 2,
    alpha = 1
  ) +
  theme_classic() +
  #labs(
  #x = "Gene"#,
  #y = "TPM"
  #) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = c("Male" = "#5d7b9f", "Female" = "#d391a0")) +
  scale_color_manual(values = c("Male" = "#5d7b9f", "Female" = "#d391a0")) #+
#coord_cartesian(ylim = c(0, 200))

boxplot_sialic

#This is for figure 2e
boxplot_sialic=ggplot(sial_acid_long, aes(x = Geneid, y = TPM, fill = Sex)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.shape = NA, alpha = 0.4
  ) +
  geom_jitter(
    aes(color = Sex),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
    size = 2,
    alpha = 1
  ) +
  theme_classic() +
  #labs(
  #x = "Gene"#,
  #y = "TPM"
  #) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = c("Male" = "#5d7b9f", "Female" = "#d391a0")) +
  scale_color_manual(values = c("Male" = "#5d7b9f", "Female" = "#d391a0")) + coord_cartesian(ylim = c(0, 200))

boxplot_sialic

#The boxplot does not contain the outlier in St3gal1, it was manually 
#added to the figure as a triangle to simbolize the outlier.
ggsave("Figure_2/boxplot_sialic_acids.png", boxplot_sialic)
