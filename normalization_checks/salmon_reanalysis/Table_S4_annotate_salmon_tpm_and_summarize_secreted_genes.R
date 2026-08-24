# Annotate Salmon TPM mastersheets and summarize secreted-gene expression.
# Run this script from the salmon_reanalysis directory with:
# Rscript annotate_salmon_tpm_and_summarize_secreted_genes.R

setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/normalization_checks/salmon_reanalysis/")

library(readr)
library(dplyr)

dir.create("salmon_annotated_secretion_orthology", showWarnings = FALSE)

annotation <- read_csv(
  "../../miscelaneous_sheets/mouse_expression/SM_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv",
  show_col_types = FALSE
) %>%
  select(Geneid, human_gene, ortholog_type, secreted, Reference, Function)

# PAR
par <- read_csv("tpm/mouse_PAR_TPM_mastersheet.csv", show_col_types = FALSE) %>%
  left_join(annotation, by = join_by(gene_name == Geneid)) %>%
  relocate(human_gene, ortholog_type, secreted, Reference, Function, .after = gene_name)

write_csv(
  par,
  "salmon_annotated_secretion_orthology/mouse_PAR_TPM_mastersheet_annotated_for_secretion_and_orthology.csv"
)

par_secreted <- par %>%
  filter(secreted %in% c("yes", "manual_annotation")) %>%
  arrange(desc(mean_tpm)) %>%
  mutate(
    secreted_rank = row_number(),
    percent_total_expression = 100 * mean_tpm / sum(par$mean_tpm, na.rm = TRUE),
    percent_secreted_expression = 100 * mean_tpm / sum(mean_tpm, na.rm = TRUE),
    cumulative_percent_total_expression = cumsum(percent_total_expression),
    cumulative_percent_secreted_expression = cumsum(percent_secreted_expression)
  )

par_top30 <- par_secreted %>% slice_head(n = 30)

#write_csv(
#  par_top30,
#  "salmon_annotated_secretion_orthology/PAR_top30_secreted_genes.csv"
#)

par_summary <- tibble(
  gland = "PAR",
  total_mean_tpm = sum(par$mean_tpm, na.rm = TRUE),
  secreted_mean_tpm = sum(par_secreted$mean_tpm, na.rm = TRUE),
  secreted_percent_of_total_expression = 100 * secreted_mean_tpm / total_mean_tpm,
  top30_secreted_mean_tpm = sum(par_top30$mean_tpm, na.rm = TRUE),
  top30_percent_of_total_expression = 100 * top30_secreted_mean_tpm / total_mean_tpm,
  top30_percent_of_secreted_expression = 100 * top30_secreted_mean_tpm / secreted_mean_tpm
)

# SM
sm <- read_csv("tpm/mouse_SM_TPM_mastersheet.csv", show_col_types = FALSE) %>%
  left_join(annotation, by = join_by(gene_name == Geneid)) %>%
  relocate(human_gene, ortholog_type, secreted, Reference, Function, .after = gene_name)

write_csv(
  sm,
  "salmon_annotated_secretion_orthology/mouse_SM_TPM_mastersheet_annotated_for_secretion_and_orthology.csv"
)

sm_secreted <- sm %>%
  filter(secreted %in% c("yes", "manual_annotation")) %>%
  arrange(desc(mean_tpm)) %>%
  mutate(
    secreted_rank = row_number(),
    percent_total_expression = 100 * mean_tpm / sum(sm$mean_tpm, na.rm = TRUE),
    percent_secreted_expression = 100 * mean_tpm / sum(mean_tpm, na.rm = TRUE),
    cumulative_percent_total_expression = cumsum(percent_total_expression),
    cumulative_percent_secreted_expression = cumsum(percent_secreted_expression)
  )


sm_top30 <- sm_secreted %>% slice_head(n = 30)

#write_csv(
#  sm_top30,
#  "salmon_annotated_secretion_orthology/SM_top30_secreted_genes.csv"
#)

sm_summary <- tibble(
  gland = "SM",
  total_mean_tpm = sum(sm$mean_tpm, na.rm = TRUE),
  secreted_mean_tpm = sum(sm_secreted$mean_tpm, na.rm = TRUE),
  secreted_percent_of_total_expression = 100 * secreted_mean_tpm / total_mean_tpm,
  top30_secreted_mean_tpm = sum(sm_top30$mean_tpm, na.rm = TRUE),
  top30_percent_of_total_expression = 100 * top30_secreted_mean_tpm / total_mean_tpm,
  top30_percent_of_secreted_expression = 100 * top30_secreted_mean_tpm / secreted_mean_tpm
)

# SL
sl <- read_csv("tpm/mouse_SL_TPM_mastersheet.csv", show_col_types = FALSE) %>%
  left_join(annotation, by = join_by(gene_name == Geneid)) %>%
  relocate(human_gene, ortholog_type, secreted, Reference, Function, .after = gene_name)

write_csv(
  sl,
  "salmon_annotated_secretion_orthology/mouse_SL_TPM_mastersheet_annotated_for_secretion_and_orthology.csv"
)

sl_secreted <- sl %>%
  filter(secreted %in% c("yes", "manual_annotation")) %>%
  arrange(desc(mean_tpm)) %>%
  mutate(
    secreted_rank = row_number(),
    percent_total_expression = 100 * mean_tpm / sum(sl$mean_tpm, na.rm = TRUE),
    percent_secreted_expression = 100 * mean_tpm / sum(mean_tpm, na.rm = TRUE),
    cumulative_percent_total_expression = cumsum(percent_total_expression),
    cumulative_percent_secreted_expression = cumsum(percent_secreted_expression)
  )

sl_top30 <- sl_secreted %>% slice_head(n = 30)

#write_csv(
#  sl_top30,
#  "salmon_annotated_secretion_orthology/SL_top30_secreted_genes.csv"
#)

sl_summary <- tibble(
  gland = "SL",
  total_mean_tpm = sum(sl$mean_tpm, na.rm = TRUE),
  secreted_mean_tpm = sum(sl_secreted$mean_tpm, na.rm = TRUE),
  secreted_percent_of_total_expression = 100 * secreted_mean_tpm / total_mean_tpm,
  top30_secreted_mean_tpm = sum(sl_top30$mean_tpm, na.rm = TRUE),
  top30_percent_of_total_expression = 100 * top30_secreted_mean_tpm / total_mean_tpm,
  top30_percent_of_secreted_expression = 100 * top30_secreted_mean_tpm / secreted_mean_tpm
)

# PANC
panc <- read_csv("tpm/mouse_PANC_TPM_mastersheet.csv", show_col_types = FALSE) %>%
  left_join(annotation, by = join_by(gene_name == Geneid)) %>%
  relocate(human_gene, ortholog_type, secreted, Reference, Function, .after = gene_name)

write_csv(
  panc,
  "salmon_annotated_secretion_orthology/mouse_PANC_TPM_mastersheet_annotated_for_secretion_and_orthology.csv"
)

panc_secreted <- panc %>%
  filter(secreted %in% c("yes", "manual_annotation")) %>%
  arrange(desc(mean_tpm)) %>%
  mutate(
    secreted_rank = row_number(),
    percent_total_expression = 100 * mean_tpm / sum(panc$mean_tpm, na.rm = TRUE),
    percent_secreted_expression = 100 * mean_tpm / sum(mean_tpm, na.rm = TRUE),
    cumulative_percent_total_expression = cumsum(percent_total_expression),
    cumulative_percent_secreted_expression = cumsum(percent_secreted_expression)
  )

panc_top30 <- panc_secreted %>% slice_head(n = 30)

#write_csv(
#  panc_top30,
#  "salmon_annotated_secretion_orthology/PANC_top30_secreted_genes.csv"
#)

panc_summary <- tibble(
  gland = "PANC",
  total_mean_tpm = sum(panc$mean_tpm, na.rm = TRUE),
  secreted_mean_tpm = sum(panc_secreted$mean_tpm, na.rm = TRUE),
  secreted_percent_of_total_expression = 100 * secreted_mean_tpm / total_mean_tpm,
  top30_secreted_mean_tpm = sum(panc_top30$mean_tpm, na.rm = TRUE),
  top30_percent_of_total_expression = 100 * top30_secreted_mean_tpm / total_mean_tpm,
  top30_percent_of_secreted_expression = 100 * top30_secreted_mean_tpm / secreted_mean_tpm
)

# LIV
liv <- read_csv("tpm/mouse_LIV_TPM_mastersheet.csv", show_col_types = FALSE) %>%
  left_join(annotation, by = join_by(gene_name == Geneid)) %>%
  relocate(human_gene, ortholog_type, secreted, Reference, Function, .after = gene_name)

write_csv(
  liv,
  "salmon_annotated_secretion_orthology/mouse_LIV_TPM_mastersheet_annotated_for_secretion_and_orthology.csv"
)

liv_secreted <- liv %>%
  filter(secreted %in% c("yes", "manual_annotation")) %>%
  arrange(desc(mean_tpm)) %>%
  mutate(
    secreted_rank = row_number(),
    percent_total_expression = 100 * mean_tpm / sum(liv$mean_tpm, na.rm = TRUE),
    percent_secreted_expression = 100 * mean_tpm / sum(mean_tpm, na.rm = TRUE),
    cumulative_percent_total_expression = cumsum(percent_total_expression),
    cumulative_percent_secreted_expression = cumsum(percent_secreted_expression)
  )

liv_top30 <- liv_secreted %>% slice_head(n = 30)

#write_csv(
#  liv_top30,
#  "salmon_annotated_secretion_orthology/LIV_top30_secreted_genes.csv"
#)

liv_summary <- tibble(
  gland = "LIV",
  total_mean_tpm = sum(liv$mean_tpm, na.rm = TRUE),
  secreted_mean_tpm = sum(liv_secreted$mean_tpm, na.rm = TRUE),
  secreted_percent_of_total_expression = 100 * secreted_mean_tpm / total_mean_tpm,
  top30_secreted_mean_tpm = sum(liv_top30$mean_tpm, na.rm = TRUE),
  top30_percent_of_total_expression = 100 * top30_secreted_mean_tpm / total_mean_tpm,
  top30_percent_of_secreted_expression = 100 * top30_secreted_mean_tpm / secreted_mean_tpm
)

secreted_expression_summary <- bind_rows(
  par_summary,
  sm_summary,
  sl_summary,
  panc_summary,
  liv_summary
)

write_csv(
  secreted_expression_summary,
  "salmon_annotated_secretion_orthology/secreted_expression_summary.csv"
)

print(secreted_expression_summary)




# HISAT secreted-gene expression summary

hisat_par <- read_csv(
  "../../miscelaneous_sheets/mouse_expression/PAR_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv",
  show_col_types = FALSE
)

hisat_par_secreted <- hisat_par %>%
  filter(secreted %in% c("yes", "manual_annotation")) %>%
  arrange(desc(Mean_TPM))

hisat_par_top30 <- hisat_par_secreted %>% slice_head(n = 30)

hisat_par_summary <- tibble(
  gland = "PAR",
  total_mean_tpm = sum(hisat_par$Mean_TPM, na.rm = TRUE),
  secreted_mean_tpm = sum(hisat_par_secreted$Mean_TPM, na.rm = TRUE),
  secreted_percent_of_total_expression = 100 * secreted_mean_tpm / total_mean_tpm,
  top30_secreted_mean_tpm = sum(hisat_par_top30$Mean_TPM, na.rm = TRUE),
  top30_percent_of_total_expression = 100 * top30_secreted_mean_tpm / total_mean_tpm,
  top30_percent_of_secreted_expression = 100 * top30_secreted_mean_tpm / secreted_mean_tpm
)

hisat_sm <- read_csv(
  "../../miscelaneous_sheets/mouse_expression/SM_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv",
  show_col_types = FALSE
)

hisat_sm_secreted <- hisat_sm %>%
  filter(secreted %in% c("yes", "manual_annotation")) %>%
  arrange(desc(Mean_TPM))

hisat_sm_top30 <- hisat_sm_secreted %>% slice_head(n = 30)

hisat_sm_summary <- tibble(
  gland = "SM",
  total_mean_tpm = sum(hisat_sm$Mean_TPM, na.rm = TRUE),
  secreted_mean_tpm = sum(hisat_sm_secreted$Mean_TPM, na.rm = TRUE),
  secreted_percent_of_total_expression = 100 * secreted_mean_tpm / total_mean_tpm,
  top30_secreted_mean_tpm = sum(hisat_sm_top30$Mean_TPM, na.rm = TRUE),
  top30_percent_of_total_expression = 100 * top30_secreted_mean_tpm / total_mean_tpm,
  top30_percent_of_secreted_expression = 100 * top30_secreted_mean_tpm / secreted_mean_tpm
)

hisat_sl <- read_csv(
  "../../miscelaneous_sheets/mouse_expression/SL_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv",
  show_col_types = FALSE
)

hisat_sl_secreted <- hisat_sl %>%
  filter(secreted %in% c("yes", "manual_annotation")) %>%
  arrange(desc(Mean_TPM))

hisat_sl_top30 <- hisat_sl_secreted %>% slice_head(n = 30)

hisat_sl_summary <- tibble(
  gland = "SL",
  total_mean_tpm = sum(hisat_sl$Mean_TPM, na.rm = TRUE),
  secreted_mean_tpm = sum(hisat_sl_secreted$Mean_TPM, na.rm = TRUE),
  secreted_percent_of_total_expression = 100 * secreted_mean_tpm / total_mean_tpm,
  top30_secreted_mean_tpm = sum(hisat_sl_top30$Mean_TPM, na.rm = TRUE),
  top30_percent_of_total_expression = 100 * top30_secreted_mean_tpm / total_mean_tpm,
  top30_percent_of_secreted_expression = 100 * top30_secreted_mean_tpm / secreted_mean_tpm
)

hisat_panc <- read_csv(
  "../../miscelaneous_sheets/mouse_expression/PANC_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv",
  show_col_types = FALSE
)

hisat_panc_secreted <- hisat_panc %>%
  filter(secreted %in% c("yes", "manual_annotation")) %>%
  arrange(desc(Mean_TPM))

hisat_panc_top30 <- hisat_panc_secreted %>% slice_head(n = 30)

hisat_panc_summary <- tibble(
  gland = "PANC",
  total_mean_tpm = sum(hisat_panc$Mean_TPM, na.rm = TRUE),
  secreted_mean_tpm = sum(hisat_panc_secreted$Mean_TPM, na.rm = TRUE),
  secreted_percent_of_total_expression = 100 * secreted_mean_tpm / total_mean_tpm,
  top30_secreted_mean_tpm = sum(hisat_panc_top30$Mean_TPM, na.rm = TRUE),
  top30_percent_of_total_expression = 100 * top30_secreted_mean_tpm / total_mean_tpm,
  top30_percent_of_secreted_expression = 100 * top30_secreted_mean_tpm / secreted_mean_tpm
)

hisat_liv <- read_csv(
  "../../miscelaneous_sheets/mouse_expression/LIV_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv",
  show_col_types = FALSE
)

hisat_liv_secreted <- hisat_liv %>%
  filter(secreted %in% c("yes", "manual_annotation")) %>%
  arrange(desc(Mean_TPM))

hisat_liv_top30 <- hisat_liv_secreted %>% slice_head(n = 30)

hisat_liv_summary <- tibble(
  gland = "LIV",
  total_mean_tpm = sum(hisat_liv$Mean_TPM, na.rm = TRUE),
  secreted_mean_tpm = sum(hisat_liv_secreted$Mean_TPM, na.rm = TRUE),
  secreted_percent_of_total_expression = 100 * secreted_mean_tpm / total_mean_tpm,
  top30_secreted_mean_tpm = sum(hisat_liv_top30$Mean_TPM, na.rm = TRUE),
  top30_percent_of_total_expression = 100 * top30_secreted_mean_tpm / total_mean_tpm,
  top30_percent_of_secreted_expression = 100 * top30_secreted_mean_tpm / secreted_mean_tpm
)

HISAT_secreted_expression_summary <- bind_rows(
  hisat_par_summary,
  hisat_sm_summary,
  hisat_sl_summary,
  hisat_panc_summary,
  hisat_liv_summary
)

write_csv(
  HISAT_secreted_expression_summary,
  "salmon_annotated_secretion_orthology/HISAT_secreted_expression_summary.csv"
)

print(HISAT_secreted_expression_summary)

