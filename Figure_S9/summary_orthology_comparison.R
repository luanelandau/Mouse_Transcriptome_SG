# Summarize Ensembl and Jackson Laboratory orthology classifications

library(readr)
library(dplyr)
library(tidyr)

# Paths -------------------------------------------------------------------

project_dir <- "/Users/luane/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG"
figure_dir <- file.path(project_dir, "Figure_S9")

ensembl_file <- file.path(figure_dir, "mouse_human_orthologs_ENS_v115.csv")
jackson_file <- file.path(project_dir, "miscelaneous_sheets", "orthologs_with_classification.csv")
sm_file <- file.path(
  project_dir,
  "miscelaneous_sheets",
  "mouse_expression",
  "SM_mastersheet_TPMs_annotated_for_secretion.csv"
)

output_file <- file.path(figure_dir, "orthology_summary_ensembl_vs_jackson.csv")

# Read and simplify the orthology classifications -------------------------

ensembl <- read_csv(ensembl_file, show_col_types = FALSE) %>%
  transmute(
    mouse_gene = external_gene_name,
    ortholog_type = recode(
      hsapiens_homolog_orthology_type,
      ortholog_one2one = "one-to-one",
      ortholog_one2many = "one-to-many",
      ortholog_many2many = "many-to-many",
      .default = NA_character_
    )
  ) %>%
  filter(!is.na(mouse_gene), mouse_gene != "", !is.na(ortholog_type)) %>%
  distinct(mouse_gene, ortholog_type)

jackson <- read_csv(jackson_file, show_col_types = FALSE) %>%
  transmute(
    mouse_gene,
    ortholog_type = recode(
      ortholog_type,
      `one-to-one` = "one-to-one",
      `one-to-many (human-to-mouse)` = "one-to-many",
      `one-to-many (mouse-to-human)` = "one-to-many",
      `many-to-many` = "many-to-many",
      .default = NA_character_
    )
  ) %>%
  filter(!is.na(mouse_gene), mouse_gene != "", !is.na(ortholog_type)) %>%
  distinct(mouse_gene, ortholog_type)

sm_genes <- read_csv(sm_file, show_col_types = FALSE) %>%
  transmute(mouse_gene = Geneid) %>%
  filter(!is.na(mouse_gene), mouse_gene != "") %>%
  distinct()

# Count unique mouse genes in each category -------------------------------

ensembl_counts <- ensembl %>%
  count(ortholog_type, name = "ensembl_genes")

jackson_counts <- jackson %>%
  count(ortholog_type, name = "jackson_genes")

concordant_counts <- inner_join(
  ensembl,
  jackson,
  by = c("mouse_gene", "ortholog_type")
) %>%
  count(ortholog_type, name = "concordant_genes")

lineage_specific <- tibble(
  ortholog_type = "lineage-specific",
  ensembl_genes = nrow(anti_join(sm_genes, ensembl, by = "mouse_gene")),
  jackson_genes = nrow(anti_join(sm_genes, jackson, by = "mouse_gene")),
  concordant_genes = nrow(
    anti_join(sm_genes, ensembl, by = "mouse_gene") %>%
      semi_join(anti_join(sm_genes, jackson, by = "mouse_gene"), by = "mouse_gene")
  )
)

# Create and save one summarized comparison sheet -------------------------

orthology_summary <- full_join(
  ensembl_counts,
  jackson_counts,
  by = "ortholog_type"
) %>%
  full_join(concordant_counts, by = "ortholog_type") %>%
  complete(
    ortholog_type = c("one-to-one", "one-to-many", "many-to-many"),
    fill = list(ensembl_genes = 0L, jackson_genes = 0L, concordant_genes = 0L)
  ) %>%
  bind_rows(lineage_specific) %>%
  mutate(
    difference_jackson_minus_ensembl = jackson_genes - ensembl_genes,
    database_with_more_genes = case_when(
      difference_jackson_minus_ensembl > 0 ~ "Jackson",
      difference_jackson_minus_ensembl < 0 ~ "Ensembl",
      TRUE ~ "Equal"
    ),
    ortholog_type = factor(
      ortholog_type,
      levels = c("one-to-one", "one-to-many", "many-to-many", "lineage-specific")
    )
  ) %>%
  arrange(ortholog_type) %>%
  mutate(ortholog_type = as.character(ortholog_type))

print(orthology_summary, n = Inf)
write_csv(orthology_summary, output_file)
