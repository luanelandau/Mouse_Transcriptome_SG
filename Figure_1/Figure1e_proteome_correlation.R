#!/usr/bin/env Rscript

# Recreate Figure 1e and the proteome/expression panels used in Figure S4.
# All inputs and outputs are restricted to Mouse_Transcriptome_SG.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(purrr)
  library(readr)
  library(stringr)
})

# ---- Paths -----------------------------------------------------------------
# This script lives in Mouse_Transcriptome_SG/Figure_1. Resolve paths from the
# script itself so it can be run from any working directory.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath(
    "~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/Figure_1/Figure1e_proteome_correlation.R",
    mustWork = TRUE
  )
}

figure_dir <- dirname(script_file)
project_dir <- dirname(figure_dir)
expression_dir <- file.path(project_dir, "miscelaneous_sheets", "mouse_expression")
output_dir <- file.path(figure_dir, "outputs_proteome_x_expression")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

proteome_file <- file.path(figure_dir, "Saliva_mouse_summary.csv")
glands <- c("PAR", "SL", "SM")
expression_files <- setNames(
  file.path(
    expression_dir,
    paste0(glands, "_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv")
  ),
  glands
)

required_files <- c(proteome_file, expression_files)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required input file(s):\n", paste(" -", missing_files, collapse = "\n"))
}
if (any(!startsWith(normalizePath(required_files), normalizePath(project_dir)))) {
  stop("All input files must be inside Mouse_Transcriptome_SG.")
}

# ---- Prepare saliva proteome ------------------------------------------------
proteome <- read_csv(proteome_file, show_col_types = FALSE)
ibaq_cols <- names(proteome)[str_detect(names(proteome), "PG\\.IBAQ$")]
if (length(ibaq_cols) < 10) {
  stop("Expected at least 10 columns ending in 'PG.IBAQ'; found ", length(ibaq_cols), ".")
}

# The original analysis used samples 1-5 as male and 6-10 as female. Samples
# 11-13 are not included in the sex-averaged saliva abundance.
male_samples <- ibaq_cols[1:5]
female_samples <- ibaq_cols[6:10]
analysis_samples <- c(male_samples, female_samples)

# Some DIA-NN cells contain semicolon-separated values. Preserve the original
# choice of using the largest numeric value in each such cell.
cell_max <- function(x) {
  vapply(x, function(cell) {
    if (is.na(cell) || str_trim(as.character(cell)) %in% c("", "NaN", "nan")) {
      return(NA_real_)
    }
    values <- as.numeric(str_extract_all(
      as.character(cell),
      "[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?"
    )[[1]])
    if (length(values)) max(values, na.rm = TRUE) else NA_real_
  }, numeric(1))
}

mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

proteome_means <- proteome %>%
  mutate(across(all_of(ibaq_cols), cell_max)) %>%
  mutate(
    gene = as.character(PG.Genes),
    ProteinAccessions = as.character(PG.ProteinAccessions)
  ) %>%
  rowwise() %>%
  mutate(
    iBAQ_mean_male = mean_or_na(c_across(all_of(male_samples))),
    iBAQ_mean_female = mean_or_na(c_across(all_of(female_samples))),
    iBAQ_mean_all_1_10 = mean_or_na(c_across(all_of(analysis_samples)))
  ) %>%
  ungroup() %>%
  select(
    ProteinAccessions, gene, PG.ProteinDescriptions, PG.Organisms,
    iBAQ_mean_male, iBAQ_mean_female, iBAQ_mean_all_1_10
  )

# Exclude ambiguous protein groups containing multiple gene identifiers, as in
# the source analysis.
proteome_unique <- proteome_means %>%
  filter(!is.na(gene), !str_detect(gene, "[;,]"))

# ---- Prepare gland expression and join --------------------------------------
secreted_values <- c("yes", "secreted", "true", "y", "1", "manual_annotation")

read_gland_expression <- function(gland) {
  read_csv(expression_files[[gland]], show_col_types = FALSE) %>%
    transmute(
      gland = gland,
      gene = as.character(Geneid),
      secreted = as.character(secreted),
      secreted_clean = str_to_lower(str_trim(secreted)),
      Mean_TPM = as.numeric(Mean_TPM),
      Mean_TPM_Male = as.numeric(Mean_TPM_Male),
      Mean_TPM_Fem = as.numeric(Mean_TPM_Fem),
      human_gene = as.character(human_gene),
      ortholog_type = as.character(ortholog_type)
    )
}

prot_expr <- map_dfr(glands, read_gland_expression) %>%
  filter(secreted_clean %in% secreted_values) %>%
  # A few gene symbols occur in multiple protein groups; retaining those rows
  # reproduces the original analysis (467 matched pairs per gland).
  left_join(proteome_unique, by = "gene", relationship = "many-to-many") %>%
  mutate(
    gland = factor(gland, levels = glands),
    log1p_TPM = log10(1 + Mean_TPM),
    log1p_iBAQ = log10(1 + iBAQ_mean_all_1_10)
  )

plot_data <- prot_expr %>%
  filter(is.finite(log1p_TPM), is.finite(log1p_iBAQ))

if (!nrow(plot_data)) stop("No finite, matched secreted gene/protein pairs were found.")

correlations <- plot_data %>%
  group_by(gland) %>%
  summarise(
    n = n(),
    spearman_rho = cor(log1p_TPM, log1p_iBAQ, method = "spearman"),
    spearman_p = cor.test(log1p_TPM, log1p_iBAQ, method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  )

write_csv(proteome_means, file.path(output_dir, "proteome_iBAQ_means_by_gene.csv"))
write_csv(plot_data, file.path(output_dir, "joined_secreted_proteome_expression_PAR_SL_SM.csv"))
write_csv(correlations, file.path(output_dir, "correlation_secreted_proteome_expression_by_gland.csv"))

# ---- Figure S4: all matched genes, top 50 labeled ---------------------------
label_top50 <- plot_data %>%
  group_by(gland) %>%
  arrange(desc(Mean_TPM), desc(iBAQ_mean_all_1_10), gene, .by_group = TRUE) %>%
  slice_head(n = 50) %>%
  ungroup()

p_s4_all <- ggplot(plot_data, aes(log1p_TPM, log1p_iBAQ)) +
  geom_point(alpha = 0.65, size = 1.6, color = "#4f8c85") +
  geom_smooth(method = "lm", se = FALSE, color = "grey50", linewidth = 0.6) +
  geom_text_repel(
    data = label_top50,
    aes(label = gene),
    size = 3, max.overlaps = Inf, box.padding = 0.35,
    point.padding = 0.2, min.segment.length = 0
  ) +
  facet_wrap(~gland, scales = "free") +
  labs(
    title = "Secreted genes (salivary glands): proteome (iBAQ) vs expression (TPM)",
    subtitle = "Labels = top 50 genes by TPM within each gland. Axes are log10(1 + value).",
    x = "log10(1 + Mean TPM) in gland",
    y = "log10(1 + mean iBAQ) in saliva proteome (samples 1-10)"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(output_dir, "proteome_vs_expression_scatter_faceted_salivary_green_labeltop50.png"),
  p_s4_all, width = 16, height = 8, dpi = 300
)

# ---- Figure S4: top 100 expressed matched genes only ------------------------
top100 <- plot_data %>%
  group_by(gland) %>%
  arrange(desc(Mean_TPM), desc(iBAQ_mean_all_1_10), gene, .by_group = TRUE) %>%
  slice_head(n = 100) %>%
  ungroup()

top100_correlations <- top100 %>%
  group_by(gland) %>%
  summarise(
    n = n(),
    spearman_rho = cor(log1p_TPM, log1p_iBAQ, method = "spearman"),
    spearman_p = cor.test(log1p_TPM, log1p_iBAQ, method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  )

p_s4_top100 <- ggplot(top100, aes(log1p_TPM, log1p_iBAQ)) +
  geom_point(alpha = 0.75, size = 1.8, color = "#4f8c85") +
  geom_smooth(method = "lm", se = FALSE, color = "grey50", linewidth = 0.6) +
  geom_text_repel(
    aes(label = gene),
    size = 2.6, max.overlaps = Inf, box.padding = 0.35,
    point.padding = 0.2, min.segment.length = 0
  ) +
  facet_wrap(~gland, scales = "free") +
  labs(
    title = "Top 100 transcripts per gland (paired protein-gene set)",
    subtitle = "Points and labels show only the top 100 genes by TPM. Axes are log10(1 + value).",
    x = "log10(1 + Mean TPM)",
    y = "log10(1 + mean iBAQ)"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(output_dir, "proteome_vs_expression_scatter_top100ONLY_salivary.png"),
  p_s4_top100, width = 16, height = 8, dpi = 300
)
#write_csv(top100, file.path(output_dir, "top100_transcripts_per_gland_paired_set.csv"))
#write_csv(top100_correlations, file.path(output_dir, "correlation_top100_transcripts_per_gland.csv"))

# ---- Figure 1e-style plots: PAR, SL and SM ----------------------------------
for (gland_name in c("PAR", "SL", "SM")) {
gland_data <- plot_data %>% filter(gland == gland_name)
gland_top50 <- gland_data %>%
  arrange(desc(Mean_TPM), desc(iBAQ_mean_all_1_10), gene) %>%
  slice_head(n = 50)
gland_top100 <- gland_data %>%
  arrange(desc(Mean_TPM), desc(iBAQ_mean_all_1_10), gene) %>%
  slice_head(n = 100)
gland_box <- gland_top100 %>%
  summarise(
    xmin = min(log1p_TPM), xmax = max(log1p_TPM),
    ymin = min(log1p_iBAQ), ymax = max(log1p_iBAQ)
  )

p_figure1e <- ggplot(gland_data, aes(log1p_TPM, log1p_iBAQ)) +
  geom_rect(
    data = gland_box,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE, color = "#4f8c85", linetype = "dashed",
    linewidth = 0.5, fill = "#a6c9bd", alpha = 0.3
  ) +
  geom_point(alpha = 0.8, size = 1.5, color = "#a6c9bd") +
  geom_smooth(method = "lm", se = FALSE, color = "grey20", linewidth = 0.5) +
  labs(
    title = paste0("Secreted genes: proteome (iBAQ) vs expression (TPM) — ", gland_name),
    subtitle = paste0(
      "Axes are log10(1 + value). Dashed box = top 100 TPM genes. ",
      "Top-50 selection saved separately. n=", nrow(gland_data)
    ),
    x = paste0("log10(1 + Mean TPM) in ", gland_name),
    y = "log10(1 + mean iBAQ) in saliva proteome"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(output_dir, paste0("scatter_secreted_TPM_vs_iBAQ_top50_box100_", gland_name, ".svg")),
  p_figure1e, width = 6.5, height = 5
)
#write_csv(gland_top50, file.path(output_dir, paste0("Figure1e_", gland_name, "_top50_labels.csv")))
#write_csv(gland_top100, file.path(output_dir, paste0("Figure1e_", gland_name, "_top100_box.csv")))
}

message("Finished. Outputs written to: ", output_dir)
