library(tidyverse)
library(ggrepel)

setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/")

# -----------------------------
# Paths
# -----------------------------
# This file is an output from the proteome analysis that was performed
# externally and returned to us by the proteomics analysis provider.
saliva_diff_path <- "Figure_2/Saliva_mouse_sex_diff.csv"
means_path <- "Figure_1/outputs_proteome_x_expression/proteome_iBAQ_means_by_gene.csv"

deg_dir <- "deseq2/deseq2_results"
par_path <- file.path(deg_dir, "DESeq2_sex_PAR.csv")
sl_path  <- file.path(deg_dir, "DESeq2_sex_SL.csv")
sm_path  <- file.path(deg_dir, "DESeq2_sex_SM.csv")

# Save all Figure 2c files directly in the Figure_2 folder.
output_dir <- "Figure_2"

# -----------------------------
# Cutoffs
# -----------------------------
q_cutoff_saliva   <- 0.05 # saliva proteome sex-dimorphic (Qvalue)
padj_cutoff_gland <- 0.05 # transcriptome DEG threshold
lfc_cutoff_gland  <- 1    # transcriptome DEG threshold: |log2FC| > 1

# -----------------------------
# Colors / shapes
# -----------------------------
color_map <- c(
  "PAR" = "#d87171",
  "SL"  = "#043a5c",
  "SM"  = "#4f8c85",
  "Multiple" = "black",
  "Saliva_only" = "grey80"
)

shape_map <- c(
  "0" = 4,   # not DEG in any gland among PAR/SL/SM
  "1" = 16,  # DEG in 1 gland
  "2" = 17,  # DEG in 2 glands
  "3" = 15   # DEG in 3 glands
)

# -----------------------------
# Helper: read DESeq2 sex CSV (first col is gene) and apply DEG definition
# DEG = padj < 0.05 AND |log2FC| > 1
# -----------------------------
read_deseq_sex <- function(path, gland_label) {
  df <- read_csv(path, show_col_types = FALSE)
  colnames(df)[1] <- "gene"

  df %>%
    transmute(
      gene = as.character(gene),
      gland = gland_label,
      log2FoldChange = suppressWarnings(as.numeric(log2FoldChange)),
      padj = suppressWarnings(as.numeric(padj))
    ) %>%
    filter(
      is.finite(padj),
      padj < padj_cutoff_gland,
      is.finite(log2FoldChange),
      abs(log2FoldChange) > lfc_cutoff_gland
    )
}

# -----------------------------
# 1) Saliva proteome sex-dimorphic genes (Qvalue < 0.05)
# -----------------------------
df_saliva_diff_raw <- read_csv(saliva_diff_path, show_col_types = FALSE)

df_saliva_DA_genes <- df_saliva_diff_raw %>%
  transmute(
    gene = str_trim(as.character(Genes)),
    avg_log2_ratio = suppressWarnings(as.numeric(`AVG Log2 Ratio`)),
    abs_avg_log2_ratio = suppressWarnings(as.numeric(`Absolute AVG Log2 Ratio`)),
    qvalue = suppressWarnings(as.numeric(Qvalue))
  ) %>%
  # Drop protein groups assigned to multiple genes (for example, Klk1;Klk1b24).
  filter(
    !is.na(gene),
    gene != "",
    gene != "NaN",
    !str_detect(gene, ";"),
    is.finite(qvalue),
    qvalue < q_cutoff_saliva
  ) %>%
  group_by(gene) %>%
  arrange(desc(abs_avg_log2_ratio), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

# -----------------------------
# 2) Gland transcriptome sex-DEGs (PAR/SL/SM)
# -----------------------------
deg_par <- read_deseq_sex(par_path, "PAR")
deg_sl  <- read_deseq_sex(sl_path,  "SL")
deg_sm  <- read_deseq_sex(sm_path,  "SM")

deg_all <- bind_rows(deg_par, deg_sl, deg_sm) %>%
  distinct(gene, gland)

# Build membership map: per gene, how many glands?
gland_membership <- deg_all %>%
  group_by(gene) %>%
  summarise(
    glands = list(sort(unique(gland))),
    n_glands = as.integer(length(unique(gland))),
    .groups = "drop"
  ) %>%
  mutate(
    single_gland = if_else(n_glands == 1L, map_chr(glands, 1), NA_character_)
  )

# -----------------------------
# 3) Proteome means (impute missing = 0; compute log10(1 + iBAQ))
# -----------------------------
df_means <- read_csv(means_path, show_col_types = FALSE)

df_corr <- df_means %>%
  transmute(
    gene = str_trim(as.character(gene)),
    iBAQ_mean_male   = suppressWarnings(as.numeric(iBAQ_mean_male)),
    iBAQ_mean_female = suppressWarnings(as.numeric(iBAQ_mean_female))
  ) %>%
  mutate(
    iBAQ_mean_male   = if_else(is.na(iBAQ_mean_male), 0, iBAQ_mean_male),
    iBAQ_mean_female = if_else(is.na(iBAQ_mean_female), 0, iBAQ_mean_female),
    log10p1_male   = log10(1 + iBAQ_mean_male),
    log10p1_female = log10(1 + iBAQ_mean_female)
  )

# -----------------------------
# 4) Final plotting table:
# ONLY saliva sex-dimorphic proteins, annotated by transcriptome gland DEGs
# -----------------------------
# Joining the DA statistics here only once avoids duplicated .x/.y columns.
df_plot <- df_corr %>%
  inner_join(df_saliva_DA_genes, by = "gene") %>%
  left_join(
    gland_membership %>% select(gene, n_glands, single_gland),
    by = "gene"
  ) %>%
  mutate(
    n_glands = if_else(is.na(n_glands), 0L, n_glands),
    overlap_class = as.character(n_glands),
    gland_color_group = case_when(
      n_glands == 0L ~ "Saliva_only",
      n_glands == 1L ~ single_gland,
      n_glands >= 2L ~ "Multiple",
      TRUE ~ "Saliva_only"
    )
  )

# Label all gland transcriptome sex-DEGs present in the plotted proteins.
label_df <- df_plot %>%
  filter(gene %in% c(deg_sm$gene, deg_par$gene, deg_sl$gene))

# -----------------------------
# Figure 1: correlation plot used in the paper
# -----------------------------
p_gland_DEGs <- ggplot(df_plot, aes(x = log10p1_male, y = log10p1_female)) +
  geom_point(
    aes(color = gland_color_group, shape = overlap_class),
    alpha = 0.80,
    size = 2.5
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    color = "grey50",
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  geom_text_repel(
    data = label_df,
    aes(x = log10p1_male, y = log10p1_female, label = gene),
    inherit.aes = FALSE,
    color = "black",
    size = 3.1,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.2,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  scale_color_manual(values = color_map, name = "Transcriptome sex-DEG") +
  scale_shape_manual(
    values = shape_map,
    name = "DEG overlap\n(PAR/SL/SM)",
    labels = c(
      "0" = "Not DEG in PAR/SL/SM",
      "1" = "DEG in 1 gland",
      "2" = "DEG in 2 glands",
      "3" = "DEG in 3 glands"
    )
  ) +
  labs(
    x = "Male saliva protein abundance (log10(1+ iBAQ))",
    y = "Female saliva protein abundance (log10(1+ iBAQ))",
    title = "Sex-dimorphic saliva proteins annotated by transcriptome sex-DEGs",
    subtitle = paste0(
      "Saliva sex-dimorphic: Qvalue < ", q_cutoff_saliva,
      "; Transcriptome DEG: padj < ", padj_cutoff_gland,
      " and |log2FC| > ", lfc_cutoff_gland
    )
  ) +
  theme_classic()

print(p_gland_DEGs)

ggsave(
  file.path(output_dir, "correlation_of_significant_male_vs_female_proteins_with_contribution_from_each_gland.png"),
  p_gland_DEGs, height = 10, width = 8.5, dpi = 300
)
ggsave(
  file.path(output_dir, "correlation_of_significant_male_vs_female_proteins_with_contribution_from_each_gland.svg"),
  p_gland_DEGs, height = 10, width = 8.5
)

# -----------------------------
# Diagnostics
# -----------------------------
cat("\n--- Diagnostics ---\n")
cat("Saliva sex-dimorphic proteins (unique genes): ", n_distinct(df_saliva_DA_genes$gene), "\n", sep = "")
cat("PAR transcriptome sex-DEGs: ", nrow(deg_par), "\n", sep = "")
cat("SL  transcriptome sex-DEGs: ", nrow(deg_sl),  "\n", sep = "")
cat("SM  transcriptome sex-DEGs: ", nrow(deg_sm),  "\n", sep = "")
cat("Plotted points (saliva DA with proteome means): ", nrow(df_plot), "\n", sep = "")
cat("Overlap counts among plotted genes:\n")
print(df_plot %>% count(overlap_class, sort = TRUE))
cat("Color-group counts among plotted genes:\n")
print(df_plot %>% count(gland_color_group, sort = TRUE))

# -----------------------------
# Figure 2: bar plot used in the assembled figure
# -----------------------------
df_bar <- df_plot %>%
  distinct(gene, .keep_all = TRUE) %>%
  filter(gland_color_group %in% c("Saliva_only", "SM", "PAR", "SL")) %>%
  count(gland_color_group) %>%
  mutate(
    x_label = factor(
      recode(
        gland_color_group,
        "Saliva_only" = "DA in protein only",
        "SM" = "SM DEG",
        "PAR" = "PAR DEG",
        "SL" = "SL DEG"
      ),
      levels = c("DA in protein only", "SM DEG", "PAR DEG", "SL DEG")
    ),
    fill_group = if_else(gland_color_group == "SM", "SM", "Other")
  )

p_bar <- ggplot(df_bar, aes(x = x_label, y = n, fill = fill_group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), vjust = -0.4, size = 3.8, color = "black") +
  scale_fill_manual(
    values = c("SM" = "#d87171", "Other" = "grey80"),
    guide = "none"
  ) +
  theme_classic()

print(p_bar)

ggsave(
  file.path(output_dir, "bars_number_of_DEGs_genes_per_gland_in_saliva.png"),
  p_bar, height = 4, width = 3.2, dpi = 300
)
ggsave(
  file.path(output_dir, "bars_number_of_DEGs_genes_per_gland_in_saliva.svg"),
  p_bar, height = 4, width = 3.2
)
