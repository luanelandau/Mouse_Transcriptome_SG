library(dplyr)
library(ggplot2)
library(scales)

gland_colors <- c(
  "PAR" = "#d87171",
  "SM" = "#f0debe",
  "SL" = "#043a5c",
  "PANC" = "#4f8c85",
  "LIV" = "#700909"
)

prepare_tpm_data <- function(data, gland) {
  data %>%
    filter(!is.na(mean_tpm), is.finite(mean_tpm), mean_tpm > 0) %>%
    mutate(Gland = gland)
}

make_annotations <- function(data) {
  data %>%
    group_by(Gland) %>%
    group_modify(~ {
      ranked_data <- .x %>%
        arrange(desc(mean_tpm)) %>%
        mutate(cumulative_fraction = cumsum(mean_tpm) / sum(mean_tpm))

      total_expressed <- nrow(ranked_data)
      genes_for_50_percent <- which(ranked_data$cumulative_fraction >= 0.50)[1]

      data.frame(
        x = quantile(.x$mean_tpm, 0.02, na.rm = TRUE),
        label = paste0(
          "Expressed genes: ", comma(total_expressed),
          "\nGenes producing 50% of TPM: ", comma(genes_for_50_percent)
        )
      )
    }) %>%
    ungroup()
}

make_tpm_plot <- function(data, annotations) {
  ggplot(data, aes(x = mean_tpm, fill = Gland)) +
    geom_histogram(bins = 50, color = "black", linewidth = 0.25) +
    geom_text(
      data = annotations,
      aes(x = x, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1.1,
      size = 3.1
    ) +
    scale_fill_manual(values = gland_colors, guide = "none") +
    scale_x_log10(breaks = breaks_log(n = 5), labels = label_number()) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    facet_wrap(~Gland, ncol = 3) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = "Distribution of mean TPM values in mouse tissues",
      subtitle = "Each panel shows the number of expressed genes and the number of top-expressed genes accounting for 50% of total TPM",
      x = "Mean TPM (log10 scale)",
      y = "Number of genes"
    )
}

# Salmon TPM data
setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/normalization_checks/salmon_reanalysis/")

salmon_par <- read.csv("tpm/mouse_PAR_TPM_mastersheet.csv")
salmon_sm <- read.csv("tpm/mouse_SM_TPM_mastersheet.csv")
salmon_sl <- read.csv("tpm/mouse_SL_TPM_mastersheet.csv")
salmon_panc <- read.csv("tpm/mouse_PANC_TPM_mastersheet.csv")
salmon_liv <- read.csv("tpm/mouse_LIV_TPM_mastersheet.csv")

salmon_tpm <- bind_rows(
  prepare_tpm_data(salmon_par, "PAR"),
  prepare_tpm_data(salmon_sm, "SM"),
  prepare_tpm_data(salmon_sl, "SL"),
  prepare_tpm_data(salmon_panc, "PANC"),
  prepare_tpm_data(salmon_liv, "LIV")
) %>%
  mutate(Gland = factor(Gland, levels = c("PAR", "SM", "SL", "PANC", "LIV")))

salmon_annotations <- make_annotations(salmon_tpm)
salmon_distribution_plot <- make_tpm_plot(salmon_tpm, salmon_annotations)

print(salmon_distribution_plot)

ggsave(
  "~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/normalization_checks/salmon_reanalysis/SALMON_mouse_TPM_distributions_by_tissue_shared_axes.png",
  salmon_distribution_plot, width = 15, height = 9, dpi = 300
)

# HISAT TPM data
setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/miscelaneous_sheets/")

hisat_par <- read.csv("mouse_expression/PAR_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv") %>% rename(mean_tpm = Mean_TPM)
hisat_sm <- read.csv("mouse_expression/SM_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv") %>% rename(mean_tpm = Mean_TPM)
hisat_sl <- read.csv("mouse_expression/SL_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv") %>% rename(mean_tpm = Mean_TPM)
hisat_panc <- read.csv("mouse_expression/PANC_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv") %>% rename(mean_tpm = Mean_TPM)
hisat_liv <- read.csv("mouse_expression/LIV_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv") %>% rename(mean_tpm = Mean_TPM)

hisat_tpm <- bind_rows(
  prepare_tpm_data(hisat_par, "PAR"),
  prepare_tpm_data(hisat_sm, "SM"),
  prepare_tpm_data(hisat_sl, "SL"),
  prepare_tpm_data(hisat_panc, "PANC"),
  prepare_tpm_data(hisat_liv, "LIV")
) %>%
  mutate(Gland = factor(Gland, levels = c("PAR", "SM", "SL", "PANC", "LIV")))

hisat_annotations <- make_annotations(hisat_tpm)
hisat_distribution_plot <- make_tpm_plot(hisat_tpm, hisat_annotations)

print(hisat_distribution_plot)

ggsave(
  "~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/normalization_checks/salmon_reanalysis/HISAT_mouse_TPM_distributions_by_tissue_shared_axes.png",
  hisat_distribution_plot, width = 15, height = 9, dpi = 300
)
