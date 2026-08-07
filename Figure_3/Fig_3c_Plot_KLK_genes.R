suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/")

df <- read_csv(
  "miscelaneous_sheets/gene_expression_matrix_C57_CD1_TPMs.csv",
  show_col_types = FALSE
)

klk_sm_sig <- c(
  "Klk1", "Klk15", "Klk1b1", "Klk1b11", "Klk1b16", "Klk1b21",
  "Klk1b22", "Klk1b24", "Klk1b26", "Klk1b27", "Klk1b3", "Klk1b4",
  "Klk1b5", "Klk1b7-ps", "Klk1b8", "Klk1b9"
)

klk_genes <- filter(df, Geneid %in% klk_sm_sig)
sm_male_cols <- grep("SM-.*Mal", names(klk_genes), value = TRUE)
sm_female_cols <- grep("SM-.*Fem", names(klk_genes), value = TRUE)

summ_log <- klk_genes %>%
  transmute(
    Geneid,
    female_mean = rowMeans(across(all_of(sm_female_cols)), na.rm = TRUE),
    female_sd = apply(across(all_of(sm_female_cols)), 1, sd, na.rm = TRUE),
    male_mean = rowMeans(across(all_of(sm_male_cols)), na.rm = TRUE),
    male_sd = apply(across(all_of(sm_male_cols)), 1, sd, na.rm = TRUE)
  ) %>%
  mutate(
    x_log = log10(female_mean + 1),
    y_log = log10(male_mean + 1),
    x_sd_log_left = log10(pmax(female_mean - female_sd, 0) + 1),
    x_sd_log_right = log10(female_mean + female_sd + 1),
    y_sd_log_low = log10(pmax(male_mean - male_sd, 0) + 1),
    y_sd_log_high = log10(male_mean + male_sd + 1)
  )

p <- ggplot(summ_log, aes(x_log, y_log, label = Geneid)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_abline(slope = 1, intercept = log10(1.5), linetype = "dotted") +
  geom_abline(slope = 1, intercept = log10(5), linetype = "dotdash") +
  geom_segment(
    aes(x = x_sd_log_left, xend = x_sd_log_right, y = y_log, yend = y_log),
    alpha = 0.6
  ) +
  geom_segment(
    aes(x = x_log, xend = x_log, y = y_sd_log_low, yend = y_sd_log_high),
    alpha = 0.6
  ) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_text_repel(
    min.segment.length = 0,
    box.padding = 0.25,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  labs(
    title = "Klk gene expression in Submandibular gland (log10 scale)",
    subtitle = paste(
      "Dots = log10(TPM+1) mean; crosshairs = +/-1 SD (log10 scale);",
      "lines = 1:1 (dashed), 1:1.5 (dotted), 1:5 (dotdash)"
    ),
    x = "Female log10(TPM + 1) (mean, SM)",
    y = "Male log10(TPM + 1) (mean, SM)"
  ) +
  theme_minimal(base_size = 12)

ggsave("Figure_3/Klk_SM_Female_vs_Male_log10_withSD.svg", p, width = 7, height = 6)
print(p)
