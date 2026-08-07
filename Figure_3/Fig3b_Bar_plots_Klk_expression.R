suppressPackageStartupMessages({
  library(tidyverse)
})

setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/")

df <- read_csv(
  "miscelaneous_sheets/mouse_expression/SM_mastersheet_TPMs_annotated_for_secretion.csv",
  show_col_types = FALSE
)

sample_cols <- grep("^Mous-", names(df), value = TRUE)

tpm_thresh <- df %>%
  select(all_of(sample_cols)) %>%
  mutate(across(everything(), ~ ifelse(.x >= 2, .x, 0)))

is_klk <- grepl("^Klk", df$Geneid, ignore.case = TRUE)
totals <- colSums(tpm_thresh, na.rm = TRUE)
klk_sums <- colSums(tpm_thresh[is_klk, , drop = FALSE], na.rm = TRUE)

prop_df <- tibble(
  sample = names(totals),
  totalTPM = as.numeric(totals),
  klkTPM = as.numeric(klk_sums)
) %>%
  filter(totalTPM > 0) %>%
  mutate(
    prop_KLK = klkTPM / totalTPM * 100,
    sex = ifelse(grepl("Mal", sample), "Male", "Female")
  )

p <- ggplot(prop_df, aes(sex, prop_KLK, fill = sex)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.1, height = 0, size = 2.2, alpha = 0.9) +
  labs(
    x = NULL,
    y = "Proportion of total SM expression (%)",
    title = "Relative expression of Klk genes (TPM >= 2) in mouse SM gland"
  ) +
  scale_fill_manual(values = c("Male" = "#5d7b9f", "Female" = "#d391a0")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

ggsave("Figure_3/KLK_proportion_box_SM_males_vs_females.svg", p, width = 3.2, height = 4.2)
print(p)
print(wilcox.test(prop_KLK ~ sex, data = prop_df, exact = FALSE))
