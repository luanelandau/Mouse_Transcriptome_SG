# =========================================================
# Compare mouse-human Spearman correlations among:
# PAR, SM, SL, PANC, and LIV
#
# Method:
# - For each pair of tissues, retain only ortholog pairs
#   present in both tissues
# - Calculate the two Spearman correlations using the same genes
# - Resample those shared genes with replacement
# - Calculate the bootstrap distribution of the difference
#
# Outputs:
# figures/ALL_5_tissues_spearman_bootstrap_comparisons.csv
# figures/ALL_5_tissues_spearman_bootstrap_values.csv
# figures/ALL_5_tissues_spearman_bootstrap_comparisons.png
# =========================================================

setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(tidyr)
})

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

set.seed(123)

N_BOOTSTRAPS <- 10000


# ---------------------------------------------------------
# Read the five unified expression tables
# ---------------------------------------------------------

par <- read_csv(
  "expression_correlations/PAR_secreted_one_to_one_TPM0_both_unified.csv",
  show_col_types = FALSE
)

sm <- read_csv(
  "expression_correlations/SM_secreted_one_to_one_TPM0_both_unified.csv",
  show_col_types = FALSE
)

sl <- read_csv(
  "expression_correlations/SL_secreted_one_to_one_TPM0_both_unified.csv",
  show_col_types = FALSE
)

panc <- read_csv(
  "expression_correlations/PANC_secreted_one_to_one_TPM0_both_unified.csv",
  show_col_types = FALSE
)

liv <- read_csv(
  "expression_correlations/LIV_secreted_one_to_one_TPM0_both_unified.csv",
  show_col_types = FALSE
)


# ---------------------------------------------------------
# Store the tables in one named list
# This is not a function; it just makes the loop simpler
# ---------------------------------------------------------

tissue_tables <- list(
  PAR = par,
  SM = sm,
  SL = sl,
  PANC = panc,
  LIV = liv
)


# ---------------------------------------------------------
# Define all 10 pairwise comparisons
# ---------------------------------------------------------

tissue_comparisons <- combn(
  names(tissue_tables),
  2,
  simplify = FALSE
)

tissue_comparisons

# The comparisons will be:
# PAR vs SM
# PAR vs SL
# PAR vs PANC
# PAR vs LIV
# SM vs SL
# SM vs PANC
# SM vs LIV
# SL vs PANC
# SL vs LIV
# PANC vs LIV


# ---------------------------------------------------------
# Empty tables to store results
# ---------------------------------------------------------

comparison_summary <- tibble()

all_bootstrap_values <- tibble()


# =========================================================
# Run each tissue comparison
# =========================================================

for (comparison_number in seq_along(tissue_comparisons)) {
  
  tissue_1 <- tissue_comparisons[[comparison_number]][1]
  tissue_2 <- tissue_comparisons[[comparison_number]][2]
  
  message(
    "Running comparison: ",
    tissue_1,
    " vs ",
    tissue_2
  )
  
  
  # -------------------------------------------------------
  # Get the two tissue tables
  # -------------------------------------------------------
  
  data_1 <- tissue_tables[[tissue_1]] %>%
    transmute(
      human_gene,
      mouse_gene,
      TPM_human_1 = TPM_human,
      TPM_mouse_1 = TPM_mouse
    )
  
  data_2 <- tissue_tables[[tissue_2]] %>%
    transmute(
      human_gene,
      mouse_gene,
      TPM_human_2 = TPM_human,
      TPM_mouse_2 = TPM_mouse
    )
  
  
  # -------------------------------------------------------
  # Keep only ortholog pairs found in both tissues
  # -------------------------------------------------------
  
  paired_data <- data_1 %>%
    inner_join(
      data_2,
      by = c("human_gene", "mouse_gene")
    ) %>%
    filter(
      complete.cases(
        TPM_human_1,
        TPM_mouse_1,
        TPM_human_2,
        TPM_mouse_2
      )
    )
  
  
  n_shared_pairs <- nrow(paired_data)
  
  message(
    "Shared ortholog pairs: ",
    n_shared_pairs
  )
  
  
  # -------------------------------------------------------
  # Observed Spearman correlations
  # -------------------------------------------------------
  
  observed_r_1 <- suppressWarnings(
    cor(
      paired_data$TPM_mouse_1,
      paired_data$TPM_human_1,
      method = "spearman"
    )
  )
  
  observed_r_2 <- suppressWarnings(
    cor(
      paired_data$TPM_mouse_2,
      paired_data$TPM_human_2,
      method = "spearman"
    )
  )
  
  
  # Difference is always tissue 1 minus tissue 2
  
  observed_difference <- observed_r_1 - observed_r_2
  
  
  # -------------------------------------------------------
  # Empty vectors for the bootstrap values
  # -------------------------------------------------------
  
  bootstrap_r_1 <- numeric(N_BOOTSTRAPS)
  bootstrap_r_2 <- numeric(N_BOOTSTRAPS)
  bootstrap_difference <- numeric(N_BOOTSTRAPS)
  
  
  # -------------------------------------------------------
  # Resample shared ortholog pairs
  # -------------------------------------------------------
  
  for (i in 1:N_BOOTSTRAPS) {
    
    sampled_rows <- sample(
      seq_len(n_shared_pairs),
      size = n_shared_pairs,
      replace = TRUE
    )
    
    sampled_data <- paired_data[sampled_rows, ]
    
    
    bootstrap_r_1[i] <- suppressWarnings(
      cor(
        sampled_data$TPM_mouse_1,
        sampled_data$TPM_human_1,
        method = "spearman"
      )
    )
    
    
    bootstrap_r_2[i] <- suppressWarnings(
      cor(
        sampled_data$TPM_mouse_2,
        sampled_data$TPM_human_2,
        method = "spearman"
      )
    )
    
    
    bootstrap_difference[i] <-
      bootstrap_r_1[i] - bootstrap_r_2[i]
  }
  
  
  # -------------------------------------------------------
  # Put bootstrap values into a table
  # -------------------------------------------------------
  
  bootstrap_results <- tibble(
    comparison = paste(tissue_1, "vs", tissue_2),
    tissue_1 = tissue_1,
    tissue_2 = tissue_2,
    bootstrap_iteration = 1:N_BOOTSTRAPS,
    r_tissue_1 = bootstrap_r_1,
    r_tissue_2 = bootstrap_r_2,
    correlation_difference = bootstrap_difference
  ) %>%
    filter(
      !is.na(correlation_difference),
      is.finite(correlation_difference)
    )
  
  
  # -------------------------------------------------------
  # Bootstrap 95% confidence interval
  # -------------------------------------------------------
  
  confidence_interval <- quantile(
    bootstrap_results$correlation_difference,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  
  # -------------------------------------------------------
  # Two-sided bootstrap p-value
  # -------------------------------------------------------
  
  proportion_zero_or_below <- (
    sum(
      bootstrap_results$correlation_difference <= 0
    ) + 1
  ) / (
    nrow(bootstrap_results) + 1
  )
  
  
  proportion_zero_or_above <- (
    sum(
      bootstrap_results$correlation_difference >= 0
    ) + 1
  ) / (
    nrow(bootstrap_results) + 1
  )
  
  
  bootstrap_p_value <- 2 * min(
    proportion_zero_or_below,
    proportion_zero_or_above
  )
  
  
  bootstrap_p_value <- min(
    bootstrap_p_value,
    1
  )
  
  
  # -------------------------------------------------------
  # Store the summary for this comparison
  # -------------------------------------------------------
  
  current_summary <- tibble(
    comparison = paste(tissue_1, "vs", tissue_2),
    tissue_1 = tissue_1,
    tissue_2 = tissue_2,
    n_shared_pairs = n_shared_pairs,
    r_tissue_1 = observed_r_1,
    r_tissue_2 = observed_r_2,
    correlation_difference = observed_difference,
    lower_95_CI = unname(confidence_interval[1]),
    upper_95_CI = unname(confidence_interval[2]),
    bootstrap_p_value = bootstrap_p_value
  )
  
  
  comparison_summary <- bind_rows(
    comparison_summary,
    current_summary
  )
  
  
  all_bootstrap_values <- bind_rows(
    all_bootstrap_values,
    bootstrap_results
  )
}


# =========================================================
# Adjust p-values across all 10 comparisons
# =========================================================

comparison_summary <- comparison_summary %>%
  mutate(
    adjusted_p_value = p.adjust(
      bootstrap_p_value,
      method = "BH"
    ),
    
    higher_correlation = case_when(
      correlation_difference > 0 ~ tissue_1,
      correlation_difference < 0 ~ tissue_2,
      TRUE ~ "equal"
    ),
    
    confidence_interval_excludes_zero =
      lower_95_CI > 0 |
      upper_95_CI < 0,
    
    significant_unadjusted =
      bootstrap_p_value < 0.05 &
      confidence_interval_excludes_zero,
    
    significant_adjusted =
      adjusted_p_value < 0.05 &
      confidence_interval_excludes_zero
  ) %>%
  arrange(
    adjusted_p_value
  )


# ---------------------------------------------------------
# Save the results
# ---------------------------------------------------------

write_csv(
  comparison_summary,
  "expression_correlations/ALL_5_tissues_spearman_bootstrap_comparisons.csv"
)

write_csv(
  all_bootstrap_values,
  "expression_correlations/ALL_5_tissues_spearman_bootstrap_values.csv"
)


# ---------------------------------------------------------
# Print the complete summary
# ---------------------------------------------------------

print(
  comparison_summary,
  n = Inf
)


# =========================================================
# Plot observed differences and bootstrap confidence intervals
# =========================================================

plot_data <- comparison_summary %>%
  mutate(
    comparison = factor(
      comparison,
      levels = rev(comparison)
    )
  )


correlation_comparison_plot <- ggplot(
  plot_data,
  aes(
    x = correlation_difference,
    y = comparison
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.6
  ) +
  geom_errorbarh(
    aes(
      xmin = lower_95_CI,
      xmax = upper_95_CI
    ),
    height = 0.15,
    linewidth = 0.7
  ) +
  geom_point(
    aes(
      shape = significant_adjusted
    ),
    size = 3,
    color = "#4f8c85"
  ) +
  scale_shape_manual(
    values = c(
      "FALSE" = 1,
      "TRUE" = 16
    ),
    labels = c(
      "FALSE" = "Not significant",
      "TRUE" = "Significant"
    )
  ) +
  labs(
    title = "Comparison of mouse–human Spearman correlations",
    subtitle = paste0(
      "Difference = correlation of first tissue − ",
      "correlation of second tissue"
    ),
    x = "Difference in Spearman correlation",
    y = NULL,
    shape = "BH-adjusted result"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

correlation_comparison_plot
ggsave(
  "expression_correlations/figures/ALL_5_tissues_spearman_bootstrap_comparisons.png",
  correlation_comparison_plot,
  width = 8,
  height = 6,
  dpi = 300
)
