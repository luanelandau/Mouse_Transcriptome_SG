# =========================================================
# Secreted one-to-one ortholog correlations:
# PAR, SM, SL, PANC, and LIV
#
# Reports only:
# - Spearman correlation (rho)
# - Bootstrap 95% confidence interval for Spearman rho
# - Spearman p-value
#
# Outputs:
# - miscelaneous_sheets/<TISSUE>_secreted_one_to_one_unified.csv
# - miscelaneous_sheets/<TISSUE>_secreted_one_to_one_labeled_topgenes.csv
# - figures/<TISSUE>_secreted_one_to_one_correlation_labeled_human_mouse.png
# - figures/ALL_secreted_one_to_one_SPEARMAN_correlation_summary_WITH_CI.csv
# =========================================================

setwd("/Users/istarr/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(purrr)
  library(patchwork)
})

dir.create("expression_correlations/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("miscelaneous_sheets", showWarnings = FALSE, recursive = TRUE)

# Store individual tissue plots for the combined panel
correlation_plots <- list()
# ---------------------------------------------------------
# Parameters
# ---------------------------------------------------------

tissues <- c("PAR", "SM", "SL", "PANC", "LIV")

TOP_N_PER_SPECIES <- 20

CONFIDENCE_LEVEL <- 0.95

N_BOOTSTRAPS <- 5000

BOOTSTRAP_SEED <- 123

MIN_TPM <- 0

# Options:
# "both"   = TPM >= MIN_TPM in both species
# "either" = TPM >= MIN_TPM in at least one species
FILTER_MODE <- "both"


# ---------------------------------------------------------
# Input paths
# ---------------------------------------------------------

human_path <- function(tissue) {
  
  # Human salivary-gland files
  if (tissue %in% c("PAR", "SM", "SL")) {
    
    return(
      sprintf(
        "miscelaneous_sheets/human expression/%s_human_mastersheet_TPMs_annotated_with_orthology.csv",
        tissue
      )
    )
  }
  
  # Human pancreas and liver files
  tissue_lower <- tolower(tissue)
  
  sprintf(
    "miscelaneous_sheets/human expression/gene_tpm_v10_%s_filteredBy_PARgenes_annotated_means.csv",
    tissue_lower
  )
}


mouse_path <- function(tissue) {
  
  sprintf(
    "miscelaneous_sheets/mouse_expression/%s_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv",
    tissue
  )
}


# ---------------------------------------------------------
# Format p-values for graph subtitles
# ---------------------------------------------------------

format_p_value <- function(p_value) {
  
  if (is.na(p_value)) {
    return("NA")
  }
  
  if (p_value < 0.001) {
    return(format(p_value, scientific = TRUE, digits = 2))
  }
  
  sprintf("%.3f", p_value)
}


# ---------------------------------------------------------
# Bootstrap confidence interval for Spearman rho
# ---------------------------------------------------------

spearman_bootstrap_ci <- function(
    x,
    y,
    n_boot = 5000,
    conf_level = 0.95,
    seed = 123) {
  
  # Keep only complete observations
  complete_rows <- complete.cases(x, y)
  
  x <- x[complete_rows]
  y <- y[complete_rows]
  
  n <- length(x)
  
  if (n < 4) {
    
    return(
      tibble(
        ci_low = NA_real_,
        ci_high = NA_real_,
        n_successful_bootstraps = 0
      )
    )
  }
  
  set.seed(seed)
  
  bootstrap_rho <- replicate(
    n_boot,
    {
      sampled_rows <- sample(
        seq_len(n),
        size = n,
        replace = TRUE
      )
      
      suppressWarnings(
        cor(
          x[sampled_rows],
          y[sampled_rows],
          method = "spearman",
          use = "complete.obs"
        )
      )
    }
  )
  
  # Remove bootstrap samples where rho could not be calculated
  bootstrap_rho <- bootstrap_rho[
    is.finite(bootstrap_rho)
  ]
  
  if (length(bootstrap_rho) == 0) {
    
    return(
      tibble(
        ci_low = NA_real_,
        ci_high = NA_real_,
        n_successful_bootstraps = 0
      )
    )
  }
  
  alpha <- 1 - conf_level
  
  confidence_limits <- quantile(
    bootstrap_rho,
    probs = c(
      alpha / 2,
      1 - alpha / 2
    ),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )
  
  tibble(
    ci_low = confidence_limits[1],
    ci_high = confidence_limits[2],
    n_successful_bootstraps = length(bootstrap_rho)
  )
}


# ---------------------------------------------------------
# Build unified human-mouse secreted one-to-one table
# ---------------------------------------------------------

# ---------------------------------------------------------
# Build unified human-mouse secreted one-to-one table
# and filter genes by mean TPM
# ---------------------------------------------------------

build_secreted_unified <- function(tissue) {
  
  human_df <- read_csv(
    human_path(tissue),
    show_col_types = FALSE
  )
  
  mouse_df <- read_csv(
    mouse_path(tissue),
    show_col_types = FALSE
  )
  
  
  # Human secreted one-to-one orthologs
  human_secreted <- human_df %>%
    filter(
      ortholog_type == "one-to-one"
    ) %>%
    filter(
      !is.na(`Secretome location`),
      grepl(
        "Secreted",
        `Secretome location`,
        ignore.case = TRUE
      )
    ) %>%
    transmute(
      human_gene = Geneid,
      mouse_gene_expected = mouse_gene,
      TPM_human = Mean_TPM
    )
  
  
  # Mouse secreted one-to-one orthologs
  mouse_secreted <- mouse_df %>%
    filter(
      ortholog_type == "one-to-one"
    ) %>%
    filter(
      !is.na(secreted),
      !secreted %in% c(
        "No_annotation",
        "Non_secreted",
        "NA"
      )
    ) %>%
    transmute(
      human_gene = human_gene,
      mouse_gene = Geneid,
      TPM_mouse = Mean_TPM
    )
  
  
  # Join the human and mouse expression tables
  unified_unfiltered <- human_secreted %>%
    inner_join(
      mouse_secreted,
      by = "human_gene"
    ) %>%
    filter(
      is.na(mouse_gene_expected) |
        mouse_gene_expected == mouse_gene
    ) %>%
    select(
      human_gene,
      mouse_gene,
      TPM_human,
      TPM_mouse
    ) %>%
    filter(
      !is.na(TPM_human),
      !is.na(TPM_mouse),
      TPM_human >= 0,
      TPM_mouse >= 0
    ) %>%
    group_by(
      human_gene,
      mouse_gene
    ) %>%
    summarise(
      TPM_human = sum(TPM_human, na.rm = TRUE),
      TPM_mouse = sum(TPM_mouse, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  # -------------------------------------------------------
  # Apply TPM threshold
  # -------------------------------------------------------
  
  if (FILTER_MODE == "both") {
    
    unified <- unified_unfiltered %>%
      filter(
        TPM_human >= MIN_TPM,
        TPM_mouse >= MIN_TPM
      )
    
  } else if (FILTER_MODE == "either") {
    
    unified <- unified_unfiltered %>%
      filter(
        TPM_human >= MIN_TPM |
          TPM_mouse >= MIN_TPM
      )
    
  } else {
    
    stop(
      "FILTER_MODE must be either 'both' or 'either'."
    )
  }
  
  
  output_file <- sprintf(
    paste0(
      "expression_correlations/",
      "%s_secreted_one_to_one_TPM%s_%s_unified.csv"
    ),
    tissue,
    MIN_TPM,
    FILTER_MODE
  )
  
  write_csv(
    unified,
    output_file
  )
  
  message(
    sprintf(
      paste0(
        "[%s] Before TPM filtering: %d gene pairs; ",
        "after TPM >= %s filtering (%s species): %d gene pairs"
      ),
      tissue,
      nrow(unified_unfiltered),
      MIN_TPM,
      FILTER_MODE,
      nrow(unified)
    )
  )
  
  message(
    sprintf(
      "[%s] Wrote: %s",
      tissue,
      output_file
    )
  )
  
  unified
}


# ---------------------------------------------------------
# Create graph and calculate Spearman statistics
# ---------------------------------------------------------
# ---------------------------------------------------------
# Create graph and calculate Spearman statistics
# ---------------------------------------------------------

format_full_numbers <- function(x) {
  format(
    round(x),
    scientific = FALSE,
    big.mark = ",",
    trim = TRUE
  )
}

plot_spearman_correlation <- function(
    tissue,
    df,
    common_x_limits,
    common_y_limits) {
  
  if (nrow(df) < 4) {
    
    warning(
      sprintf(
        "[%s] Fewer than four gene pairs. Analysis skipped.",
        tissue
      )
    )
    
    return(
      tibble(
        tissue = tissue,
        n_pairs = nrow(df),
        spearman_rho = NA_real_,
        spearman_ci_low = NA_real_,
        spearman_ci_high = NA_real_,
        spearman_p_value = NA_real_,
        n_successful_bootstraps = NA_integer_
      )
    )
  }
  
  
  # -------------------------------------------------------
  # Spearman correlation and p-value
  # -------------------------------------------------------
  
  spearman_test <- suppressWarnings(
    cor.test(
      df$TPM_mouse,
      df$TPM_human,
      method = "spearman",
      exact = FALSE
    )
  )
  
  spearman_rho <- unname(
    spearman_test$estimate
  )
  
  spearman_p_value <- spearman_test$p.value
  
  
  # -------------------------------------------------------
  # Bootstrap confidence interval for Spearman rho
  # -------------------------------------------------------
  
  spearman_ci <- spearman_bootstrap_ci(
    x = df$TPM_mouse,
    y = df$TPM_human,
    n_boot = N_BOOTSTRAPS,
    conf_level = CONFIDENCE_LEVEL,
    seed = BOOTSTRAP_SEED
  )
  
  
  # -------------------------------------------------------
  # Regression line and confidence interval
  # -------------------------------------------------------
  
  regression_model <- lm(
    log1p(TPM_human) ~ log1p(TPM_mouse),
    data = df
  )
  
  prediction_input <- tibble(
    TPM_mouse = expm1(
      seq(
        min(log1p(df$TPM_mouse), na.rm = TRUE),
        max(log1p(df$TPM_mouse), na.rm = TRUE),
        length.out = 300
      )
    )
  )
  
  predictions <- predict(
    regression_model,
    newdata = prediction_input,
    interval = "confidence",
    level = CONFIDENCE_LEVEL
  )
  
  regression_df <- tibble(
    TPM_mouse = prediction_input$TPM_mouse,
    
    fitted = pmax(
      expm1(predictions[, "fit"]),
      0
    ),
    
    confidence_low = pmax(
      expm1(predictions[, "lwr"]),
      0
    ),
    
    confidence_high = pmax(
      expm1(predictions[, "upr"]),
      0
    )
  )
  
  
  # -------------------------------------------------------
  # Select top genes for labeling
  # -------------------------------------------------------
  
  number_to_label <- min(
    TOP_N_PER_SPECIES,
    nrow(df)
  )
  
  top_human <- df %>%
    slice_max(
      order_by = TPM_human,
      n = number_to_label,
      with_ties = FALSE
    )
  
  top_mouse <- df %>%
    slice_max(
      order_by = TPM_mouse,
      n = number_to_label,
      with_ties = FALSE
    )
  
  label_df <- bind_rows(
    top_human,
    top_mouse
  ) %>%
    distinct(
      human_gene,
      mouse_gene,
      .keep_all = TRUE
    ) %>%
    mutate(
      label_combined = paste0(
        human_gene,
        " | ",
        mouse_gene
      )
    )
  
  
  #write_csv(
  #  label_df,
  #  sprintf(
  #    paste0(
  #      "miscelaneous_sheets/",
  #      "%s_secreted_one_to_one_TPM%s_%s_labeled_topgenes.csv"
  #    ),
  #    tissue,
  #    MIN_TPM,
  #    FILTER_MODE
  #  )
  #)
  
  
  # -------------------------------------------------------
  # Graph subtitle
  # -------------------------------------------------------
  
  subtitle_text <- paste0(
    "Spearman rho = ",
    sprintf("%.3f", spearman_rho),
    #" | ",
    #round(CONFIDENCE_LEVEL * 100),
    #"% bootstrap CI [",
    #sprintf("%.3f", spearman_ci$ci_low),
    #", ",
    #sprintf("%.3f", spearman_ci$ci_high),
    #"]",
    " | p = ",
    format_p_value(spearman_p_value),
    " | n = ",
    nrow(df)
  )
  
  
  # -------------------------------------------------------
  # Create graph
  # -------------------------------------------------------
  
  correlation_plot <- ggplot(
    df,
    aes(
      x = TPM_mouse,
      y = TPM_human
    )
  ) +
    
    #geom_ribbon(
    #  data = regression_df,
    #  aes(
    #    x = TPM_mouse,
    #    ymin = confidence_low,
    #    ymax = confidence_high
    #  ),
    #  inherit.aes = FALSE,
    #  fill = "#4f8c85",
    #  alpha = 0.22
    #) +
    
    #geom_line(
    #  data = regression_df,
    #  aes(
    #    x = TPM_mouse,
    #    y = fitted
    #  ),
    #  inherit.aes = FALSE,
    #  color = "#4f8c85",
    #  linewidth = 0.9
    #) +
    
    geom_point(
      alpha = 0.8,
      size = 1.8,
      color = "#a6c9bd"
    ) +
    
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey45",
      linewidth = 0.5
    ) +
    
    ggrepel::geom_text_repel(
      data = label_df,
      aes(
        label = label_combined
      ),
      size = 2.4,
      color = "#111111",
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.15,
      min.segment.length = 0
    ) +
    
    scale_x_continuous(
      trans = "log1p",
      limits = common_x_limits,
      labels = format_full_numbers
    ) +
    
    scale_y_continuous(
      trans = "log1p",
      limits = common_y_limits,
      labels = format_full_numbers
    ) +
    
    labs(
      title = tissue,
      subtitle = subtitle_text,
      x = "Mouse Mean TPM",
      y = "Human Mean TPM"
    ) +
    
    theme_classic(
      base_size = 10
    ) +
    
    theme(
      plot.title = element_text(
        size = 12,
        face = "bold"
      ),
      plot.subtitle = element_text(
        size = 8.5,
        lineheight = 1.1
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      )
    )
  
  
  # Store plot for combined panel
  correlation_plots[[tissue]] <<- correlation_plot
  
  
  # Save individual plot
  output_plot <- sprintf(
    paste0(
      "figures/",
      "%s_secreted_one_to_one_TPM%s_%s_",
      "correlation_labeled_human_mouse.png"
    ),
    tissue,
    MIN_TPM,
    FILTER_MODE
  )
  
  #ggsave(
  #  filename = output_plot,
  #  plot = correlation_plot,
  #  width = 7.5,
  #  height = 6,
  #  dpi = 300
  #)
  #
  #message(
  #  sprintf(
  #    "[%s] Wrote: %s",
  #    tissue,
  #    output_plot
  #  )
  #)
  
  
  # -------------------------------------------------------
  # Return Spearman statistics
  # -------------------------------------------------------
  
  tibble(
    tissue = tissue,
    n_pairs = nrow(df),
    spearman_rho = spearman_rho,
    spearman_ci_low = spearman_ci$ci_low,
    spearman_ci_high = spearman_ci$ci_high,
    spearman_p_value = spearman_p_value,
    n_successful_bootstraps =
      spearman_ci$n_successful_bootstraps
  )
}

# ---------------------------------------------------------
# Save combined summary
# ---------------------------------------------------------

#summary_output <- sprintf(
#  paste0(
#    "miscelaneous_sheets/",
#    "ALL_secreted_one_to_one_TPM%s_%s_",
#    "SPEARMAN_correlation_summary_WITH_CI.csv"
#  ),
#  MIN_TPM,
#  FILTER_MODE
#)
#
#write_csv(
#  results,
#  summary_output
#)
#
#message(
#  paste0(
#    "\nCombined summary written to: ",
#    summary_output
#  )
#)
#
#print(results)
#
#
#
# ---------------------------------------------------------
# Build filtered unified tables for all tissues
# ---------------------------------------------------------

unified_tables <- set_names(
  map(
    tissues,
    build_secreted_unified
  ),
  tissues
)

#
# ---------------------------------------------------------
# Determine shared x-axis and y-axis limits
# ---------------------------------------------------------

all_unified_data <- bind_rows(
  unified_tables,
  .id = "tissue"
)

common_x_limits <- c(
  MIN_TPM,
  max(
    all_unified_data$TPM_mouse,
    na.rm = TRUE
  )
)

common_y_limits <- c(
  MIN_TPM,
  max(
    all_unified_data$TPM_human,
    na.rm = TRUE
  )
)

message(
  sprintf(
    "Shared mouse TPM range: %.2f–%.2f",
    common_x_limits[1],
    common_x_limits[2]
  )
)

message(
  sprintf(
    "Shared human TPM range: %.2f–%.2f",
    common_y_limits[1],
    common_y_limits[2]
  )
)


# ---------------------------------------------------------
# Analyze and plot all tissues
# ---------------------------------------------------------

results <- map_dfr(
  tissues,
  function(tissue) {
    
    message(
      paste0(
        "\nProcessing ",
        tissue,
        "..."
      )
    )
    
    plot_spearman_correlation(
      tissue = tissue,
      df = unified_tables[[tissue]],
      common_x_limits = common_x_limits,
      common_y_limits = common_y_limits
    )
  }
)


# ---------------------------------------------------------
# Create vertically stacked panel with shared axes
# ---------------------------------------------------------

combined_correlation_panel <- wrap_plots(
  correlation_plots[tissues],
  ncol = 1,
  guides = "collect"
) +
  
  plot_annotation(
    title = paste0(
      "Secreted one-to-one ortholog expression: Mouse vs Human"
    ),
    subtitle = paste0(
      "Mean TPM ≥ ",
      MIN_TPM,
      " in ",
      ifelse(
        FILTER_MODE == "both",
        "both species",
        "at least one species"
      ),
      "; all tissues use the same x- and y-axis ranges"
    )
  ) &
  
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )


combined_panel_output <- sprintf(
  paste0(
    "figures/",
    "ALL_secreted_one_to_one_TPM%s_%s_",
    "correlations_stacked_shared_axes.png"
  ),
  MIN_TPM,
  FILTER_MODE
)

#ggsave(
#  filename = combined_panel_output,
#  plot = combined_correlation_panel,
#  width = 8,
#  height = 26,
#  dpi = 300,
#  limitsize = FALSE
#)

message(
  paste0(
    "\nCombined stacked panel written to: ",
    combined_panel_output
  )
)



# ---------------------------------------------------------
# Create combined panel: two tissue plots per row
# ---------------------------------------------------------

combined_correlation_panel <- wrap_plots(
  correlation_plots[tissues],
  ncol = 2,
  byrow = TRUE
) +
  
  plot_annotation(
    title = paste0(
      "Secreted one-to-one ortholog expression: Mouse vs Human"
    ),
    subtitle = paste0(
      "Mean TPM ≥ ",
      MIN_TPM,
      " in ",
      ifelse(
        FILTER_MODE == "both",
        "both species",
        "at least one species"
      ),
      "; all tissues use the same x- and y-axis ranges"
    )
  ) &
  
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )


combined_panel_output <- sprintf(
  paste0(
    "expression_correlations/figures/",
    "ALL_secreted_one_to_one_TPM%s_%s_",
    "correlations_two_per_row_shared_axes.png"
  ),
  MIN_TPM,
  FILTER_MODE
)

ggsave(
  filename = combined_panel_output,
  plot = combined_correlation_panel,
  width = 15,
  height = 18,
  dpi = 300,
  limitsize = FALSE
)

message(
  paste0(
    "\nCombined two-column panel written to: ",
    combined_panel_output
  )
)

