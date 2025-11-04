# =========================================================
# Secreted one-to-one orthologs: build tables + correlations for PAR, SM, SL
# - Writes: figures/<GLAND>_secreted_one_to_one_unified.csv
# - Plots:  figures/<GLAND>_secreted_one_to_one_correlation_labeled_human_mouse.png
# - Summary correlations: figures/ALL_secreted_one_to_one_correlation_summary.csv
# =========================================================

setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
  library(purrr)
})

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# ---------- Parameters ----------
glands <- c("PAR","SM","SL")
TOP_N_PER_SPECIES <- 20

# ---------- Helpers ----------
human_path <- function(g) sprintf("miscelaneous_sheets/human expression/%s_human_mastersheet_TPMs_annotated_with_orthology.csv", g)
mouse_path <- function(g) sprintf("miscelaneous_sheets/mouse_expression/%s_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv", g)

# Build unified secreted 1:1 for a gland
build_secreted_unified <- function(g) {
  human_df <- read_csv(human_path(g), show_col_types = FALSE)
  mouse_df <- read_csv(mouse_path(g), show_col_types = FALSE)
  
  # Human: one-to-one & secreted (string contains "Secreted")
  human_1to1_sec <- human_df %>%
    filter(ortholog_type == "one-to-one") %>%
    filter(!is.na(`Secretome location`) & grepl("Secreted", `Secretome location`, ignore.case = TRUE)) %>%
    transmute(
      human_gene          = Geneid,
      mouse_gene_expected = mouse_gene,
      TPM_human           = Mean_TPM
    )
  
  # Mouse: one-to-one & secreted (exclude "No_annotation"/"Non_secreted"/NA)
  mouse_1to1_sec <- mouse_df %>%
    filter(ortholog_type == "one-to-one") %>%
    filter(!is.na(secreted) & !secreted %in% c("No_annotation", "Non_secreted", "NA")) %>%
    transmute(
      human_gene = human_gene,   # key
      mouse_gene = Geneid,
      TPM_mouse  = Mean_TPM
    )
  
  # Join & enforce expected mouse match if provided
  unified <- human_1to1_sec %>%
    inner_join(mouse_1to1_sec, by = "human_gene") %>%
    filter(is.na(mouse_gene_expected) | mouse_gene_expected == mouse_gene) %>%
    select(human_gene, mouse_gene, TPM_human, TPM_mouse) %>%
    group_by(human_gene, mouse_gene) %>%
    summarise(
      TPM_human = sum(TPM_human, na.rm = TRUE),
      TPM_mouse = sum(TPM_mouse, na.rm = TRUE),
      .groups = "drop"
    )
  
  out_csv <- sprintf("figures/%s_secreted_one_to_one_unified.csv", g)
  write_csv(unified, out_csv)
  message(sprintf("[%s] Wrote: %s (n=%d)", g, out_csv, nrow(unified)))
  unified
}

# Plot correlation with labels for top-N human and top-N mouse
plot_corr <- function(g, df) {
  # Correlations
  r_pearson_raw  <- suppressWarnings(cor(df$TPM_mouse, df$TPM_human,  method = "pearson"))
  r_spearman_raw <- suppressWarnings(cor(df$TPM_mouse, df$TPM_human,  method = "spearman"))
  r_pearson_log  <- suppressWarnings(cor(log1p(df$TPM_mouse), log1p(df$TPM_human), method = "pearson"))
  
  # Fit regression in log1p space and back-transform for plotting
  fit <- lm(log1p(TPM_human) ~ log1p(TPM_mouse), data = df)
  coef_fit <- coef(fit)
  xlog_seq <- seq(min(log1p(df$TPM_mouse), na.rm = TRUE),
                  max(log1p(df$TPM_mouse), na.rm = TRUE),
                  length.out = 300)
  ylog_hat <- coef_fit[1] + coef_fit[2] * xlog_seq
  pred_df  <- data.frame(
    x = pmax(expm1(xlog_seq), 0),
    y = pmax(expm1(ylog_hat), 0)
  )
  
  # Top labels per species
  top_human <- df %>% slice_max(order_by = TPM_human, n = TOP_N_PER_SPECIES) %>% mutate(label_source = "human_top")
  top_mouse <- df %>% slice_max(order_by = TPM_mouse, n = TOP_N_PER_SPECIES) %>% mutate(label_source = "mouse_top")
  
  label_df <- bind_rows(top_human, top_mouse) %>%
    distinct(human_gene, .keep_all = TRUE) %>%
    mutate(label_combined = paste0(human_gene, " | ", mouse_gene))
  
  # Save which genes were labeled
  write_csv(label_df, sprintf("figures/%s_secreted_one_to_one_labeled_topgenes.csv", g))
  
  # Plot
  p <- ggplot(df, aes(x = TPM_mouse, y = TPM_human)) +
    geom_point(alpha = 0.8, size = 1.8, color = "#a6c9bd") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_line(data = pred_df, aes(x = x, y = y), color = "grey50", linewidth = 0.7) +  # orange regression
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = label_combined),
      size = 2.4,
      color = "#111111",
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.15,
      min.segment.length = 0
    ) +
    scale_x_continuous(trans = "log1p") +
    scale_y_continuous(trans = "log1p") +
    labs(
      title = sprintf("%s secreted one-to-one orthologs: Mouse vs Human", g),
      subtitle = paste0(
        "Dashed line = y=x (equal expression)\n",
        "Orange line = linear regression (Pearson correlation)\n",
        "r(Pearson raw)=", round(r_pearson_raw, 3),
        " | r(Spearman raw)=", round(r_spearman_raw, 3),
        " | r(Pearson log1p)=", round(r_pearson_log, 3)
      ),
      x = "Mouse Mean TPM (log1p)",
      y = "Human Mean TPM (log1p)"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9)
    )
  
  out_png <- sprintf("figures/%s_secreted_one_to_one_correlation_labeled_human_mouse.png", g)
  ggsave(out_png, p, width = 7.5, height = 6, dpi = 300)
  message(sprintf("[%s] Wrote: %s", g, out_png))
  
  list(
    plot = p,
    r_pearson_raw = r_pearson_raw,
    r_spearman_raw = r_spearman_raw,
    r_pearson_log = r_pearson_log
  )
}

# ---------- Run for all glands ----------
results <- map(glands, function(g) {
  df <- build_secreted_unified(g)
  stats <- plot_corr(g, df)
  tibble(
    gland = g,
    n_pairs = nrow(df),
    r_pearson_raw = stats$r_pearson_raw,
    r_spearman_raw = stats$r_spearman_raw,
    r_pearson_log = stats$r_pearson_log
  )
})

summary_tbl <- bind_rows(results)
write_csv(summary_tbl, "figures/ALL_secreted_one_to_one_correlation_summary.csv")
print(summary_tbl)
