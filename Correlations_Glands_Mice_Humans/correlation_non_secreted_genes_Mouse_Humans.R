# =========================================================
# Non-secreted one-to-one orthologs: build tables + correlations for PAR, SM, SL
# - Writes:  figures/<GLAND>_NONsecreted_one_to_one_unified.csv
# - Plots:   figures/<GLAND>_NONsecreted_one_to_one_correlation_labeled_human_mouse.png
# - Summary: figures/ALL_NONsecreted_one_to_one_correlation_summary.csv
#
# NOTE (mouse filter): we keep strictly "Non_secreted" and EXCLUDE "No_annotation".
# To include unannotated too, replace the mouse filter with: secreted %in% c("Non_secreted","No_annotation")
# NOTE (human filter): we keep genes whose Secretome location is NA or does NOT contain "Secreted".
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
# ---------- Helpers (robust non-secreted classification) ----------
human_path <- function(g) sprintf("miscelaneous_sheets/human expression/%s_human_mastersheet_TPMs_annotated_with_orthology.csv", g)
mouse_path <- function(g) sprintf("miscelaneous_sheets/mouse_expression/%s_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv", g)

# Human: NON-secreted if Secretome location is NA OR does NOT contain "Secreted"
is_human_nonsecreted <- function(x) {
  is.na(x) | !grepl("Secreted", x, ignore.case = TRUE)
}

# Mouse: normalize and classify non-secreted
is_mouse_nonsecreted <- function(x) {
  x2 <- tolower(trimws(as.character(x)))
  # include explicit non-secreted and unknowns as "non-secreted" for this analysis
  is.na(x2) | x2 %in% c("non_secreted","non-secreted","no_annotation","no annotation","na","unknown")
}

# Build unified NON-secreted 1:1 for a gland (includes unknowns)
build_nonsecreted_unified <- function(g) {
  human_df <- readr::read_csv(human_path(g), show_col_types = FALSE)
  mouse_df <- readr::read_csv(mouse_path(g), show_col_types = FALSE)
  
  human_1to1_nonsec <- human_df %>%
    dplyr::filter(ortholog_type == "one-to-one") %>%
    dplyr::filter(is_human_nonsecreted(`Secretome location`)) %>%
    dplyr::transmute(
      human_gene          = Geneid,
      mouse_gene_expected = mouse_gene,
      TPM_human           = Mean_TPM
    )
  
  mouse_1to1_nonsec <- mouse_df %>%
    dplyr::filter(ortholog_type == "one-to-one") %>%
    dplyr::filter(is_mouse_nonsecreted(secreted)) %>%
    dplyr::transmute(
      human_gene = human_gene,
      mouse_gene = Geneid,
      TPM_mouse  = Mean_TPM
    )
  
  unified <- human_1to1_nonsec %>%
    dplyr::inner_join(mouse_1to1_nonsec, by = "human_gene") %>%
    dplyr::filter(is.na(mouse_gene_expected) | mouse_gene_expected == mouse_gene) %>%
    dplyr::select(human_gene, mouse_gene, TPM_human, TPM_mouse) %>%
    dplyr::group_by(human_gene, mouse_gene) %>%
    dplyr::summarise(
      TPM_human = sum(TPM_human, na.rm = TRUE),
      TPM_mouse = sum(TPM_mouse, na.rm = TRUE),
      .groups = "drop"
    )
  
  out_csv <- sprintf("figures/%s_NONsecreted_one_to_one_unified.csv", g)
  readr::write_csv(unified, out_csv)
  message(sprintf("[%s] Wrote: %s (n=%d)", g, out_csv, nrow(unified)))
  unified
}

# Plot correlation with safety for small n
plot_corr_nonsec <- function(g, df, TOP_N_PER_SPECIES = 20) {
  if (nrow(df) < 2) {
    message(sprintf("[%s] Fewer than 2 pairs; skipping regression. Writing scatter only.", g))
    p <- ggplot2::ggplot(df, ggplot2::aes(TPM_mouse, TPM_human)) +
      ggplot2::geom_point(alpha = 0.9, size = 2, color = "grey35") +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
      ggplot2::scale_x_continuous(trans = "log1p") +
      ggplot2::scale_y_continuous(trans = "log1p") +
      ggplot2::labs(
        title = sprintf("%s non-secreted one-to-one orthologs: Mouse vs Human", g),
        subtitle = "Insufficient pairs for regression",
        x = "Mouse Mean TPM (log1p)", y = "Human Mean TPM (log1p)"
      ) + ggplot2::theme_classic(base_size = 10)
    out_png <- sprintf("figures/%s_NONsecreted_one_to_one_correlation_labeled_human_mouse.png", g)
    ggplot2::ggsave(out_png, p, width = 7.5, height = 6, dpi = 300)
    return(list(plot = p, n_pairs = nrow(df), r_pearson_raw = NA, r_spearman_raw = NA, r_pearson_log = NA))
  }
  
  r_pearson_raw  <- suppressWarnings(cor(df$TPM_mouse, df$TPM_human,  method = "pearson"))
  r_spearman_raw <- suppressWarnings(cor(df$TPM_mouse, df$TPM_human,  method = "spearman"))
  r_pearson_log  <- suppressWarnings(cor(log1p(df$TPM_mouse), log1p(df$TPM_human), method = "pearson"))
  
  fit <- lm(log1p(TPM_human) ~ log1p(TPM_mouse), data = df)
  coef_fit <- coef(fit)
  xlog_seq <- seq(min(log1p(df$TPM_mouse), na.rm = TRUE),
                  max(log1p(df$TPM_mouse), na.rm = TRUE), length.out = 300)
  ylog_hat <- coef_fit[1] + coef_fit[2] * xlog_seq
  pred_df  <- data.frame(x = pmax(expm1(xlog_seq), 0), y = pmax(expm1(ylog_hat), 0))
  
  top_human <- df %>% dplyr::slice_max(order_by = TPM_human, n = TOP_N_PER_SPECIES) %>% dplyr::mutate(label_source = "human_top")
  top_mouse <- df %>% dplyr::slice_max(order_by = TPM_mouse, n = TOP_N_PER_SPECIES) %>% dplyr::mutate(label_source = "mouse_top")
  label_df  <- dplyr::bind_rows(top_human, top_mouse) %>%
    dplyr::distinct(human_gene, .keep_all = TRUE) %>%
    dplyr::mutate(label_combined = paste0(human_gene, " | ", mouse_gene))
  readr::write_csv(label_df, sprintf("figures/%s_NONsecreted_one_to_one_labeled_topgenes.csv", g))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = TPM_mouse, y = TPM_human)) +
    ggplot2::geom_point(alpha = 0.85, size = 1.8, color = "#a6c9bd") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
    ggplot2::geom_line(data = pred_df, ggplot2::aes(x = x, y = y), color = "grey15", linewidth = 0.7) +
    ggrepel::geom_text_repel(
      data = label_df, ggplot2::aes(label = label_combined),
      size = 2.4, color = "#111111",
      max.overlaps = Inf, box.padding = 0.35, point.padding = 0.15, min.segment.length = 0
    ) +
    ggplot2::scale_x_continuous(trans = "log1p") +
    ggplot2::scale_y_continuous(trans = "log1p") +
    ggplot2::labs(
      title = sprintf("%s non-secreted one-to-one orthologs: Mouse vs Human", g),
      subtitle = paste0(
        "Dashed line = y=x (equal expression)\n",
        "Dark line = linear regression (Pearson correlation)\n",
        "r(Pearson raw)=", round(r_pearson_raw, 3),
        " | r(Spearman raw)=", round(r_spearman_raw, 3),
        " | r(Pearson log1p)=", round(r_pearson_log, 3)
      ),
      x = "Mouse Mean TPM (log1p)", y = "Human Mean TPM (log1p)"
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 9))
  
  out_png <- sprintf("figures/%s_NONsecreted_one_to_one_correlation_labeled_human_mouse.png", g)
  ggplot2::ggsave(out_png, p, width = 7.5, height = 6, dpi = 300)
  message(sprintf("[%s] Wrote: %s", g, out_png))
  
  list(plot = p, n_pairs = nrow(df),
       r_pearson_raw = r_pearson_raw, r_spearman_raw = r_spearman_raw, r_pearson_log = r_pearson_log)
}

results <- purrr::map(glands, function(g) {
  df <- build_nonsecreted_unified(g)
  stats <- plot_corr_nonsec(g, df)
  tibble::tibble(
    gland = g,
    n_pairs = stats$n_pairs,
    r_pearson_raw = stats$r_pearson_raw,
    r_spearman_raw = stats$r_spearman_raw,
    r_pearson_log = stats$r_pearson_log
  )
})

summary_tbl <- dplyr::bind_rows(results)
readr::write_csv(summary_tbl, "figures/ALL_NONsecreted_one_to_one_correlation_summary.csv")
print(summary_tbl)
