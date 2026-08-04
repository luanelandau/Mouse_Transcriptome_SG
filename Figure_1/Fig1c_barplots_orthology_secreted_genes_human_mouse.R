# ============================
# % of total expression per orthology bucket (secreted, TPM >= 2)
# - Loads all human + mouse mastersheets
# - Buckets orthology: one-to-one / lineage-specific / other
# - Computes % of total TPM per bucket
# - Saves one bar plot per gland (PAR, SM, SL), faceted by species
# ============================

setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(stringr)
  library(patchwork)
})

#dir.create("figures", showWarnings = FALSE, recursive = TRUE)
#dir.create("derived", showWarnings = FALSE, recursive = TRUE)

# ---------- Palette (as requested) ----------
orthology3_colors <- c(
  "one-to-one"       = "#a6c9bd",  # greenish
  "lineage-specific" = "#934e65",  # plum/brown
  "other"            = "grey80"    # grey
)

# ---------- Helpers ----------
# Normalize raw orthology labels coming from sheets
normalize_orth <- function(x) {
  lab <- tolower(gsub("[^a-z0-9]+", "_", x))
  dplyr::case_when(
    lab %in% c("one_to_one","1_to_1","one_2_one","one-one") ~ "one-to-one",
    lab %in% c("human_specific","human-only","human_only","no_mouse_ortholog") ~ "human-specific",
    lab %in% c("mouse_specific","mouse-only","mouse_only","no_human_ortholog") ~ "mouse-specific",
    lab %in% c("one_to_many","many_to_one","many_to_many","one_2_many","many_2_one","many_2_many") ~ "other",
    lab %in% c("no_ortholog","none","no_match") ~ "other",
    TRUE ~ lab
  )
}

# Collapse to 3 buckets, species-aware
to_three_buckets <- function(fine_label, species) {
  if (species == "Human") {
    dplyr::case_when(
      fine_label == "one-to-one"     ~ "one-to-one",
      fine_label == "human-specific" ~ "lineage-specific",
      TRUE                           ~ "other"
    )
  } else {
    dplyr::case_when(
      fine_label == "one-to-one"     ~ "one-to-one",
      fine_label == "mouse-specific" ~ "lineage-specific",
      TRUE                           ~ "other"
    )
  }
}

# Generic processor for a single CSV (one gland × one species)
# secreted_col: "Secretome location" (human) or "secreted" (mouse)
# tpm_col: "Mean_TPM"
# orth_col: "ortholog_type"
process_sheet <- function(path, species, gland, secreted_col, tpm_col, orth_col) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    # Secreted filter (exclude NA and explicit "No_annotation")
    filter(!is.na(.data[[secreted_col]]), .data[[secreted_col]] != "No_annotation") %>%
    mutate(
      TPM       = suppressWarnings(as.numeric(.data[[tpm_col]])),
      orth_fine = normalize_orth(.data[[orth_col]]),
      bucket    = to_three_buckets(orth_fine, species)
    ) %>%
    filter(!is.na(TPM), TPM >= 2) %>%
    summarize(total_TPM = sum(TPM, na.rm = TRUE), .by = bucket) %>%
    mutate(
      percent = 100 * total_TPM / sum(total_TPM),
      species = species,
      gland   = gland
    ) %>%
    select(gland, species, bucket, percent)
}

# ---------- File paths ----------
hum_dir   <- file.path("miscelaneous_sheets","human expression")
mouse_dir <- file.path("miscelaneous_sheets","mouse_expression")

# Human (your exact files)
PAR_h_path <- file.path(hum_dir, "PAR_human_mastersheet_TPMs_annotated_with_orthology.csv")
SM_h_path  <- file.path(hum_dir, "SM_human_mastersheet_TPMs_annotated_with_orthology.csv")
SL_h_path  <- file.path(hum_dir, "SL_human_mastersheet_TPMs_annotated_with_orthology.csv")

# Mouse (your mouse mastersheets)
PAR_m_path <- file.path(mouse_dir, "PAR_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv")
SM_m_path  <- file.path(mouse_dir, "SM_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv")
SL_m_path  <- file.path(mouse_dir, "SL_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv")

# ---------- Build tidy summary for each gland × species ----------
# Human: secreted = "Secretome location", TPM = "Mean_TPM", orthology = "ortholog_type"
PAR_h_sum <- process_sheet(PAR_h_path, "Human", "PAR", "Secretome location", "Mean_TPM", "ortholog_type")
SM_h_sum  <- process_sheet(SM_h_path,  "Human", "SM",  "Secretome location", "Mean_TPM", "ortholog_type")
SL_h_sum  <- process_sheet(SL_h_path,  "Human", "SL",  "Secretome location", "Mean_TPM", "ortholog_type")

# Mouse: secreted = "secreted", TPM = "Mean_TPM", orthology = "ortholog_type"
PAR_m_sum <- process_sheet(PAR_m_path, "Mouse", "PAR", "secreted", "Mean_TPM", "ortholog_type")
SM_m_sum  <- process_sheet(SM_m_path,  "Mouse", "SM",  "secreted", "Mean_TPM", "ortholog_type")
SL_m_sum  <- process_sheet(SL_m_path,  "Mouse", "SL",  "secreted", "Mean_TPM", "ortholog_type")

orth_pct <- bind_rows(PAR_h_sum, SM_h_sum, SL_h_sum,
                      PAR_m_sum, SM_m_sum, SL_m_sum) %>%
  mutate(
    bucket  = factor(bucket,  levels = c("one-to-one","lineage-specific","other")),
    species = factor(species, levels = c("Human","Mouse")),
    gland   = factor(gland,   levels = c("PAR","SM","SL"))
  )

# Save summary table
#readr::write_csv(orth_pct, file.path("derived", "percent_total_expression_by_orthology_secreted_TPM10.csv"))

# ---------- Plotter: one figure per gland; facets = species ----------
plot_pct_by_gland <- function(df, gland_code, title_map = c(PAR="Parotid", SM="Submandibular", SL="Sublingual")) {
  ggplot(df %>% filter(gland == gland_code),
         aes(x = bucket, y = percent, fill = bucket)) +
    geom_col(width = 0.7, color = NA) +
    facet_wrap(~ species, nrow = 1) +
    scale_fill_manual(values = orthology3_colors, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = paste0(title_map[[gland_code]], ": % total expression by orthology (secreted, TPM \u2265 10)"),
      x = NULL, y = "Percent of total TPM"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(face = "bold"),
      strip.text  = element_text(face = "bold")
    )
}

# ---------- Make & save figures ----------
p_par_pct <- plot_pct_by_gland(orth_pct, "PAR")
p_sm_pct  <- plot_pct_by_gland(orth_pct, "SM")
p_sl_pct  <- plot_pct_by_gland(orth_pct, "SL")

#ggsave("figures/PAR_percent_total_expression_orthology_TPM10_secreted.svg", p_par_pct, width = 10, height = 5, dpi = 300)
#ggsave("figures/SM_percent_total_expression_orthology_TPM10_secreted.svg",  p_sm_pct,  width = 10, height = 5, dpi = 300)
#ggsave("figures/SL_percent_total_expression_orthology_TPM10_secreted.svg",  p_sl_pct,  width = 10, height = 5, dpi = 300)

# Optional combined grid (vertical stack)
#combined_pct <- p_par_pct + p_sm_pct + p_sl_pct
#ggsave("Figure_1/percent_total_expression_orthology_TPM10_secreted_combined.png", combined_pct, width = 10, height = 6, dpi = 300)
#ggsave("Figure_1/percent_total_expression_orthology_TPM10_secreted_combined.svg", combined_pct, width = 10, height = 6, dpi = 300)

# Done.



# ---------- Shared y-limit across all panels ----------
global_ymax <- ceiling(max(orth_pct$percent, na.rm = TRUE) / 10) * 10
if (!is.finite(global_ymax) || global_ymax <= 0) global_ymax <- 100

# ---------- plotter: one figure per gland; facets = species; shared y via global_ymax ----------
plot_pct_by_gland <- function(df, gland_code, title_map = c(PAR="Parotid", SM="Submandibular", SL="Sublingual"),
                              show_y_title = TRUE) {
  p <- ggplot(df %>% filter(gland == gland_code),
              aes(x = bucket, y = percent, fill = bucket)) +
    geom_col(width = 0.7, color = NA) +
    facet_wrap(~ species, nrow = 1) +
    scale_fill_manual(values = orthology3_colors, drop = FALSE) +
    scale_y_continuous(limits = c(0, global_ymax), expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = paste0(title_map[[gland_code]], ": % total expression by orthology (secreted, TPM \u2265 10)"),
      x = NULL, y = if (show_y_title) "Percent of total TPM" else NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(face = "bold"),
      strip.text  = element_text(face = "bold"),
      # hide y ticks/text when not showing y title (for center/right panels)
      axis.text.y = if (show_y_title) element_text() else element_blank(),
      axis.ticks.y = if (show_y_title) element_line() else element_blank()
    )
  p
}

# ---------- Make figures (left shows y-axis label; middle/right hide it) ----------
p_par_pct <- plot_pct_by_gland(orth_pct, "PAR", show_y_title = TRUE)
p_sm_pct  <- plot_pct_by_gland(orth_pct, "SM",  show_y_title = FALSE)
p_sl_pct  <- plot_pct_by_gland(orth_pct, "SL",  show_y_title = FALSE)

# ---------- Combined: side-by-side, matched y, single y label on the left ----------
combined_pct_h <- p_par_pct + p_sm_pct + p_sl_pct + patchwork::plot_layout(nrow = 1)

ggsave("Figure_1/percent_total_expression_orthology_TPM10_secreted_combined_horizontal.png",
       combined_pct_h, width = 10, height = 6, dpi = 300)
ggsave("Figure_1/percent_total_expression_orthology_TPM10_secreted_combined_horizontal.svg",
       combined_pct_h, width = 10, height = 6)

