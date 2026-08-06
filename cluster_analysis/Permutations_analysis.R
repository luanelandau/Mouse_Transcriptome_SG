#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
})

# ===========================
# Paths (adjust if needed)
# ===========================
setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/cluster_analysis/")

# Summaries generated from the 100-kb clustered permutation and observed BEDs.
perm_tsv_100kb <- "size_distribution_100kb.tsv"                  # header: perm size n_clusters
obs_tsv_100kb  <- "SM_DEGs_observed_size_distribution_100kb.tsv" # headerless: size n_clusters

if (!file.exists(perm_tsv_100kb)) {
  stop("Permutation summary not found: ", perm_tsv_100kb)
}
if (!file.exists(obs_tsv_100kb)) {
  stop("Observed summary not found: ", obs_tsv_100kb)
}

out_dir <- "100kb_analysis"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ===========================
# 1) Load data
# ===========================
# Permutations: with header
perm_df <- read_tsv(
  perm_tsv_100kb,
  col_types = cols(
    perm = col_character(),
    size = col_integer(),
    n_clusters = col_integer()
  )
)

# Observed: no header
obs_df <- read_tsv(
  obs_tsv_100kb,
  col_names = c("size","n_clusters"),
  col_types = cols(size = col_integer(), n_clusters = col_integer())
)

# Complete missing (perm,size) with zeros; same for observed sizes
all_sizes <- sort(unique(c(perm_df$size, obs_df$size)))
all_perms <- sort(unique(perm_df$perm))

perm_df_complete <- perm_df %>%
  tidyr::complete(perm = all_perms, size = all_sizes, fill = list(n_clusters = 0L))

obs_df <- obs_df %>%
  tidyr::complete(size = all_sizes, fill = list(n_clusters = 0L)) %>%
  arrange(size)

n_perm <- length(all_perms)

# ===========================
# 2) Per-size empirical p-values (one-sided: observed >= null)
# ===========================
per_size_obs <- obs_df %>% rename(obs_n_clusters = n_clusters)

per_size_join <- perm_df_complete %>%
  left_join(per_size_obs, by = "size") %>%
  group_by(size) %>%
  summarise(
    obs_n_clusters = first(obs_n_clusters),
    null_ge_obs    = sum(n_clusters >= first(obs_n_clusters)),
    n_perm         = n(),
    p_empirical    = (null_ge_obs + 1) / (n_perm + 1),  # +1 pseudocount
    null_mean      = mean(n_clusters),
    null_sd        = sd(n_clusters),
    .groups = "drop"
  ) %>%
  arrange(size) %>%
  mutate(p_adj_fdr = p.adjust(p_empirical, method = "fdr"))

# ===========================
# 3) Key tests: size==6 and size>=6
# ===========================
# Exactly size == 6
obs_6 <- per_size_join %>% filter(size == 6) %>% pull(obs_n_clusters)
if (length(obs_6) == 0) obs_6 <- 0L
null_6 <- perm_df_complete %>%
  filter(size == 6) %>%
  group_by(perm) %>%
  summarise(n_clusters = sum(n_clusters), .groups = "drop")
if (nrow(null_6) == 0) { null_6 <- tibble(n_clusters = integer(0)) }
p_6 <- if (nrow(null_6) > 0) (sum(null_6$n_clusters >= obs_6) + 1) / (nrow(null_6) + 1) else 1

# Cumulative size >= 6
obs_ge6 <- obs_df %>% filter(size >= 6) %>% summarise(total = sum(n_clusters)) %>% pull(total)
null_ge6 <- perm_df_complete %>%
  filter(size >= 6) %>%
  group_by(perm) %>%
  summarise(total = sum(n_clusters), .groups = "drop")
p_ge6 <- (sum(null_ge6$total >= obs_ge6) + 1) / (nrow(null_ge6) + 1)

# ===========================
# 4) Save tables
# ===========================
summary_tbl <- per_size_join %>%
  transmute(
    size,
    observed_n_clusters = obs_n_clusters,
    null_mean, null_sd,
    p_empirical, p_adj_fdr
  )
write_csv(summary_tbl, file.path(out_dir, "cluster_size_empirical_pvalues_100kb.csv"))

readr::write_lines(
  c(
    sprintf("size==6: observed=%d, p_empirical=%.6f", obs_6, p_6),
    sprintf("size>=6: observed=%d, p_empirical=%.6f", obs_ge6, p_ge6),
    sprintf("n_permutations=%d", n_perm)
  ),
  file.path(out_dir, "cluster_size_key_tests_100kb.txt")
)

# Optional: permutation cutoffs (95%/99%) per size
cutoffs_tbl <- perm_df_complete %>%
  group_by(size) %>%
  summarise(
    q95 = as.integer(quantile(n_clusters, 0.95, type = 7)),
    q99 = as.integer(quantile(n_clusters, 0.99, type = 7)),
    .groups = "drop"
  )
write_csv(cutoffs_tbl, file.path(out_dir, "cluster_size_permutation_cutoffs_100kb.csv"))

# ===========================
# 5) Plot
# ===========================
sig_sizes <- summary_tbl %>%
  filter(p_adj_fdr < 0.05) %>%
  pull(size)

p <- ggplot(perm_df_complete, aes(x = factor(size), y = n_clusters)) +
  geom_boxplot(outlier.shape = NA, fill = "grey85", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.08, size = 0.5) +
  geom_point(
    data = obs_df,
    aes(x = factor(size), y = n_clusters),
    color = "red", size = 1.6, inherit.aes = FALSE
  ) +
  geom_point(
    data = obs_df %>% dplyr::filter(size %in% sig_sizes),
    aes(x = factor(size), y = n_clusters),
    shape = 21, stroke = 0.5, size = 3.2, fill = NA, color = "red",
    inherit.aes = FALSE
  ) +
  labs(
    title = "Clusters per size (100 kb): observed vs permutation null",
    subtitle = "Red = observed; boxes/jitter = permutations. Open circles: FDR < 0.05.",
    x = "Cluster size (# genes)",
    y = "# clusters"
  ) +
  theme_minimal(base_size = 12)
p

ggsave(file.path(out_dir, "cluster_size_distribution_100kb.png"), p, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "cluster_size_distribution_100kb.svg"), p, width = 8, height = 6)

# ===========================
# 6) Expected-versus-observed line plot
# ===========================
# Include every integer cluster size so that absent sizes are shown as zero
# instead of connecting lines across missing values.
plot_sizes <- seq.int(min(all_sizes), max(all_sizes))

expected_line <- perm_df_complete %>%
  tidyr::complete(perm = all_perms, size = plot_sizes, fill = list(n_clusters = 0L)) %>%
  group_by(size) %>%
  summarise(
    expected_mean = mean(n_clusters),
    expected_sd = sd(n_clusters),
    .groups = "drop"
  ) %>%
  mutate(
    expected_lower = pmax(expected_mean - expected_sd, 0),
    expected_upper = expected_mean + expected_sd
  )

observed_line <- obs_df %>%
  tidyr::complete(size = plot_sizes, fill = list(n_clusters = 0L)) %>%
  arrange(size)

p_line <- ggplot() +
  geom_ribbon(
    data = expected_line,
    aes(x = size, ymin = expected_lower, ymax = expected_upper, fill = "Expected +/- 1 SD"),
    alpha = 0.22
  ) +
  geom_line(
    data = expected_line,
    aes(x = size, y = expected_mean, color = "Expected"),
    linewidth = 1
  ) +
  geom_point(
    data = expected_line,
    aes(x = size, y = expected_mean, color = "Expected"),
    size = 1.5
  ) +
  geom_line(
    data = observed_line,
    aes(x = size, y = n_clusters, color = "Observed"),
    linewidth = 1
  ) +
  geom_point(
    data = observed_line,
    aes(x = size, y = n_clusters, color = "Observed"),
    size = 1.8
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Expected" = "#3366AA", "Observed" = "#CC3311")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("Expected +/- 1 SD" = "#3366AA")
  ) +
  scale_x_continuous(breaks = plot_sizes) +
  scale_y_continuous(trans = "log1p") +
  labs(
    title = "Observed versus expected clusters by size (100 kb)",
    subtitle = "Expected = mean across permutations; shading = +/- 1 SD",
    x = "Cluster size (# genes)",
    y = "# clusters (log1p scale)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

print(p_line)

ggsave(
  file.path(out_dir, "cluster_size_expected_vs_observed_100kb.png"),
  p_line, width = 9, height = 6, dpi = 300
)
ggsave(
  file.path(out_dir, "cluster_size_expected_vs_observed_100kb.svg"),
  p_line, width = 9, height = 6
)

# Console summary
cat("\n--- 100 kb key tests ---\n")
cat(sprintf("Exactly size==6: observed=%d | empirical p=%.6f\n", obs_6, p_6))
cat(sprintf("Cumulative size>=6: observed=%d | empirical p=%.6f\n", obs_ge6, p_ge6))
cat(sprintf("Per-size results → %s\n", file.path(out_dir, "cluster_size_empirical_pvalues_100kb.csv")))
