library(dplyr)
library(readr)
library(ggplot2)
library(stringr)
library(tidyr)

setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/")


# ------------------------------------------------
# 1. Load permutation files
# ------------------------------------------------
perm_dir <- "cluster_analysis/tmp_perms"
perm_files <- list.files(perm_dir, pattern = "^clustered_\\d+\\.bed$", full.names = TRUE)

read_cluster_file <- function(file, perm_id) {
  df <- read_tsv(file, col_names = FALSE, show_col_types = FALSE)
  colnames(df)[1:5] <- c("chr", "start", "end", "gene", "cluster")
  df %>%
    count(cluster, name = "size") %>%
    count(size, name = "n_clusters") %>%
    mutate(perm = perm_id)
}

perm_df <- purrr::map_dfr(seq_along(perm_files), function(i) {
  if (i %% 50 == 0) message("Processing permutation ", i, "/", length(perm_files))
  read_cluster_file(perm_files[i], i)
})


# ------------------------------------------------
# 2. Load observed file
# ------------------------------------------------
obs_file <- "cluster_analysis/SM_sex_DEGs_clusters_1mb_update.bed"
obs_df <- read_tsv(obs_file, col_names = FALSE, show_col_types = FALSE)
colnames(obs_df)[1:5] <- c("chr", "start", "end", "gene", "cluster")

obs_df <- obs_df %>%
  count(cluster, name = "size") %>%
  count(size, name = "n_clusters")

# ------------------------------------------------
# 3. Calculate 99% percentile of observed n_clusters per size
# ------------------------------------------------
obs_cluster_sizes <- obs_df %>%
  uncount(n_clusters) %>%  # Repeat rows by n_clusters
  pull(size)
obs_size_p99 <- quantile(obs_cluster_sizes, 0.99)
cat("99th percentile of observed cluster sizes:", obs_size_p99, "\n")

# ------------------------------------------------
# 4. Plot
# ------------------------------------------------
p_size_dist <- ggplot(perm_df, aes(x = factor(size), y = n_clusters)) +
  geom_boxplot(outlier.shape = NA, fill = "grey80", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
  geom_point(data = obs_df, aes(x = factor(size), y = n_clusters),
             color = "red", size = 1, inherit.aes = FALSE) +
  geom_vline(xintercept = as.numeric(obs_size_p99), linetype = "dashed", color = "grey30") +
  labs(
    title = "Distribution of #clusters per size (permutations vs observed)",
    subtitle = "Red dots = observed; boxplots = null distribution; blue line = 99%ile of observed",
    x = "Cluster size (# genes)",
    y = "# Clusters"
  ) +
  theme_minimal()

print(p_size_dist)

ggsave("cluster_analysis/cluster_size_distribution.png",
       plot = p_size_dist,
       width = 8, height = 6, dpi = 300)
ggsave("cluster_analysis/cluster_size_distribution.svg",
       plot = p_size_dist,
       width = 8, height = 6)



