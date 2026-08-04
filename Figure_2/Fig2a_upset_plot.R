## =========================================================
## UpSet: ALL sex comparisons (PAR, SM, SL, PANC, LIV)
## =========================================================
setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ComplexUpset)
  library(tidyr)
  library(purrr)
})

#dir.create("figures", showWarnings = FALSE, recursive = TRUE)

## Params
l2fc_min <- 1
padj_cut <- 0.05

files <- c(
  PAR  = "deseq2/deseq2_results/DESeq2_sex_PAR.csv",
  SM   = "deseq2/deseq2_results/DESeq2_sex_SM.csv",
  SL   = "deseq2/deseq2_results/DESeq2_sex_SL.csv",
  PANC = "deseq2/deseq2_results/DESeq2_sex_PANC.csv",
  LIV  = "deseq2/deseq2_results/DESeq2_sex_LIV.csv"
)
glands <- names(files)

## Helpers
read_deseq_table <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE, col_names = TRUE)
  if (is.na(names(df)[1]) || names(df)[1] == "" || names(df)[1] %in% c("...1","X1","X")) {
    names(df)[1] <- "mouse_gene"
  }
  if (!"mouse_gene" %in% names(df)) names(df)[1] <- "mouse_gene"
  df %>% mutate(mouse_gene = as.character(mouse_gene)) %>% distinct(mouse_gene, .keep_all = TRUE)
}
sig_genes <- function(df, l2fc_min = 1, padj_cut = 0.05) {
  df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    filter(padj < padj_cut & abs(log2FoldChange) > l2fc_min) %>%
    pull(mouse_gene) %>% unique()
}

## Read all sex comparisons and stack
sets <- imap(files, ~ {
  tib <- read_deseq_table(.x)
  genes <- sig_genes(tib, l2fc_min, padj_cut)
  tibble(mouse_gene = genes, set = .y, present = TRUE)
})

## Build membership by widening the stacked table
membership <- bind_rows(sets) %>%
  mutate(set = factor(set, levels = glands)) %>%
  tidyr::pivot_wider(
    names_from  = set,
    values_from = present,
    values_fill = FALSE
  )

## Ensure all gland columns exist (in case a set has 0 DEGs)
for (g in glands) {
  if (!g %in% names(membership)) membership[[g]] <- FALSE
}

## Keep rows present in at least one set
up_df <- membership %>%
  filter(rowSums(across(all_of(glands))) > 0)

## Quick sanity print
msg_counts <- colSums(select(up_df, all_of(glands)))
print(msg_counts)

## Plot
p_up <- ComplexUpset::upset(
  up_df,
  intersect = glands,
  min_size=5,
  name = "Sex DEGs",
  base_annotations = list("Intersection size" = intersection_size(text = list(size = 3))),
  set_sizes = upset_set_size(),
  encode_sets = FALSE,
  width_ratio = 0.12
) + theme(text = element_text(size = 11))

print(p_up)

## Save
ggsave("Figure_2/Fig2a_ALLglands_upset_SEX.png", p_up, width = 3, height = 6, dpi = 300)
ggsave("Figure_2/Fig2a_ALLglands_upset_SEX.svg", p_up, width = 3, height = 6)


## ============================================
## Get the genes in each (exclusive) intersection
## ============================================
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# all non-empty combinations of glands (PAR, SM, SL, PANC, LIV)
all_combos <- unlist(
  lapply(seq_along(glands), function(k) combn(glands, k, simplify = FALSE)),
  recursive = FALSE
)

# helper: genes that are TRUE in exactly S (and FALSE in others)
genes_in_combo <- function(S) {
  in_S  <- rowSums(up_df[, S, drop = FALSE]) == length(S)
  others <- setdiff(glands, S)
  out_S <- if (length(others) == 0) TRUE else rowSums(up_df[, others, drop = FALSE]) == 0
  up_df$mouse_gene[in_S & out_S]
}

# compute list of genes per intersection
intersection_list <- lapply(all_combos, genes_in_combo)
names(intersection_list) <- sapply(all_combos, paste, collapse = "&")

# tidy table: one row per intersection, with a list-column of genes
intersection_tbl <- tibble::tibble(
  intersection = names(intersection_list),
  degree       = lengths(strsplit(names(intersection_list), "&")),
  genes        = intersection_list
) %>%
  mutate(n_genes = lengths(genes)) %>%
  arrange(desc(n_genes), degree, intersection)

# (Optional) match your plot filter: keep intersections with at least 5 genes
intersection_tbl_plot <- intersection_tbl %>% filter(n_genes >= 5)

# View in RStudio
View(intersection_tbl_plot)   # or View(intersection_tbl) for all

# Export a long CSV: one gene per row per intersection (nice for Supplementary)
intersection_long <- intersection_tbl_plot %>%
  tidyr::unnest_longer(genes, values_to = "mouse_gene") %>%
  select(intersection, degree, n_genes, mouse_gene)

#readr::write_csv(intersection_long, "figures/ALLglands_upset_SEX_intersections_genes_min5.csv")

# If you also want a compact CSV (genes joined in one cell):
intersection_compact <- intersection_tbl_plot %>%
  mutate(genes_joined = vapply(genes, function(v) paste(v, collapse = ", "), character(1))) %>%
  select(intersection, degree, n_genes, genes_joined)

#readr::write_csv(intersection_compact, "figures/ALLglands_upset_SEX_intersections_genes_min5_compact.csv")

