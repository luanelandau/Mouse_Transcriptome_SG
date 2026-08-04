# =========================================================
# Volcano plots for sex differences in human & mouse SM glands
# =========================================================
setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

# ---------- Load data ----------
hum_file <- "deseq2_results_hum/DESeq2_sex_SM.csv"
mouse_file <- "deseq2_results/DESeq2_sex_SM.csv"

hum <- read_csv(hum_file, show_col_types = FALSE)
mouse <- read_csv(mouse_file, show_col_types = FALSE)

# ---------- Prepare data ----------
if (!"Gene" %in% names(hum)) hum <- hum %>% rename(Gene = 1)
if (!"Gene" %in% names(mouse)) mouse <- mouse %>% rename(Gene = 1)

hum <- hum %>% mutate(species = "Human")
mouse <- mouse %>% mutate(species = "Mouse")

# ---------- Parameters ----------
lfc_thr <- 1
padj_thr <- 0.05

# ---------- Count summary ----------
count_sig <- function(df) {
  df %>%
    mutate(
      sig = case_when(
        padj < padj_thr & log2FoldChange > lfc_thr  ~ "Up in males",
        padj < padj_thr & log2FoldChange < -lfc_thr ~ "Up in females",
        TRUE ~ "Not significant"
      )
    ) %>%
    count(sig)
}

cat("\nMouse significant sex-biased genes:\n")
print(count_sig(mouse))

cat("\nHuman significant sex-biased genes:\n")
print(count_sig(hum))

# ---------- Determine shared Y-axis limits ----------
# compute -log10(pvalue) for both to find max
hum$neglogp <- -log10(hum$pvalue)
mouse$neglogp <- -log10(mouse$pvalue)
shared_ymax <- max(hum$neglogp, mouse$neglogp, na.rm = TRUE)


# ---------- Volcano plot function ----------
plot_volcano <- function(df, species_name, lfc_thr = 1, padj_thr = 0.05, shared_ymax = NULL) {
  df <- df %>%
    mutate(
      sig = case_when(
        padj < padj_thr & log2FoldChange >  lfc_thr ~ "Up in males",
        padj < padj_thr & log2FoldChange < -lfc_thr ~ "Up in females",
        TRUE ~ "Not significant"
      ),
      neglogp = -log10(pvalue),
      is_klk = case_when(
        species_name == "Mouse" & grepl("^Klk", Gene, ignore.case = FALSE) ~ TRUE,
        species_name == "Human" & grepl("^KLK", Gene, ignore.case = FALSE) ~ TRUE,
        TRUE ~ FALSE
      ),
      is_sig = padj < padj_thr
    )
  
  base_cols <- c(
    "Up in females"   = "#d391a0",
    "Up in males"     = "#5d7b9f",
    "Not significant" = "grey80"
  )
  
  ggplot(df, aes(x = log2FoldChange, y = neglogp)) +
    # all points
    geom_point(aes(color = sig), alpha = 0.7, size = 1.8, na.rm = TRUE) +
    scale_color_manual(values = base_cols) +
    # ONLY significant KLK/Klk overlay (filled by class, thick black outline)
    geom_point(
      data = df %>% dplyr::filter(is_klk & is_sig),
      aes(fill = sig),
      shape = 21, size = 2.8, color = "black", stroke = 1.2, na.rm = TRUE
    ) +
    scale_fill_manual(values = base_cols, guide = "none") +
    # labels for ONLY significant KLK/Klk
    ggrepel::geom_text_repel(
      data = df %>% dplyr::filter(is_klk & is_sig),
      aes(label = Gene),
      size = 3,
      box.padding = 0.3,
      max.overlaps = Inf,
      segment.color = "grey40",
      force = 6
    ) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    labs(
      title = paste0(species_name, " submandibular gland"),
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~pvalue)
    ) +
    coord_cartesian(ylim = if (is.null(shared_ymax)) NULL else c(0, shared_ymax)) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 14)
    )
}

# ---------- Build plots ----------
volc_hum <- plot_volcano(hum, "Human")
volc_mouse <- plot_volcano(mouse, "Mouse")

combined <- volc_mouse + volc_hum + plot_layout(ncol = 2, guides = "collect")

ggsave("figures/Volcano_SexDifferences_Mouse_Human_SM_Klk_outline.pdf",
       combined, width = 15, height = 5, useDingbats = FALSE)
ggsave("figures/Volcano_SexDifferences_Mouse_Human_SM_Klk_outline.svg",
       combined, width = 15, height = 5)

ggsave("figures/Volcano_SexDifferences_Mouse_Human_SM_Klk_outline.png",
       combined, width = 15, height = 5)

combined









# =========================================================
# Identify and inspect extremely significant genes
# =========================================================

# Threshold in -log10 scale
cutoff <- 50

# Create filtered tables
top_mouse <- mouse %>%
  mutate(neglogp = -log10(pvalue)) %>%
  filter(neglogp > cutoff) %>%
  arrange(desc(neglogp)) %>%
  select(Gene, log2FoldChange, pvalue, padj, neglogp)

top_human <- hum %>%
  mutate(neglogp = -log10(pvalue)) %>%
  filter(neglogp > cutoff) %>%
  arrange(desc(neglogp)) %>%
  select(Gene, log2FoldChange, pvalue, padj, neglogp)

# Print summaries
cat("\nGenes with -log10(pvalue) > 50 in MOUSE:\n")
print(top_mouse)

cat("\nGenes with -log10(pvalue) > 50 in HUMAN:\n")
print(top_human)

# Optional: quick label plot to visualize them
ggplot(top_mouse, aes(x = log2FoldChange, y = neglogp, label = Gene)) +
  geom_point(color = "#5d7b9f", size = 2.5) +
  geom_text_repel(size = 3) +
  theme_classic() +
  labs(title = "Mouse genes with -log10(pvalue) > 50",
       x = expression(Log[2]~Fold~Change),
       y = expression(-Log[10]~pvalue))

###Plot volcano with the top genes. 
#Change label_cutoff for which genes you want labelled.
plot_volcano <- function(df, species_name, lfc_thr = 1, padj_thr = 0.05,
                         shared_ymax = NULL, label_cutoff = 50) {
  df <- df %>%
    mutate(
      sig = case_when(
        padj < padj_thr & log2FoldChange >  lfc_thr ~ "Up in males",
        padj < padj_thr & log2FoldChange < -lfc_thr ~ "Up in females",
        TRUE ~ "Not significant"
      ),
      neglogp = -log10(pvalue),
      is_klk = case_when(
        species_name == "Mouse" & grepl("^Klk", Gene, ignore.case = FALSE) ~ TRUE,
        species_name == "Human" & grepl("^KLK", Gene, ignore.case = FALSE) ~ TRUE,
        TRUE ~ FALSE
      ),
      is_sig = padj < padj_thr
    )
  
  base_cols <- c(
    "Up in females"   = "#d391a0",
    "Up in males"     = "#5d7b9f",
    "Not significant" = "grey80"
  )
  
  # data for labels
  lab_highsig <- df %>% dplyr::filter(neglogp > label_cutoff, !is_klk)   # all very significant, excluding KLKs
  lab_klk     <- df %>% dplyr::filter(is_klk & is_sig)                   # significant KLKs
  
  ggplot(df, aes(x = log2FoldChange, y = neglogp)) +
    # base points
    geom_point(aes(color = sig), alpha = 0.7, size = 1.8, na.rm = TRUE) +
    scale_color_manual(values = base_cols) +
    # KLK overlay: same fill as class, thick black outline, on top
    geom_point(
      data = df %>% dplyr::filter(is_klk & is_sig),
      aes(fill = sig),
      shape = 21, size = 2.8, color = "black", stroke = 1.2, na.rm = TRUE
    ) +
    scale_fill_manual(values = base_cols, guide = "none") +
    # Labels: first non-KLK very significant genes
    ggrepel::geom_text_repel(
      data = lab_highsig,
      aes(label = Gene),
      size = 3,
      box.padding = 0.3,
      max.overlaps = Inf,
      segment.color = "grey25",
      force = 7,
      na.rm = TRUE
    ) +
    # Then KLK/Klk significant labels
    ggrepel::geom_text_repel(
      data = lab_klk,
      aes(label = Gene),
      size = 3,
      box.padding = 0.3,
      max.overlaps = Inf,
      segment.color = "grey25",
      force = 7,
      na.rm = TRUE
    ) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    labs(
      title = paste0(species_name, " submandibular gland"),
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~pvalue)
    ) +
    coord_cartesian(ylim = if (is.null(shared_ymax)) c(0, max(df$neglogp, na.rm = TRUE)) else c(0, shared_ymax)) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 14)
    )
}

# ---------- Build plots ----------
volc_hum <- plot_volcano(hum, "Human")
volc_mouse <- plot_volcano(mouse, "Mouse")

combined <- volc_mouse + volc_hum + plot_layout(ncol = 2, guides = "collect")

ggsave("figures/Volcano_SexDifferences_Mouse_Human_SM_Klk_outline_top_genes.pdf",
       combined, width = 15, height = 5, useDingbats = FALSE)

combined







#March 2026 - volcano without labels 

plot_volcano <- function(df, species_name, lfc_thr = 1, padj_thr = 0.05, shared_ymax = NULL) {
  df <- df %>%
    mutate(
      sig = case_when(
        padj < padj_thr & log2FoldChange >  lfc_thr ~ "Up in males",
        padj < padj_thr & log2FoldChange < -lfc_thr ~ "Up in females",
        TRUE ~ "Not significant"
      ),
      neglogp = -log10(pvalue),
      is_sig = padj < padj_thr
    )
  
  base_cols <- c(
    "Up in females"   = "#d391a0",
    "Up in males"     = "#5d7b9f",
    "Not significant" = "grey80"
  )
  
  ggplot(df, aes(x = log2FoldChange, y = neglogp)) +
    # all points
    geom_point(aes(color = sig), alpha = 0.7, size = 1.8, na.rm = TRUE) +
    scale_color_manual(values = base_cols) +
    scale_fill_manual(values = base_cols, guide = "none") +
    # labels for ONLY significant KLK/Klk
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    labs(
      title = paste0(species_name, " submandibular gland"),
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~pvalue)
    ) +
    coord_cartesian(ylim = if (is.null(shared_ymax)) NULL else c(0, shared_ymax)) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 14)
    )
}

volc_mouse <- plot_volcano(mouse, "Mouse")






















# =========================================================
# Identify and inspect extremely significant genes
# =========================================================

# Threshold in -log10 scale

cutoff <- 50

# Create filtered tables
top_mouse <- mouse %>%
  mutate(neglogp = -log10(pvalue)) %>%
  filter(
    neglogp > 50 |
      (log2FoldChange > 5 & neglogp > 20)
  ) %>%
  arrange(desc(neglogp)) %>%
  select(Gene, log2FoldChange, pvalue, padj, neglogp)

top_human <- hum %>%
  mutate(neglogp = -log10(pvalue)) %>%
  filter(
    neglogp > 20 ) %>%
  arrange(desc(neglogp)) %>%
  select(Gene, log2FoldChange, pvalue, padj, neglogp)

# Print summaries
cat("\nGenes with -log10(pvalue) > 50 in MOUSE:\n")
print(top_mouse)

cat("\nGenes with -log10(pvalue) > 50 in HUMAN:\n")
print(top_human)

###Plot volcano with the top genes. 
#Change label_cutoff for which genes you want labelled.
plot_volcano <- function(df, species_name,
                         lfc_thr = 1,
                         padj_thr = 0.05,
                         shared_ymax = NULL,
                         label_df = NULL) {
  
  df <- df %>%
    mutate(
      sig = case_when(
        padj < padj_thr & log2FoldChange >  lfc_thr ~ "Up in males",
        padj < padj_thr & log2FoldChange < -lfc_thr ~ "Up in females",
        TRUE ~ "Not significant"
      ),
      neglogp = -log10(pvalue)
    )
  
  base_cols <- c(
    "Up in females"   = "#d391a0",
    "Up in males"     = "#5d7b9f",
    "Not significant" = "grey80"
  )
  
  # If a label table is provided, merge to get coordinates
  if (!is.null(label_df)) {
    lab_data <- df %>% 
      dplyr::filter(Gene %in% label_df$Gene)
  } else {
    lab_data <- NULL
  }
  
  ggplot(df, aes(x = log2FoldChange, y = neglogp)) +
    geom_point(aes(color = sig), alpha = 0.7, size = 1.8, na.rm = TRUE) +
    scale_color_manual(values = base_cols) +
    ggrepel::geom_text_repel(
      data = lab_data,
      aes(label = Gene),
      size = 3,
      box.padding = 0.3,
      max.overlaps = Inf,
      segment.color = "grey25",
      na.rm = TRUE
    ) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr),
               linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(padj_thr),
               linetype = "dashed", color = "grey50") +
    labs(
      title = paste0(species_name, " submandibular gland"),
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~pvalue)
    ) +
    coord_cartesian(
      ylim = if (is.null(shared_ymax))
        c(0, max(df$neglogp, na.rm = TRUE))
      else
        c(0, shared_ymax)
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 14)
    )
}

# ---------- Build plots ----------
volc_hum   <- plot_volcano(hum,   "Human", label_df = top_human)
volc_mouse <- plot_volcano(mouse, "Mouse", label_df = top_mouse)

combined <- volc_mouse + volc_hum + plot_layout(ncol = 2, guides = "collect")

combined
ggsave("figures/Volcano_SexDifferences_Mouse_Human_SM_Klk_outline_top_genes.png",
       combined, width = 15, height = 5, dpi=300)#, useDingbats = FALSE)

combined
