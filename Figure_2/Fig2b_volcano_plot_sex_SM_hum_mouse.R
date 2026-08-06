# Volcano plots for sex differences in human and mouse submandibular glands

setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# Load DESeq2 results
hum <- read_csv("deseq2_hum/DESeq2_sex_SM.csv", show_col_types = FALSE)
mouse <- read_csv("deseq2/deseq2_results/DESeq2_sex_SM.csv", show_col_types = FALSE)

if (!"Gene" %in% names(hum)) hum <- rename(hum, Gene = 1)
if (!"Gene" %in% names(mouse)) mouse <- rename(mouse, Gene = 1)

# Genes to label
top_mouse <- mouse %>%
  mutate(neglogp = -log10(pvalue)) %>%
  filter(neglogp > 50 | (log2FoldChange > 5 & neglogp > 20))

top_human <- hum %>%
  mutate(neglogp = -log10(pvalue)) %>%
  filter(neglogp > 20)

plot_volcano <- function(df, species_name, label_df,
                         lfc_thr = 1, padj_thr = 0.05) {
  plot_data <- df %>%
    mutate(
      sig = case_when(
        padj < padj_thr & log2FoldChange > lfc_thr  ~ "Up in males",
        padj < padj_thr & log2FoldChange < -lfc_thr ~ "Up in females",
        TRUE ~ "Not significant"
      ),
      neglogp = -log10(pvalue)
    )
  
  label_data <- filter(plot_data, Gene %in% label_df$Gene)
  
  ggplot(plot_data, aes(log2FoldChange, neglogp)) +
    geom_point(aes(color = sig), alpha = 0.7, size = 1.8, na.rm = TRUE) +
    scale_color_manual(values = c(
      "Up in females" = "#d391a0",
      "Up in males" = "#5d7b9f",
      "Not significant" = "grey80"
    )) +
    geom_text_repel(
      data = label_data,
      aes(label = Gene),
      size = 3,
      box.padding = 0.3,
      max.overlaps = Inf,
      segment.color = "grey25",
      na.rm = TRUE
    ) +
    geom_vline(
      xintercept = c(-lfc_thr, lfc_thr),
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_hline(
      yintercept = -log10(padj_thr),
      linetype = "dashed",
      color = "grey50"
    ) +
    labs(
      title = paste0(species_name, " submandibular gland"),
      x = expression(Log[2] ~ Fold ~ Change),
      y = expression(-Log[10] ~ pvalue)
    ) +
    coord_cartesian(ylim = c(0, max(plot_data$neglogp, na.rm = TRUE))) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 14)
    )
}

volc_hum <- plot_volcano(hum, "Human", label_df = top_human)
volc_mouse <- plot_volcano(mouse, "Mouse", label_df = top_mouse)

combined <- volc_mouse + volc_hum + plot_layout(ncol = 2, guides = "collect")

ggsave(
  "Figure_2/Volcano_SexDifferences_Mouse_Human_SM.png",
  combined,
  width = 15,
  height = 5,
  dpi = 300
)

combined
