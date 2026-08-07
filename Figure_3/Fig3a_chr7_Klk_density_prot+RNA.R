suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})

setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/")

f_prot <- "Figure_2/Saliva_mouse_sex_diff.csv"
f_rna <- "deseq2/deseq2_results/DESeq2_sex_SM.csv"
f_map <- "miscelaneous_sheets/gene_expression_matrix_C57_CD1.csv"

map_clean <- read.csv(f_map, check.names = FALSE) %>%
  transmute(
    gene = Geneid,
    Chr = str_split(Chr, ";", simplify = TRUE)[, 1],
    Start = as.numeric(str_split(Start, ";", simplify = TRUE)[, 1]),
    End = as.numeric(str_split(End, ";", simplify = TRUE)[, 1]),
    pos_mb = floor((Start + End) / 2) / 1e6
  ) %>%
  filter(str_detect(Chr, "^NC_000073\\.7"), is.finite(pos_mb))

rna <- read.csv(f_rna, check.names = FALSE)
colnames(rna)[1] <- "gene"

rna_chr <- rna %>%
  transmute(gene, log2FC = log2FoldChange, padj) %>%
  filter(!is.na(padj), padj < 0.05, is.finite(log2FC), abs(log2FC) >= 1) %>%
  inner_join(map_clean, by = "gene") %>%
  mutate(datatype = "Transcriptome (SM DEGs)")

prot_chr <- read.csv(f_prot, check.names = FALSE) %>%
  transmute(
    gene = as.character(Genes),
    log2FC = -suppressWarnings(as.numeric(`AVG Log2 Ratio`)),
    Qvalue = suppressWarnings(as.numeric(Qvalue))
  ) %>%
  filter(
    !is.na(gene), gene != "", gene != "NaN", !str_detect(gene, ";"),
    is.finite(log2FC), !is.na(Qvalue), Qvalue < 0.05, abs(log2FC) >= 1
  ) %>%
  inner_join(map_clean, by = "gene") %>%
  mutate(datatype = "Proteome (signif. proteins)")

raw_df <- bind_rows(rna_chr, prot_chr)

overlap_genes <- raw_df %>%
  distinct(gene, datatype) %>%
  count(gene) %>%
  filter(n >= 2) %>%
  pull(gene)

raw_df <- raw_df %>%
  mutate(
    bias = if_else(log2FC > 0, "Up in males", "Up in females"),
    overlap = gene %in% overlap_genes,
    color_group = if_else(overlap, bias, "Not significant")
  )

label_df <- filter(raw_df, overlap)
ymax <- ceiling(max(abs(raw_df$log2FC), na.rm = TRUE) * 1.05)

color_map <- c(
  "Up in females" = "#d391a0",
  "Up in males" = "#5d7b9f",
  "Not significant" = "grey80"
)

p_points <- ggplot(raw_df, aes(pos_mb, log2FC)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(aes(color = color_group), alpha = 0.9, size = 1.2) +
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = 2.2,
    max.overlaps = 20,
    box.padding = 0.25,
    point.padding = 0.15,
    segment.size = 0.2,
    segment.alpha = 0.4
  ) +
  scale_color_manual(values = color_map, name = "Overlap bias") +
  facet_wrap(~datatype, ncol = 1) +
  coord_cartesian(ylim = c(-ymax, ymax), expand = FALSE) +
  labs(x = NULL, y = "log2 fold-change (Male/Female)") +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold")
  )

p_density <- ggplot(raw_df, aes(pos_mb)) +
  geom_density(adjust = 0.1, linewidth = 0.8) +
  facet_wrap(~datatype, ncol = 1, scales = "free_y") +
  labs(x = "chr7 position (Mb)", y = "Spatial density\n(of significant genes)") +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold")
  )

p_two_panels <- p_points / p_density + plot_layout(heights = c(3, 1))

ggsave(
  "Figure_3/chr7_two_panels_points_and_density_colored_overlap.png",
  p_two_panels, width = 12, height = 10, dpi = 220
)
ggsave(
  "Figure_3/chr7_two_panels_points_and_density_colored_overlap.svg",
  p_two_panels, width = 12, height = 10
)

p_two_panels
