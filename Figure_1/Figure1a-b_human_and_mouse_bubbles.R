## ================================================================
## Mouse (top) + Human (bottom) bubble figure across 3 salivary glands
## (Parotid, Submandibular, Sublingual) — NO PANCREAS
## - Secreted only
## - Top 30 by Mean TPM per gland
## - Colors by orthology (one-to-one / lineage-specific / other)
## - Bubble area ~ Mean TPM, comparable across all panels
## ================================================================

setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(packcircles)
  library(ggplot2)
  library(ggforce)
  library(patchwork)
  library(stringr)
})

#dir.create("figures", showWarnings = FALSE, recursive = TRUE)

## ---------- Helpers ----------
normalize_orth_mouse <- function(x) {
  x <- ifelse(is.na(x), "", x)
  lab <- tolower(gsub("[^a-z0-9]+", "_", x))
  dplyr::case_when(
    lab %in% c("one_to_one","1_to_1","one_2_one","one-one") ~ "one-to-one",
    lab %in% c("mouse_specific","mouse_only","lineage_specific","no_human_ortholog") ~ "mouse-specific",
    lab %in% c("no_ortholog","none","no_match") ~ "mouse-specific",
    TRUE ~ "other"
  )
}

normalize_orth_human <- function(x) {
  x <- ifelse(is.na(x), "", x)
  lab <- tolower(gsub("[^a-z0-9]+", "_", x))
  dplyr::case_when(
    lab %in% c("one_to_one","1_to_1","one_2_one","one-one") ~ "one-to-one",
    lab %in% c("human_specific","human_only","lineage_specific","no_mouse_ortholog") ~ "human-specific",
    lab %in% c("no_ortholog","none","no_match") ~ "human-specific",
    TRUE ~ "other"
  )
}

# Consistent 3-class palette for both species
orthology3_colors <- c(
  "one-to-one"      = "#a6c9bd",  # green
  "mouse-specific"  = "#934e65",  # brown (mouse lineage-specific)
  "human-specific"  = "#934e65",  # brown (human lineage-specific)
  "other"           = "grey80"
)

# Circle-packing plot with shared TPM scaling
make_bubble <- function(df, gene_col, fill_col, title_text) {
  df <- df %>% arrange(desc(mean_tpm))
  
  # packcircles interprets 'size' as AREA; use mean_tpm for cross-panel comparability
  packing  <- circleProgressiveLayout(df$mean_tpm, sizetype = "area")
  stopifnot(nrow(df) == nrow(packing))
  df <- bind_cols(df, packing)
  
  vertices <- circleLayoutVertices(packing, npoints = 60)
  # Attach fill to vertices by id (id indexes rows of df)
  vertices$fill <- df[[fill_col]][vertices$id]
  
  ggplot() +
    geom_polygon(
      data = vertices,
      aes(x, y, group = id, fill = fill),
      color = NA, alpha = 0.9
    ) +
    geom_text(
      data = df,
      aes(x = x, y = y, label = .data[[gene_col]]),
      size = 2.5, color = "black", fontface = "bold.italic", lineheight = 0.9
    ) +
    scale_fill_manual(
      values = c(
        "one-to-one"     = "#a6c9bd",
        "mouse-specific" = "#934e65",
        "human-specific" = "#934e65",
        "other"          = "grey80"
      ),
      drop = FALSE,
      name = "Orthology"
    ) +
    coord_equal() +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold")
      # keep legends on; patchwork will collect them
    ) +
    labs(title = title_text)
}


## ---------- Mouse: load + prep (PAR, SM, SL) ----------
mouse_files <- c(
  PAR = "miscelaneous_sheets/mouse_expression/PAR_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv",
  SM  = "miscelaneous_sheets/mouse_expression/SM_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv",
  SL  = "miscelaneous_sheets/mouse_expression/SL_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv"
)

prep_mouse <- function(path, gland_code) {
  df <- readr::read_csv(path, show_col_types = FALSE) %>%
    filter(!is.na(secreted), secreted != "No_annotation") %>%
    mutate(
      mouse_gene   = if ("Geneid" %in% names(.)) Geneid else coalesce(.data[["mouse_gene"]], .data[["Geneid"]]),
      mean_tpm     = suppressWarnings(as.numeric(Mean_TPM)),
      orthology3   = normalize_orth_mouse(if ("ortholog_type" %in% names(.)) ortholog_type else NA_character_)
    ) %>%
    filter(!is.na(mean_tpm)) %>%
    arrange(desc(mean_tpm)) %>%
    distinct(mouse_gene, .keep_all = TRUE) %>%        # dedupe by gene ID
    mutate(rank = row_number()) %>%
    filter(rank <= 30) %>%
    transmute(gland = gland_code, gene = mouse_gene, mean_tpm, orthology3)
  df
}

mouse_list <- lapply(names(mouse_files), function(g) prep_mouse(mouse_files[[g]], g))
mouse_all  <- bind_rows(mouse_list)
mouse_all$orthology3 <- factor(mouse_all$orthology3, levels = c("one-to-one","mouse-specific","other"))

################### Added Feb 2026: to count how much of total TPM the 
#top 30 genes make ####################################

PAR = read.csv("miscelaneous_sheets/mouse_expression/PAR_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv")
SM  = read.csv("miscelaneous_sheets/mouse_expression/SM_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv")
SL  = read.csv("miscelaneous_sheets/mouse_expression/SL_mastersheet_TPMs_annotated_for_secretion_and_orthology.csv")

total_tpm_PAR=sum(PAR$Mean_TPM)
total_tpm_SM=sum(SM$Mean_TPM)
total_tpm_SL=sum(SL$Mean_TPM)

PAR_sec <- subset(PAR, !is.na(secreted) & secreted != "No_annotation")
SM_sec <- subset(SM, !is.na(secreted) & secreted != "No_annotation")
SL_sec <- subset(SL, !is.na(secreted) & secreted != "No_annotation")

PAR_sec_top <- PAR_sec %>% arrange(desc(Mean_TPM)) %>% slice_head(n = 30)
SM_sec_top <- SM_sec %>% arrange(desc(Mean_TPM)) %>% slice_head(n = 30)
SL_sec_top <- SL_sec %>% arrange(desc(Mean_TPM)) %>% slice_head(n = 30)

#total of all TPMs
sum(PAR_sec_top$Mean_TPM)/total_tpm_PAR
sum(SM_sec_top$Mean_TPM)/total_tpm_SM
sum(SL_sec_top$Mean_TPM)/total_tpm_SL

#total of all secreted genes
sum(PAR_sec_top$Mean_TPM)/sum(PAR_sec$Mean_TPM)
sum(SM_sec_top$Mean_TPM)/sum(SM_sec$Mean_TPM)
sum(SL_sec_top$Mean_TPM)/sum(SL_sec$Mean_TPM)

#####now for number of genes coding for one-to-one orth

PAR_sec_one <- subset(PAR_sec, !is.na(ortholog_type) & ortholog_type == "one-to-one",)
SM_sec_one <- subset(SM_sec, !is.na(ortholog_type) & ortholog_type == "one-to-one",)
SL_sec_one <- subset(SL_sec, !is.na(ortholog_type) & ortholog_type == "one-to-one",)

#this is to count number of genes that make up for secreted and one-to-one orthologs. 
#I didnt filter for higher TPMs because I wanted to know of total expression... 
PAR_sec_one_=subset(PAR_sec_one, PAR_sec_one$Mean_TPM>0,)
SM_sec_one_=subset(SM_sec_one, SM_sec_one$Mean_TPM>0,)
SL_sec_one_=subset(SL_sec_one, SL_sec_one$Mean_TPM>0,)

sum(PAR_sec_one$Mean_TPM)/total_tpm_PAR
sum(SM_sec_one$Mean_TPM)/total_tpm_SM
sum(SL_sec_one$Mean_TPM)/total_tpm_SL

sum(PAR_sec_one$Mean_TPM)/sum(PAR_sec$Mean_TPM)
sum(SM_sec_one$Mean_TPM)/sum(SM_sec$Mean_TPM)
sum(SL_sec_one$Mean_TPM)/sum(SL_sec$Mean_TPM)

nrow(PAR_sec_one_)
nrow(SM_sec_one_)
nrow(SL_sec_one_)

nrow(PAR_sec_one)
nrow(SM_sec_one)
nrow(SL_sec_one)
######################################################
######################################################






## ---------- Human: load + prep (PAR, SM, SL) ----------
human_files <- c(
  PAR = "miscelaneous_sheets/human expression/PAR_human_mastersheet_TPMs_annotated_with_orthology.csv",
  SM  = "miscelaneous_sheets/human expression/SM_human_mastersheet_TPMs_annotated_with_orthology.csv",
  SL  = "miscelaneous_sheets/human expression/SL_human_mastersheet_TPMs_annotated_with_orthology.csv"
)

prep_human <- function(path, gland_code) {
  df <- readr::read_csv(path, show_col_types = FALSE) %>%
    filter(!is.na(`Secretome location`), `Secretome location` != "No_annotation") %>%
    mutate(
      human_gene   = if ("Geneid" %in% names(.)) Geneid else coalesce(.data[["human_gene"]], .data[["Geneid"]]),
      mean_tpm     = suppressWarnings(as.numeric(Mean_TPM)),
      orthology3   = normalize_orth_human(if ("ortholog_type" %in% names(.)) ortholog_type else NA_character_)
    ) %>%
    filter(!is.na(mean_tpm)) %>%
    arrange(desc(mean_tpm)) %>%
    distinct(human_gene, .keep_all = TRUE) %>%        # collapse paralogs to top TPM
    mutate(rank = row_number()) %>%
    filter(rank <= 30) %>%
    transmute(gland = gland_code, gene = human_gene, mean_tpm, orthology3)
  df
}

human_list <- lapply(names(human_files), function(g) prep_human(human_files[[g]], g))
human_all  <- bind_rows(human_list)
human_all$orthology3 <- factor(human_all$orthology3, levels = c("one-to-one","human-specific","other"))

## ---------- Order glands and build panels ----------
gland_order  <- c("PAR","SM","SL")
gland_titles <- c(PAR = "Parotid", SM = "Submandibular", SL = "Sublingual")

mouse_all$gland <- factor(mouse_all$gland, levels = gland_order)
human_all$gland <- factor(human_all$gland, levels = gland_order)

# Build mouse (top row)
mouse_plots <- lapply(levels(mouse_all$gland), function(g) {
  df_g <- dplyr::filter(mouse_all, gland == g)
  make_bubble(
    df   = df_g,
    gene_col = "gene",
    fill_col = "orthology3",
    title_text = paste0(gland_titles[g], " (Mouse)")
  )
})

# Build human (bottom row)
human_plots <- lapply(levels(human_all$gland), function(g) {
  df_g <- dplyr::filter(human_all, gland == g)
  make_bubble(
    df   = df_g,
    gene_col = "gene",
    fill_col = "orthology3",
    title_text = paste0(gland_titles[g], " (Human)")
  )
})

# Assemble 3x2 grid with a single collected legend
top_row    <- wrap_plots(mouse_plots, nrow = 1)
bottom_row <- wrap_plots(human_plots, nrow = 1)
combined   <- (top_row / bottom_row) + plot_layout(guides = "collect") &
  theme(legend.position = "right")

## ---------- Save ----------
ggsave("Figure_1/bubbles_mouse_top_human_bottom_3glands.svg", combined, width = 14, height = 8)
ggsave("Figure_1/bubbles_mouse_top_human_bottom_3glands.png", combined, width = 14, height = 8, dpi = 300)

# Also save species-only rows if useful
ggsave("Figure_1/bubbles_mouse_row_3glands.png",  top_row  + plot_layout(guides = "collect"),  width = 18, height = 6,  dpi = 300)
ggsave("Figure_1/bubbles_human_row_3glands.png",  bottom_row + plot_layout(guides = "collect"), width = 18, height = 6,  dpi = 300)

message("Done. Wrote: Figure_1/bubbles_mouse_top_human_bottom_3glands.(svg|png)")
