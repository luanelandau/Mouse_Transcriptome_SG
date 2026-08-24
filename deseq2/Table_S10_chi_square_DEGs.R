setwd("~/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/")

deseq_par=read_csv("deseq2/deseq2_results/DESeq2_sex_PAR.csv")
deseq_sm=read_csv("deseq2/deseq2_results/DESeq2_sex_SM.csv")
deseq_sl=read_csv("deseq2/deseq2_results/DESeq2_sex_SL.csv")
deseq_panc=read_csv("deseq2/deseq2_results/DESeq2_sex_PANC.csv")
deseq_liv=read_csv("deseq2/deseq2_results/DESeq2_sex_LIV.csv")

s_par=subset(deseq_par, deseq_par$padj<0.05 & abs(deseq_par$log2FoldChange)>1,)
s_sm=subset(deseq_sm, deseq_sm$padj<0.05 & abs(deseq_sm$log2FoldChange)>1,)
s_sl=subset(deseq_sl, deseq_sl$padj<0.05 & abs(deseq_sl$log2FoldChange)>1,)
s_panc=subset(deseq_panc, deseq_panc$padj<0.05 & abs(deseq_panc$log2FoldChange)>1,)
s_liv=subset(deseq_liv, deseq_liv$padj<0.05 & abs(deseq_liv$log2FoldChange)>1,)

nrow(s_par)
nrow(s_sm)
nrow(s_sl)
nrow(s_panc)
nrow(s_liv)


deg_table <- matrix(
  c(
    nrow(s_par), sum(!is.na(deseq_par$padj))-nrow(s_par),
    nrow(s_sm), sum(!is.na(deseq_sm$padj))-nrow(s_sm),
    nrow(s_sl), sum(!is.na(deseq_sl$padj))-nrow(s_sl),
    nrow(s_panc), sum(!is.na(deseq_panc$padj))-nrow(s_panc),
    nrow(s_liv), sum(!is.na(deseq_liv$padj))-nrow(s_liv)
  ),
  ncol = 2,
  byrow = TRUE
)

rownames(deg_table) <- c("PAR", "SM", "SL", "PANC", "LIV")
colnames(deg_table) <- c("DEG", "Not_DEG")

chisq.test(deg_table)

deg_counts <- c(
  PAR = nrow(s_par),
  SM = nrow(s_sm),
  SL = nrow(s_sl),
  PANC = nrow(s_panc),
  LIV = nrow(s_liv)
)

total_tested <- c(
  PAR = sum(!is.na(deseq_par$padj)),
  SM = sum(!is.na(deseq_sm$padj)),
  SL = sum(!is.na(deseq_sl$padj)),
  PANC = sum(!is.na(deseq_panc$padj)),
  LIV = sum(!is.na(deseq_liv$padj))
)

# DEG proportions and fold differences reported in the manuscript:
# SM is approximately 3.8-fold higher than LIV and 52.8-fold higher than PANC.
deg_proportions <- deg_counts / total_tested

sm_fold_higher <- c(
  LIV = unname(deg_proportions["SM"] / deg_proportions["LIV"]),
  PANC = unname(deg_proportions["SM"] / deg_proportions["PANC"])
)

round(sm_fold_higher, digits = 1)
###

pairwise_prop_test <- pairwise.prop.test(
  x = deg_counts,
  n = total_tested,
  p.adjust.method = "BH"
)

pairwise_p_values <- pairwise_prop_test$p.value
pairwise_p_values_formatted <- matrix(
  format.pval(
    as.vector(pairwise_p_values),
    digits = 2,
    eps = 2e-16,
    na.form = ""
  ),
  nrow = nrow(pairwise_p_values),
  dimnames = dimnames(pairwise_p_values)
)

write.csv(
  pairwise_p_values_formatted,
  "deseq2/Table_S13_chi_square_DEGs.csv",
  row.names = TRUE,
  quote = FALSE
)
