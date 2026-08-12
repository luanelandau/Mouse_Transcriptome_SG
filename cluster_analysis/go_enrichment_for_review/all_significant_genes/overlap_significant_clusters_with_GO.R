# Overlap genes from significant genomic clusters with genes in significant
# GO Biological Process terms for PAR, SM, and SL.
#
# Significant clusters are defined as every cluster containing 3 or more genes.
#
# This script intentionally uses no for loops, user-defined helper functions,
# or if statements.

output_dir <- "/Users/llandau/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/go_enrichment/all_significant_genes"
cluster_file <- "/Users/llandau/Library/CloudStorage/Box-Box/SalivaryGlands_LL/Mouse_Transcriptome_SG/go_enrichment/clusters_SM.csv"

# Read cluster assignments and identify all clusters containing 3 or more genes.
clusters <- read.csv(cluster_file, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
cluster_sizes <- table(clusters$Cluster)
significant_cluster_ids <- as.integer(names(cluster_sizes)[cluster_sizes >= 3])
significant_cluster_genes <- clusters[clusters$Cluster %in% significant_cluster_ids, ]
significant_cluster_genes$Cluster_size <- as.integer(cluster_sizes[as.character(significant_cluster_genes$Cluster)])
significant_cluster_genes <- significant_cluster_genes[order(significant_cluster_genes$Cluster,
                                                               significant_cluster_genes$Chr,
                                                               significant_cluster_genes$Start), ]

#write.table(significant_cluster_genes,
#            file.path(output_dir, "significant_clusters_size3_or_more.tsv"),
#            sep = "\t", quote = FALSE, row.names = FALSE)

# PAR: retain significant GO terms and split slash-separated gene symbols.
#par_go <- read.csv(file.path(output_dir, "PAR_sex_all_significant_absLFC1_padj0.05_GO_BP.csv"),
#                   stringsAsFactors = FALSE)
#par_go_significant <- par_go[!is.na(par_go$p.adjust) & par_go$p.adjust < 0.05, ]
#par_go_genes <- sort(unique(unlist(strsplit(as.character(par_go_significant$geneID), "/", fixed = TRUE))))
#par_overlap <- significant_cluster_genes[significant_cluster_genes$Geneid %in% par_go_genes,
#                                         c("Geneid", "Cluster", "Cluster_size", "Chr", "Start", "End")]
#par_overlap <- par_overlap[order(par_overlap$Cluster, par_overlap$Start, par_overlap$Geneid), ]
#
#write.table(par_overlap,
#            file.path(output_dir, "PAR_significant_cluster_GO_overlap.tsv"),
#            sep = "\t", quote = FALSE, row.names = FALSE)
#writeLines(unique(par_overlap$Geneid),
#           file.path(output_dir, "PAR_significant_cluster_GO_overlap_genes.txt"))
#
# SM: retain significant GO terms and split slash-separated gene symbols.
sm_go <- read.csv(file.path(output_dir, "SM_sex_all_significant_absLFC1_padj0.05_GO_BP.csv"),
                  stringsAsFactors = FALSE)
sm_go_significant <- sm_go[!is.na(sm_go$p.adjust) & sm_go$p.adjust < 0.05, ]
sm_go_genes <- sort(unique(unlist(strsplit(as.character(sm_go_significant$geneID), "/", fixed = TRUE))))
sm_overlap <- significant_cluster_genes[significant_cluster_genes$Geneid %in% sm_go_genes,
                                        c("Geneid", "Cluster", "Cluster_size", "Chr", "Start", "End")]
sm_overlap <- sm_overlap[order(sm_overlap$Cluster, sm_overlap$Start, sm_overlap$Geneid), ]

#write.table(sm_overlap,
#            file.path(output_dir, "SM_significant_cluster_GO_overlap.tsv"),
#            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(unique(sm_overlap$Geneid),
           file.path(output_dir, "SM_significant_cluster_GO_overlap_genes.txt"))

## SL: the current GO file has no enriched terms, so these outputs are empty.
#sl_go <- read.csv(file.path(output_dir, "SL_sex_all_significant_absLFC1_padj0.05_GO_BP.csv"),
#                  stringsAsFactors = FALSE)
#sl_go_significant <- sl_go[!is.na(sl_go$p.adjust) & sl_go$p.adjust < 0.05, ]
#sl_go_genes <- sort(unique(unlist(strsplit(as.character(sl_go_significant$geneID), "/", fixed = TRUE))))
#sl_overlap <- significant_cluster_genes[significant_cluster_genes$Geneid %in% sl_go_genes,
#                                        c("Geneid", "Cluster", "Cluster_size", "Chr", "Start", "End")]
#sl_overlap <- sl_overlap[order(sl_overlap$Cluster, sl_overlap$Start, sl_overlap$Geneid), ]
#
#write.table(sl_overlap,
#            file.path(output_dir, "SL_significant_cluster_GO_overlap.tsv"),
#            sep = "\t", quote = FALSE, row.names = FALSE)
#writeLines(unique(sl_overlap$Geneid),
#           file.path(output_dir, "SL_significant_cluster_GO_overlap_genes.txt"))
#
# Compact summary of the results.
#overlap_summary <- data.frame(
#  Gland = c("PAR", "SM", "SL"),
#  Significant_GO_terms = c(nrow(par_go_significant), nrow(sm_go_significant), nrow(sl_go_significant)),
#  Unique_GO_genes = c(length(par_go_genes), length(sm_go_genes), length(sl_go_genes)),
#  Overlapping_cluster_genes = c(nrow(par_overlap), nrow(sm_overlap), nrow(sl_overlap))
#)
#
#write.table(overlap_summary,
#            file.path(output_dir, "significant_cluster_GO_overlap_summary.tsv"),
#            sep = "\t", quote = FALSE, row.names = FALSE)
#
#print(overlap_summary)
#