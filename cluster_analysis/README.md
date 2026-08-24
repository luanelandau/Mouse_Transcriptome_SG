# Genomic cluster analyses

This folder contains analyses of genomic clustering among sex-associated genes in mouse submandibular gland (SM), including the results summarized in Table S12 and the kallikrein promoter-motif analysis reported in Table S14.

## Contents

- `clusters_SM.csv` lists the observed SM gene clusters.
- `make_clusters_SM_summary.R` summarizes genes by cluster in `clusters_SM_summary.csv`.
- `permutations/` contains the observed-versus-random genomic-clustering workflow and its inputs and outputs; see `permutations/README.md`.
- `cluster_klk/` contains the HOMER promoter-motif analysis of the `Klk` cluster; see `cluster_klk/README.md`.
- `go_enrichment_for_review/` contains GO Biological Process enrichment and cluster-summary scripts used for the review analysis.

Run scripts from the directory expected by each script and review project-specific paths before execution. Detailed requirements and commands are provided in the nested READMEs.
