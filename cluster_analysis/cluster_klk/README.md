# HOMER motif analysis of the Klk cluster

This analysis searches for DNA sequence motifs enriched in promoters of genes in the Klk genomic cluster.

## Input

`SM_cluster530.txt` is a one-column list of 22 gene identifiers in cluster 530, which contains the kallikrein (`Klk`) locus. The cluster was identified from submandibular-gland sex-associated DEGs using `bedtools cluster -d 100000`; therefore, consecutive genes in the cluster are separated by no more than 100 kb. The clustering procedure and permutation analysis are documented in `../permutations/README.md`.

HOMER recognized 19 of the 22 identifiers as target genes in the saved analysis. It searched the region from 400 bp upstream to 100 bp downstream of each transcription start site and compared these promoters with background mouse promoters.

## Script

`homer.sh` shows the essential HOMER de novo motif-enrichment command. Cluster-specific computing-environment setup and downstream motif-scanning commands were omitted to keep the analysis easy to read.

## Main output

`homerResults.html` is the main de novo motif report. It ranks enriched motifs and shows their sequence logos, enrichment statistics, target/background frequencies, and best matches to known motifs; the top motif was matched to GRHL2.
