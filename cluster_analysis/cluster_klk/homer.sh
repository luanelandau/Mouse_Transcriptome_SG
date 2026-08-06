#!/bin/bash

# Identify DNA motifs enriched near the transcription start sites of genes in
# the Klk genomic cluster. HOMER uses promoters from other mouse genes as the
# background and writes the results to the output directory.

GENE_LIST="SM_cluster530.txt"
OUTPUT_DIR="homer_results_SM_cluster530"

findMotifs.pl "$GENE_LIST" mouse "$OUTPUT_DIR" \
    -start -400 \
    -end 100

# Main report: homer_results_SM_cluster530/homerResults.html
