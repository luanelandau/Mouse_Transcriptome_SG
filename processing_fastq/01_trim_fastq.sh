#!/bin/bash

# Run this script from the directory containing the FASTQ files.
# The same command was used for all samples from both mouse strains.

mkdir -p trimmed_fastq

for R1 in *_R1_*.fastq.gz
do
    R2=${R1/_R1_/_R2_}

    trim_galore --paired --fastqc --quality 20 --cores 8 \
        --output_dir trimmed_fastq \
        "$R1" "$R2"
done
