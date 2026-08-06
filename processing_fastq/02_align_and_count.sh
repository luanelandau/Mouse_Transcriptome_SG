#!/bin/bash

# Replace the reference index and annotation paths below.
# The same commands were used for all samples from both mouse strains.

HISAT2_INDEX=/path/to/hisat2/mouse_reference
GTF=/path/to/mouse_reference_annotation.gtf

mkdir -p bam featurecounts

for R1 in trimmed_fastq/*_R1_*_val_1.fq.gz
do
    R2=${R1/_R1_/_R2_}
    R2=${R2/_val_1.fq.gz/_val_2.fq.gz}
    SAMPLE=$(basename "$R1" | sed 's/_R1_.*//')

    hisat2 --threads 8 -x "$HISAT2_INDEX" -1 "$R1" -2 "$R2" \
        | samtools sort -@ 8 -o "bam/${SAMPLE}.sorted.bam"

    samtools index "bam/${SAMPLE}.sorted.bam"

    featureCounts -T 8 -p -a "$GTF" \
        -o "featurecounts/${SAMPLE}.counts.txt" \
        "bam/${SAMPLE}.sorted.bam"
done
