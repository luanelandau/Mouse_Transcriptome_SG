#!/bin/bash

###############################################################################
# fastp + Salmon commands used for the mouse RNA-seq analysis
#
# fastp was used to trim and clean the reads, followed by Salmon for transcript
# quantification.
###############################################################################


############################ 1. FASTP #########################################

# Trim adapters and low-quality reads and generate HTML and JSON QC reports.

fastp -i $sample_R1 -I $sample_R2 -o ${sample}.R1.fq.gz -O ${sample}.R2.fq.gz -w "${SLURM_CPUS_PER_TASK}" -h ${sample}.fastp.html -j ${sample}.fastp.json


#################### 2. SALMON DECOY-AWARE INDEX #############################

# Create a decoy-aware index so reads matching genomic regions are not
# incorrectly assigned to transcripts.
gunzip Mus_musculus.GRCm39.cdna.all.fa.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

grep "^>" Mus_musculus.GRCm39.dna.primary_assembly.fa \
  | cut -d " " -f 1 \
  | sed 's/>//' \
  > decoys.txt

cat Mus_musculus.GRCm39.cdna.all.fa \
    Mus_musculus.GRCm39.dna.primary_assembly.fa \
    > gentrome.fa

salmon index \
  -t gentrome.fa \
  -d decoys.txt \
  -p ${SLURM_CPUS_PER_TASK} \
  -i salmon_index \
  --gencode


######################## 3. SALMON QUANTIFICATION #############################

# Quantify the trimmed reads; Salmon infers the library type and corrects for
# sequence and GC bias.
salmon quant -i salmon_index -l A -1 ${sample}.R1.fq.gz -2 ${sample}.R2.fq.gz --validateMappings --seqBias --gcBias -o ${sample} -p ${SLURM_CPUS_PER_TASK}
