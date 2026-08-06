# RNA-seq FASTQ processing

These scripts summarize the pipeline used to process paired-end mouse RNA-seq
data. The same workflow was applied separately to all samples from both mouse
strains. Large sequencing and alignment files are not included in GitHub.

1. `01_trim_fastq.sh`: Trim Galore removed adapters and low-quality sequence
   and ran FastQC.
2. `02_align_and_count.sh`: HISAT2 aligned the reads to the mouse reference
   genome, SAMtools created sorted BAM files, and featureCounts generated
   gene-level read counts.

Required programs: Trim Galore, FastQC, HISAT2, SAMtools, and featureCounts.
The example reference paths in `02_align_and_count.sh` must be replaced before
running the script.
