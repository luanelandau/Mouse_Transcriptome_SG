#!/bin/bash
#SBATCH --qos=omergokc
#SBATCH --partition=omergokc
#SBATCH --cluster=faculty
#SBATCH --account=omergokc
#SBATCH --time=4:00:00
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --job-name="HiC"
#SBATCH --output=hic_sergra_new.out
#SBATCH --error=hic_sergra_new.err
#SBATCH --export=NONE
#SBATCH --mail-user=luanejan@buffalo.edu
#SBATCH --mail-type=ALL

eval "$(/projects/academic/omergokc/Luane/softwares/anaconda_2025/conda/bin/conda shell.bash hook)"
conda activate herro

# Run hicFindTADs at 25 kb
hicFindTADs \
  -m granulosa_25kb.h5 \
  --outPrefix granulosa_25kb_TADs_new \
  --minDepth 125000 \
  --maxDepth 1000000 \
  --step 25000 \
  --thresholdComparisons 0.05 \
  --correctForMultipleTesting fdr \
  --numberOfProcessors 8

hicFindTADs \
  -m sertoli_25kb.h5 \
  --outPrefix sertoli_25kb_TADs_new \
  --minDepth 75000 \
  --maxDepth 1500000 \
  --step 25000 \
  --thresholdComparisons 0.05 \
  --correctForMultipleTesting fdr \
  --numberOfProcessors 8



WIN="chr7:43700000-44250000"

hicPlotTADs \
  --tracks ser_gra_svg_25kb_new.ini \
  --region "$WIN" \
  --outFileName sertoli_gra_25kb_new.png \
  --dpi 300
