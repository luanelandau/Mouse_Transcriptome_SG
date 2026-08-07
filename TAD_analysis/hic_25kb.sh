#!/bin/bash
#SBATCH --qos=omergokc
#SBATCH --partition=omergokc
#SBATCH --cluster=faculty
#SBATCH --account=omergokc
#SBATCH --time=72:00:00
#SBATCH --nodes=1
##SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --job-name="10kb_CH12_TADs"
#SBATCH --output=hic_10kb.out
#SBATCH --error=hic_10kb.err
#SBATCH --export=NONE

eval "$(/projects/academic/omergokc/Luane/softwares/anaconda_new/bin/conda shell.bash hook)"
conda activate herro

# 1) Balance the 1kb resolution (in place)
#cooler balance CH12_25kbp.mcool::/resolutions/25000 

# 2) Call TADs directly on the balanced 25 kb group
hicFindTADs \
  --matrix CH12_25kbp.mcool::/resolutions/25000 \
  --chromosomes chr7 \
  --minDepth 150000 \
  --maxDepth 2500000 \
  --step 50000 \
  --correctForMultipleTesting fdr \
  --thresholdComparisons 0.05 \
  --outPrefix TADs_CH12_chr7_25kb


#WIN="chr7:43700000-44250000"
WIN2="chr7:43650000-44300000"
  
hicPlotTADs --tracks CH12_25kb_select.ini --region $WIN2 \
  --outFileName CH12_TADs_25kb_select2.svg --dpi 300