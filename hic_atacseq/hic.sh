#!/bin/bash
#SBATCH --qos=omergokc
#SBATCH --partition=omergokc
#SBATCH --cluster=faculty
#SBATCH --account=omergokc
#SBATCH --time=00:05:00
#SBATCH --nodes=1
##SBATCH --cpus-per-task=32
#SBATCH --mem=1G
#SBATCH --job-name="heart_cells_TADs"
#SBATCH --output=hic_callTADs.out
#SBATCH --error=hic_callTADs.err
#SBATCH --mail-user=luanejan@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

eval "$(/projects/academic/omergokc/Luane/softwares/anaconda_2025/conda/bin/conda shell.bash hook)"
#conda activate hictk
conda activate herro

#hictk dump -t resolutions *.hic
#10
#20
#50
#100
#200
#500
#1000
#2000
#5000
#10000
#25000
#50000
#100000
#250000
#500000
#1000000
#2500000

## 2. Re-run conversion without norms
#hictk convert \
#  --normalization-methods NONE \
#  ENCFF651SGY.hic \
#  ENCFF651SGY.mcool
#
## Example: ICE on the mcool
#hictk balance ice ENCFF651SGY.mcool
## or SCALE / VC similarly:
## hictk balance scale ENCFF651SGY.mcool
## hictk balance vc ENCFF651SGY.mcool
#
#cooler ls ENCFF651SGY.mcool
#
#hicFindTADs \
#  --matrix ENCFF651SGY.mcool::/resolutions/25000 \
#  --chromosomes chr7 \
#  --minDepth 75000 \
#  --maxDepth 250000 \
#  --step 50000 \
#  --correctForMultipleTesting fdr \
#  --thresholdComparisons 0.05 \
#  --outPrefix TADs_heart_chr7_25kb_Local
#
#
# 3) Call TADs on the balanced 50 kb resolution matrix
#hicFindTADs \
#  --matrix ENCFF651SGY.mcool::/resolutions/50000 \
#  --chromosomes chr7 \
#  --minDepth 150000 \
#  --maxDepth 1500000 \
#  --step 50000 \
#  --correctForMultipleTesting fdr \
#  --thresholdComparisons 0.05 \
#  --outPrefix TADs_heart_chr7_50kb_Canonical
#  
#hicFindTADs \
#  --matrix ENCFF651SGY.mcool::/resolutions/10000 \
#  --chromosomes chr7 \
#  --minDepth 30000 \
#  --maxDepth 100000 \
#  --step 50000 \
#  --correctForMultipleTesting fdr \
#  --thresholdComparisons 0.05 \
#  --outPrefix TADs_heart_chr7_10kb_Local
#
#
#WIN="chr7:43700000-44250000"
#WIN2="chr7:30963161-33800262"
WIN2="chr7:31000000-34000000"

#hicPlotTADs --tracks CH12_25kb.ini --region $WIN2 \
#  --outFileName CH12_TADs_25kb_4.png --dpi 300
  
#hicPlotTADs --tracks 10kbTADs.ini --region $WIN2 \
#  --outFileName heart_TADs_10kb.png --dpi 300
#  
#hicPlotTADs --tracks 25kbTADs.ini --region $WIN2 \
#  --outFileName heart_TADs_25kb.png --dpi 300
#
#hicPlotTADs --tracks 50kbTADs.ini --region $WIN2 \
#  --outFileName heart_TADs_50kb.png --dpi 300


#hicPlotTADs --tracks 25kbTADs_selec.ini --region $WIN2 \
#  --outFileName 25kbTADs_scgb.png --dpi 300
  
#WIN3="chr7:31000000-34000000"
WIN4="chr7:31200000-32000000"

hicPlotTADs --tracks 25kbTADs_selec.ini --region $WIN4 \
  --outFileName 25kbTADs_scgb_try.png --dpi 300


hicPlotMatrix   --matrix /projects/academic/omergokc/Luane/Mouse_SG/heart_cells/ENCFF651SGY.mcool::/resolutions/25000   --region chr7:31000000-31700000   --outFileName test_matrix_25kb_31-31.7mb.png 
hicPlotMatrix   --matrix /projects/academic/omergokc/Luane/Mouse_SG/heart_cells/ENCFF651SGY.mcool::/resolutions/100000   --region chr7:31000000-34000000   --outFileName test_matrix_100kb.png 
hicPlotMatrix   --matrix /projects/academic/omergokc/Luane/Mouse_SG/heart_cells/ENCFF651SGY.mcool::/resolutions/250000   --region chr7:31000000-31700000   --outFileName test_matrix_250kb_31-31.7mb.png 


WIN4="chr7:31200000-32000000"

hicPlotTADs --tracks 25kbTADs_selec.ini --region $WIN4 \
  --outFileName 5kbTADs_scgb_try.png --dpi 300
