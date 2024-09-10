#!/bin/bash

#SBATCH --job-name=PlotCorr
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --nodes=1  
#SBATCH --output=logs/Plot_Corr.out
#SBATCH --error=logs/Plot_Corr.err

#=========================
# User defined parameters: relevant paths and analysis type
#=========================

# SPECIFY Root directory in the cluster (usually /projects/cancer)
ROOTDIR="/projects/cancer"

# SPECIFY your project working directory. 
WKD=$ROOTDIR'/IKBa_SOF_mutants_LEspinosa/ChIPseq_Ikba_CRC_cell_lines'

#=========================
# General configuration
#=========================
START=$(date +%s)
# Enable Singularity image to look into the general path (equivalent to -B)
export SINGULARITY_BIND=$ROOTDIR 
# Path to images folder in cluster
IMAGES_PATH=$ROOTDIR"/images"
# Path to databases folder in cluster
DB_PATH=$ROOTDIR"/db_files"

################################################################################
##       Plot Correlation between samples with plotCorrelation from deepTools
################################################################################

# First, compute the read coverages for genomic regions for two or more BAM files. 
# The analysis can be performed for the entire genome by running the program in ‘bins’ mode.  

# Correlation or PCA is computed over previous matrix

# Link to deepTools > multiBamSummary
# https://deeptools.readthedocs.io/en/develop/content/tools/multiBamSummary.html

# Link to deepTools > plotCorrelation
# https://deeptools.readthedocs.io/en/develop/content/tools/plotCorrelation.html

# Link to deepTools > plotPCA
# https://deeptools.readthedocs.io/en/develop/content/tools/plotPCA.html

###########################
## 1. Other relevant paths
###########################

# Folder where input BAM files are available 
# Deduplicated BAM files. OR Alternatively, one can set the 
# -ignoreDuplicates parameter on
DATA=$WKD'/Bowtie_align/BAM_Markdup' # BAM with marked duplicates so can be ignored

# Folder where results will be stored: 
OUT=$WKD'/Other'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
DEEPTOOLS='deepTools_v3.5.1.simg'  #This image inludes deepTools v3.5.1

# Specify any particular tool parameters
# Number of processors
T=4

################################################################################
## 3. Compute Read coverages along genome
################################################################################

bamfiles=$(ls $DATA/*R1_001_trimmed.sorted.unique.markdup.filtered_blacklisted.bam)
labels=$(for f in $DATA/*R1_001_trimmed.sorted.unique.markdup.filtered_blacklisted.bam; do basename $f | rev | cut -c94- | rev ;done)

singularity exec $IMAGES_PATH/$DEEPTOOLS multiBamSummary bins \
            --bamfiles $bamfiles --labels $labels -p max \
            --ignoreDuplicates -out $OUT/readCounts.npz

################################################################################
## 4. Plot Correlation in a heatmap manner and PCA
################################################################################

# Plot Spearman Correlation
# singularity exec $IMAGES_PATH/$DEEPTOOLS plotCorrelation \
#             -in $OUT/readCounts.npz \
#             --corMethod spearman \
#             --skipZeros \
#             --plotTitle "Spearman Correlation of Read Counts" \
#             --whatToPlot heatmap \
#             --colorMap RdYlBu --plotNumbers \
#             -o $OUT/Heatmap_SpearmanCorr_readCounts.png
            
# Plot PEarson Correlation
# singularity exec $IMAGES_PATH/$DEEPTOOLS plotCorrelation \
#             -in $OUT/readCounts.npz \
#             --corMethod pearson \
#             --skipZeros \
#             --plotTitle "Pearson Correlation of Read Counts" \
#             --whatToPlot heatmap \
#             --colorMap RdYlBu --plotNumbers \
#             -o $OUT/Heatmap_PearsonCorr_readCounts.png
            
# Plot PCA. By default, it is based on the top1k more variable bins (10kb size)
# singularity exec $IMAGES_PATH/$DEEPTOOLS plotPCA \
#             -in $OUT/readCounts.npz \
#             --transpose \
#             --plotTitle "PCA of Read Counts" \
#             -o $OUT/PCA_readCounts_Top1k_bins.png
            
            
# Plot PCA. Plot it with 100k
singularity exec $IMAGES_PATH/$DEEPTOOLS plotPCA \
            -in $OUT/readCounts.npz \
            --transpose \
            --ntop 100000 \
            --plotTitle "PCA of Read Counts (all)" \
            -o $OUT/PCA_readCounts_Top100k_bins.png

################################################################################
## 5. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'BigWig creation completed' 
echo "Processing Time: $DIFF seconds"


