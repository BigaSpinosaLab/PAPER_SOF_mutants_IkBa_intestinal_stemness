#!/bin/bash

#SBATCH --job-name=Black_Regions
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --nodes=1  
#SBATCH --output=logs/Black_Regions.out
#SBATCH --error=logs/Black_Regions.err
#SBATCH --array=1-8%3

# REMARK!!!! Adapt the number of array tasks to the number of samples i.e. if you have
# 12 samples to analyze you need to specify 12 as indicated. %X means refers to the number 
# of tasks will be sent to cluster execution simultaneously. Each task meaning one sample
# to be analyzed

#=========================
# User defined parameters: relevant paths
#=========================
# SPECIFY Root directory in the cluster (usually /projects/cancer)
ROOTDIR="/projects/cancer"

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
##       Remove BlackListed regions: ENCODE Black list
## https://github.com/Boyle-Lab/Blacklist
## These includes: High Signal Regions and Low Mappability region
################################################################################

###########################
## 1. Other relevant paths
###########################

# Folder where the BAM data is (there must be corresponding index file)
DATA=$WKD'/Bowtie_align/BAM_Markdup'

# Folder where to store BAM files without blacklisted regions
OUTBAM=$WKD'/Bowtie_align/BAM_Markdup'

# BED file including blacklisted regions
BL=$DB_PATH'/Blacklisted_Regions/hg38-blacklist.v2_woCHR.bed'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
SAMTOOLS='samtools_v1.15.sif'   
BEDTOOLS='bedtools_v2.30.0.sif'

# Specify any particular tool parameters
T=4  # Number of threads for samtools

################################################################################
## 3. Remove black listed regions: create command
################################################################################

# Command for samtools execution -> easier to read for below command
SAMTOOLS_exec="singularity exec $IMAGES_PATH/$SAMTOOLS samtools"

for FILENAME in $DATA/*.markdup.bam 
do
        NAME=${FILENAME%.bam}
        SAMPLE=$(basename $NAME)
        
        echo "singularity exec $IMAGES_PATH/$BEDTOOLS bedtools intersect -abam $FILENAME -b $BL -v | $SAMTOOLS_exec sort -O BAM | tee $OUTBAM/$SAMPLE.filtered_blacklisted.bam | $SAMTOOLS_exec index - $OUTBAM/$SAMPLE.filtered_blacklisted.bam.bai"
        
done > $WKD'/scripts/cmds/Remove_BL.cmd'

################################################################################
## 4. Execute removing black-listed regions
################################################################################

echo "-----------------------------------------"
echo "Creating BAM index"
echo "-----------------------------------------"

DATE=$(date +%m-%d-%Y--%T)
echo "  Creating BAM index in array mode: $DATE"
echo " "

SEEDFILE=$WKD'/scripts/cmds/Remove_BL.cmd'
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
eval $SEED

DATE=$(date +%m-%d-%Y--%T)
echo "Blacklisted regions removed: $DATE"

################################################################################
## 5. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'BL removed' 
echo "Processing Time: $DIFF seconds"
