#!/bin/bash

#SBATCH --job-name=FingerPrints
#SBATCH --partition=fast
#SBATCH --cpus-per-task=4 
#SBATCH --mem=24G
#SBATCH --nodes=1  
#SBATCH --output=logs/FingerPrints.out
#SBATCH --error=logs/FingerPrints.err

# REMARK!!!! Adapt the number of array tasks to the number of samples i.e. if you have
# 12 samples to analyze you need to specify 12 as indicated. %X means refers to the number 
# of tasks will be sent to cluster execution simultaneously. Each task meaning one sample
# to be analyzed

#=========================
# User defined parameters: relevant paths
#=========================
# SPECIFY Root directory in the cluster (usually /projects/cancer)
ROOTDIR="/projects/cancer"

# SPECIFY your project working directory. 
# NOTE: It is assumed that Raw data is within 'raw data' folder
#       and Trimmed data is within 'trimmed_data' folder

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
##       QC on BAM files from ChIPseq experiment: plotFingerPrints (deepTools)
################################################################################

# Cumulative enrichment, aka BAM fingerprint, is yet another way of checking the 
# quality of ChIP-seq signal. It determines how well the signal in the ChIP-seq 
# sample can be differentiated from the background distribution of reads in the 
# control input sample.

# Link to plotFingerPrints (deepTools) 
# https://deeptools.readthedocs.io/en/latest/content/tools/plotFingerprint.html

###########################
## 1. Other relevant paths
###########################

# Folder where BAM files (duplicates are automatically ignored)
DATA=$WKD'/Bowtie_align/BAM_Markdup' 

# Folder where phantompeaktools results are stored
OUTPLOT=$WKD'/Bowtie_align/QC/FingerPrints'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
DEEPTOOLS='deepTools_v3.5.1.simg ' # This image includes deepTools v3.5.1

# Specify any particular tool parameters


################################################################################
## 3. Command file preparation: to execute batch array
################################################################################

bamfiles=$(ls $DATA/*R1_001_trimmed.sorted.unique.markdup.filtered_blacklisted.bam)
labels=$(for f in $DATA/*R1_001_trimmed.sorted.unique.markdup.filtered_blacklisted.bam; do basename $f | rev | cut -c94- | rev ;done)

# Plot FingerPrints from deepTools:   # --extendReads
# Extension value is the most representative for all samples

singularity exec $IMAGES_PATH/$DEEPTOOLS plotFingerprint \
          --ignoreDuplicates --bamfiles $bamfiles \
          --labels $labels \
          --extendReads 160 \
          -T 'Fingerprints' -plot $OUTPLOT'/Fingerprints_ChIP_IkBa_HCT_HT29_ext160.pdf'

################################################################################
## 4. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'plot FingerPrints completed' 
echo "Processing Time: $DIFF seconds"