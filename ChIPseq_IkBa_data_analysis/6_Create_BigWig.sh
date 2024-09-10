#!/bin/bash

#SBATCH --job-name=BigWigSOF
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --nodes=1  
#SBATCH --output=logs/BigWig_creation.out
#SBATCH --error=logs/BigWig_creation.err
#SBATCH --array=1-8%2

# REMARK!!!! Adapt the number of array tasks to the number of samples i.e. if you have
# 12 samples to analyze you need to specify 12 as indicated. %X means refers to the number 
# of tasks will be sent to cluster execution simultaneously. 
# Each task meaning one sample to be analyzed

#=========================
# User defined parameters: relevant paths and analysis type
#=========================

# SPECIFY Root directory in the cluster (usually /projects/cancer)
ROOTDIR="/projects/cancer"

# SPECIFY your project working directory. 
WKD=$ROOTDIR'/IKBa_SOF_mutants_LEspinosa/ChIPseq_Ikba_CRC_cell_lines'

# SPECIFY the file name where the 'sample;input' is included. 
# This is required if each sample wants to be normalized with its input
# Always left an empty line at the end of the file
SAMPLESHEET=$WKD"/BigWig/Samples_Input_IkBa_BW.txt"

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
##       Create BigWig files with bamCoverage from deepTools
################################################################################

# Take alignment of reads or fragments (BAM format) and generate a coverage track
# in bigWig format.
# BigWig will be individually created (not compared i.e. to an input control)

# Link to deepTools > bamCoverage 
# https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html

###########################
## 1. Other relevant paths
###########################

# Folder where input BAM files are available 
DATA=$WKD'/Bowtie_align/BAM_Markdup' # BAM with marked duplicates so can be ignored

# Folder where BigWig files will be stored: 
OUTBW=$WKD'/BigWig'

# (OPTIONAL) BlackList FileName (BED or GTF)
# Not needed if BAM files have already been filtered of blacklisted regions
#BLACK=

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
DEEPTOOLS='deepTools_v3.5.1.simg'  #This image inludes deepTools v3.5.1

# Specify any particular tool parameters
# Recommended to use these ones adapted to read length 
# https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
GSIZE=2701495711 # hg38 50bp read length 

# Type of normalization to be used: Not used, we will use scaleFactor (readCount)
# RPGC normalizes to 1x coverage
NORM=RPGC  # Other options: CPM, RPKM, BPM, RPGC (not supported for BamCompare) or None

# Number of processors
T=4

# BinSize (by default is 50)
BS=10

# Smoothing (should be larger than BinSize)
SMOOTH=30

################################################################################
## 3. Command file preparation: to execute batch array
################################################################################

while IFS=";" read -r sample fragment; 
do
  # Sample name
  NAME=${sample%_trimmed.sorted.unique.markdup.filtered_blacklisted.bam}

    # Create BigWig files:
    # --ignoreDuplicates # Include this parameter if duplicates to be ignored
    # --normalizeUsing $NORM -> Not valid for bamCompare (sample to input)
    # --blackListFileName $BLACK #Include this parameter if blacklisted regions have not already filtered
    # --centerReads # We discard this option -> maybe wiser for TFs to visualize a sharper peak

    # FineTuned option for individual Coverage BigWig files
    echo "singularity exec $IMAGES_PATH/$DEEPTOOLS bamCoverage -b $DATA/$sample --smoothLength $SMOOTH  --extendReads $fragment  --normalizeUsing $NORM --effectiveGenomeSize $GSIZE --binSize $BS -of bigwig -p $T -o $OUTBW/$NAME.ext_bs10_smooth30_woCenter_wDups.bw"

done < $SAMPLESHEET > $WKD'/scripts/cmds/Create_BigWig_wDups_bs10_smooth30.cmd'


################################################################################
## 4. Create BigWig files
################################################################################

echo "-----------------------------------------"
echo "Starting BigWig creation"
echo "-----------------------------------------"

DATE=$(date +%m-%d-%Y--%T)
echo "  BigWig creation in array mode: $DATE"
echo " "

SEEDFILE=$WKD'/scripts/cmds/Create_BigWig_wDups_bs10_smooth30.cmd'
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
eval $SEED

DATE=$(date +%m-%d-%Y--%T)
echo "  All bigwig files created: $DATE"


################################################################################
## 5. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'BigWig creation completed' 
echo "Processing Time: $DIFF seconds"


