#!/bin/bash

#SBATCH --job-name=BAM_RemoveDups
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --nodes=1  
#SBATCH --output=logs/BAM_RemoveDUPS.out
#SBATCH --error=logs/BAM_RemoveDUPS.err
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
##       Remove Duplicates, already marked, in BAM files using Samtools
################################################################################

###########################
## 1. Other relevant paths
###########################

# Folder where the mark-duplicated BAM data is located
DATA=$WKD'/Bowtie_align/BAM_Markdup'

# Folder where BAM files Remove DUP will be stored
OUTBAM=$WKD'/Bowtie_align/BAM_NoDups'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
SAMTOOLS='samtools_v1.15.sif' # This image includes SAMTOOLS 1.15  

# Specify any particular tool parameters
T=4  # Number of threads for samtools

################################################################################
## 3. Mark Duplicates: create command and execute
################################################################################

# Command for samtools execution -> easier to read for below command
SAMTOOLS_exec="singularity exec $IMAGES_PATH/$SAMTOOLS samtools"

for FILENAME in $DATA/*.markdup.filtered_blacklisted.bam
do
    NAME=${FILENAME%.markdup.filtered_blacklisted.bam}
    SAMPLE=$(basename $NAME)

    # Create a new bam file without duplicates for QC purposes
    echo "$SAMTOOLS_exec view -F 1024 -bh $FILENAME > $OUTBAM/$SAMPLE.NoDups.filtered_blacklisted.bam"

done > $WKD'/scripts/cmds/BAM_remove_dups.cmd'

################################################################################
## 4. Execute marking duplicates
################################################################################

echo "-----------------------------------------"
echo "Starting Removing Duplicates in BAM Files"
echo "-----------------------------------------"

DATE=$(date +%m-%d-%Y--%T)
echo "  Samples Removing Duplicates in array mode: $DATE"
echo " "

SEEDFILE=$WKD'/scripts/cmds/BAM_remove_dups.cmd'
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
eval $SEED

DATE=$(date +%m-%d-%Y--%T)
echo "  All samples with duplicates removes: $DATE"

################################################################################
## 5. Create index for BAM mark-duplicates and remove-duplicates
################################################################################

echo "-----------------------------------------"
echo "Create index from BAM Files"
echo "-----------------------------------------"

DATE=$(date +%m-%d-%Y--%T)
echo "  Create index: $DATE"
echo " "

for FILENAME in $OUTBAM/*.NoDups.bam
do
    echo "$SAMTOOLS_exec index -@ $T $FILENAME"
    
done > $WKD'/scripts/cmds/BAM_index.cmd'

SEEDFILE=$WKD'/scripts/cmds/BAM_index.cmd'
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
eval $SEED

DATE=$(date +%m-%d-%Y--%T)
echo "  All samples with index created: $DATE"


################################################################################
## 6. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'Samtools Remove Duplicates completed' 
echo "Processing Time: $DIFF seconds"
