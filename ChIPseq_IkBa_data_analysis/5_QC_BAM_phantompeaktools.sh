#!/bin/bash

#SBATCH --job-name=PhantomPeak
#SBATCH --partition=long
#SBATCH --cpus-per-task=2 
#SBATCH --mem=8G
#SBATCH --nodes=1  
#SBATCH --output=logs/Phantom.out
#SBATCH --error=logs/Phantom.err
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
##       QC on BAM files from ChIPseq experiment: PhantomPeakTools
################################################################################

# QC metrics based on Strand cross-correlation. They are based on the fact that 
# a high-quality ChIP-seq experiment produces significant clustering of enriched 
# DNA sequence tags at locations bound by the protein of interest, and that the 
# sequence tag density accumulates on forward and reverse strands centered around 
# the binding site. 

# BAM files should not include: unmapped, multimappers, low quality mappers or 
# duplicates

# Link to PhantomPeakTools github
# https://github.com/kundajelab/phantompeakqualtools

###########################
## 1. Other relevant paths
###########################

# Folder where BAM files already WITHOUT duplicates are stored
DATA=$WKD'/Bowtie_align/BAM_NoDups'

# Folder where phantompeaktools results are stored
OUTPHANTOM=$WKD'/Bowtie_align/QC/PhantomPeakTools'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
PHANTOM='phantompeakqualtools_v1.2.sif' # This image includes PhantomPeakTools v1.2 

# Specify any particular tool parameters

################################################################################
## 3. Command file preparation: to execute batch array
################################################################################

# Command for samtools execution -> easier to read for below command
PHANTOMPEAK_exec="singularity exec $IMAGES_PATH/$PHANTOM Rscript /usr/bin/phantompeakqualtools-1.2/run_spp.R"

for FILENAME in $DATA/*.NoDups.filtered_blacklisted.bam
do
    NAME=${FILENAME%.bam}
    SAMPLE=$(basename $NAME)
    
    # Construct the phantompeaktools execution
    
    echo "$PHANTOMPEAK_exec -c=$FILENAME -savp=$OUTPHANTOM/$SAMPLE.spp.pdf -out=$OUTPHANTOM/$SAMPLE.spp.out"
                            
done > $WKD'/scripts/cmds/PhantomPeakTools.cmd'

################################################################################
## 4. Samtools Mark Duplicates
################################################################################

echo "-----------------------------------------"
echo "Starting QC on ChIPseq BAM files: PhantomPeakTools"
echo "-----------------------------------------"

DATE=$(date +%m-%d-%Y--%T)
echo "  PhantomPeakTools in array mode: $DATE"
echo " "

SEEDFILE=$WKD'/scripts/cmds/PhantomPeakTools.cmd'
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
eval $SEED

DATE=$(date +%m-%d-%Y--%T)
echo "  All samples have been analyzed: $DATE"


################################################################################
## 5. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'BAM phantompeaktools completed' 
echo "Processing Time: $DIFF seconds"
