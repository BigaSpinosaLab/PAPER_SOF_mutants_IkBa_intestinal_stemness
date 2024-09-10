#!/bin/bash

#SBATCH --job-name=Bowtie_alignment 
#SBATCH --partition=long
#SBATCH --cpus-per-task=4 
#SBATCH --mem=32G
#SBATCH --nodes=1  
#SBATCH --output=logs/Bowtie_IkBa.out
#SBATCH --error=logs/Bowtie_IkBa.err
#SBATCH --array=1-8%2

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
##       Bowtie2 alignment
################################################################################

# Bowtie2 alignment requires a genome INDEX (previously computed) - It should have
# been  computed for you with Bowtie2 build.

# Link to Bowtie2 aligner manual
# http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-aligner

###########################
## 1. Other relevant paths
###########################

# Folder where proper Bowtie2 index is stored - INCLUDING BASENAME
INDEX=$DB_PATH'/Genomes/Ensembl/human/hg38/release-106/Bowtie2_index/hg38_Ensembl_r106'

# Folder where data to be aligned is located (i.e 'trimmed_data')
DATA=$WKD'/trimmed_data'

# Folder where Bowtie2 alignment results will be stored
OUT=$WKD'/Bowtie_align/Other_results'

# Folder where BAM files will be finally stored
OUTBAM=$WKD'/Bowtie_align/BAM'

# Folder where Alignment stats files will be finally stored
OUTSTATS=$WKD'/Bowtie_align/Other_results'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
BOWTIE='bowtie2_v2.4.4.sif '  #This image inludes Bowtie2 2.4.4
SAMTOOLS='samtools_v1.15.sif' # This image includes SAMTOOLS 1.15  

# Specify any particular tool parameters
# Number of threads                         
T='8' 

################################################################################
## 3. Command file preparation: to execute batch array
################################################################################

# File suffix to be adapted!

# Command for samtools execution -> easier to read piped command below
SAMTOOLS_exec="singularity exec $IMAGES_PATH/$SAMTOOLS samtools"

for FILENAME in $DATA/*.fq.gz 
do
    NAME=${FILENAME%.fq.gz}
    SAMPLE=$(basename $NAME)
    
    # Construct the full execution command for Bowtie2 alignment + remove those with MAPQ < 20 or unmapped + unique mappers + create BAM index
    echo "singularity exec $IMAGES_PATH/$BOWTIE bowtie2 --threads $T -x $INDEX -U $FILENAME --no-unal 2> $OUTSTATS/$SAMPLE.align.stats.txt | $SAMTOOLS_exec view -F 4 -h -q20 | grep -v 'XS:i:' | $SAMTOOLS_exec sort -O BAM | tee $OUTBAM/$SAMPLE.sorted.unique.bam | $SAMTOOLS_exec index - $OUTBAM/$SAMPLE.sorted.unique.bam.bai"

done > $WKD'/scripts/cmds/Bowtie2_align_samples.cmd'

################################################################################
## 4. Bowtie2 alignment
################################################################################

echo "-----------------------------------------"
echo "Starting Alignment to reference genome"
echo "-----------------------------------------"

DATE=$(date +%m-%d-%Y--%T)
echo "  Samples alignment in array mode: $DATE"
echo " "

SEEDFILE=$WKD'/scripts/cmds/Bowtie2_align_samples.cmd'
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
eval $SEED

DATE=$(date +%m-%d-%Y--%T)
echo "  All samples aligned: $DATE"


################################################################################
## 5. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'Bowtie2 alignment completed' 
echo "Processing Time: $DIFF seconds"
