#!/bin/bash

#SBATCH --job-name=STAR_alignment 
#SBATCH --partition=long
#SBATCH --cpus-per-task=4 
#SBATCH --mem=32G
#SBATCH --nodes=1  
#SBATCH --output=STAR_alignment.out
#SBATCH --error=STAR_alignment.err
#SBATCH --array=1-12%1

# REMARK!!!! Adapt the number of array tasks to the number of samples i.e. if you have
# 12 samples to trim (24 fastq files) you need to specify 12 as indicated

# NOTE: %X means that how many tasks will be sent to cluster execution simultaneously
# Each task meaning one sample to be aligned

#=========================
# User defined parameters: relevant paths
#=========================
# SPECIFY Root directory in the cluster (usually /projects/cancer)
ROOTDIR="/projects/cancer"

# SPECIFY your project working directory. 
# NOTE: It is assumed that Raw data is within 'raw data' folder
#       and Trimmed data is within 'trimmed_data' folder

WKD=$ROOTDIR'/RNAseq_INDs_ikBa_DANI'

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
##       STAR alignment
################################################################################

# STAR alignment requires a genome INDEX (previously computed) - see STAR_Build_Index script.
# STAR default parameters are optimized for mammalian genomes, however specific parameters
# can be tuned. Check below section (Tool Parametrization) as an example.

# Link to STAR manual
# https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

###########################
## 1. Other relevant paths
###########################

# Folder where data to be aligned is located (expected 'trimmed_data')
DATA=$WKD'/trimmed_data'

# Folder where STAR index is stored
INDEX=$WKD'/STAR_align/Index'

# Folder where STAR alignment results will be stored
OUT=$WKD'/STAR_align/Other_results'

# Folder where BAM files will be finally stored
OUTBAM=$WKD'/STAR_align/BAM'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
STAR='star_2.7.8a.sif'  #This image inludes STAR 2.7.8
SAMTOOLS='samtools_v1.15.sif' # This image includes SAMTOOLS 1.15  

# Specify any particular tool parameters
# Number of threads                         
T='8' 

# Number of allowed mismatches 
MISMATCH='1' # It could be reduce to 1 or even 0 for human/mouse genomes

# Number of multimappers allowed (1==unique mapping reads)
MULTIMAP="1"

################################################################################
## 3. Command file preparation: to execute batch array
################################################################################

# File suffix to be adapted!

for FILENAME in $DATA/*_R1_001_val_1.fq.gz 
do
    NAME=${FILENAME%_R1_001_val_1.fq.gz}
    SAMPLE=$(basename $NAME)

    # Forward and Reverse Reads for that sample
	READ1=$NAME'_R1_001_val_1.fq.gz'  
	READ2=$NAME'_R2_001_val_2.fq.gz'  
    
    # Construct the full execution command for STAR alignment
    echo singularity exec $IMAGES_PATH/$STAR STAR --runThreadN $T \
                                                --genomeDir $INDEX \
                                                --genomeLoad LoadAndKeep \
						--limitBAMsortRAM 32000000000 \
                                                --readFilesIn $READ1 $READ2 \
                                                --readFilesCommand zcat \
                                                --outSAMtype BAM SortedByCoordinate \
                                                --outFilterMismatchNmax $MISMATCH \
                                                --outFilterMultimapNmax $MULTIMAP \
                                                --quantMode GeneCounts  \
                                                --outFileNamePrefix $OUT/$SAMPLE'_'

# Other interesting parameters that could be added
#--outReadsUnmapped Fastx \  #Uncomment if you want to store unmapped reads

done > $WKD'/scripts/STAR_align_samples.cmd'


################################################################################
## 4. STAR alignment: Genome load and samples alignment
################################################################################

echo "-----------------------------------------"
echo "Starting Alignment to reference genome: Loading genome index"
echo "-----------------------------------------"

singularity exec $IMAGES_PATH/$STAR STAR --genomeDir $INDEX --genomeLoad LoadAndExit

DATE=$(date +%m-%d-%Y--%T)
echo "  Samples alignment in array mode: $DATE"
echo " "

SEEDFILE=$WKD'/scripts/STAR_align_samples.cmd'
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
eval $SEED

DATE=$(date +%m-%d-%Y--%T)
echo "  All samples aligned: $DATE"

################################################################################
## 5. Clean up: remove index, move BAM files and create respective bam index
################################################################################

echo "  Remove index from memory"
echo " "
singularity exec $IMAGES_PATH/$STAR STAR --genomeDir $INDEX --genomeLoad Remove

echo " Move BAM files to specific folder"
echo " "

for BAM in $OUT/*.bam
do
	mv $BAM $OUTBAM
done

echo " Create BAM index"
for BAM in $OUTBAM/*.bam
do
	singularity exec $IMAGES_PATH/$SAMTOOLS samtools index $BAM
done


################################################################################
## 4. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'STAR alignment completed' 
echo "Processing Time: $DIFF seconds"
