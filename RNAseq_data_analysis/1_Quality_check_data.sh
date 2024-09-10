#!/bin/bash

#SBATCH --job-name=QC_Trim
#SBATCH --partition=long
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1  
#SBATCH --output=QC.Trim.out
#SBATCH --error=QC.Trim.err

#=========================
# User defined parameters: relevant paths
#=========================
# SPECIFY Root directory in the cluster (usually /projects/cancer)
ROOTDIR="/projects/cancer"

# Specify project name. It will be used in final reports title
PROJECT='RNAseq_INDs_ikBa'

# SPECIFY your project working directory
WKD=$ROOTDIR'/RNAseq_INDs_ikBa_DANI'

# SPECIFY the location to your raw data or trimmed data. 
# Comment or uncomment the relevant one

#TYPE="Raw"
#DATA=$WKD/'raw_data'

TYPE="Trimmed"
DATA=$WKD/'trimmed_data'


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

##################
## 1. Other relevant Paths
##################

# Location where to add the output results (FASTQC and MULTIQC)
OUTFASTQC=$DATA/'FASTQC'
OUTMULTIQC=$DATA/'MULTIQC'

#####################
## 2. Singularity image and Tool Parametrization
#####################

# Specify image/s name to be used (tool-related)
FASTQC='fastqc_0.11.9.sif'  #This image inludes FASTQC 0.11.9
MULTIQC='multiqc_v1.12.sif' # This images includes MultiQC 1.12

# Specify any particular tool parameters
T='10'            # Number of threads

##################
## 3. Quality Assessment
##################

echo ''
echo 'Begin Quality Assessment:FASTQC.............................. '`date`
echo ''

# WARNING: If there are other file types inside $DATA, FASTQC will output an error for them. But
# relevant files (.fastq or .fq will be perfectly inspected)

singularity exec $IMAGES_PATH/$FASTQC fastqc $DATA/* \
				     -t $T \
                                    -o $OUTFASTQC \
                                    --noextract 

echo ''
echo 'Run MULTIQC.............................. '`date`
echo ''

cd $OUTMULTIQC

singularity exec $IMAGES_PATH/$MULTIQC multiqc $OUTFASTQC \
					--title $PROJECT': Quality Assessment of '$TYPE' Reads' \
					--filename $PROJECT'_Summary_QC_'$TYPE'_reads.html'

#######################
## 4. End Preprocessing
#######################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Processing Time: $DIFF seconds"
