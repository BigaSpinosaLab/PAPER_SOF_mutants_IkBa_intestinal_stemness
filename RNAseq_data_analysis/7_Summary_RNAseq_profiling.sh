#!/bin/bash

#SBATCH --job-name=Summary_RNAseq
#SBATCH --partition=long
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1  
#SBATCH --output=Summary_RNAseq.out
#SBATCH --error=Summary_RNAseq.err

#=========================
# User defined parameters: relevant paths
#=========================
# SPECIFY Root directory in the cluster (usually /projects/cancer)
ROOTDIR="/projects/cancer"

# Specify project name. It will be used in final reports title
PROJECT='RNAseq_INDs_ikBa'

# SPECIFY your project working directory
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

##################
## 1. Other relevant Paths
##################

# Location where to add the summary output results
OUTMULTIQC=$WKD/'summary'

#####################
## 2. Singularity image and Tool Parametrization
#####################

# Specify image/s name to be used (tool-related)
MULTIQC='multiqc_v1.12.sif' # This images includes MultiQC 1.12

#####################
## 3. Create file with results files
#####################

echo $WKD/'FASTQC' > $OUTMULTIQC/'results.paths.txt'
echo $WKD/'trimmed_data' >> $OUTMULTIQC/'results.paths.txt'
#echo $WKD/'trimmed_data/FASTQC' >> $OUTMULTIQC/'results.paths.txt'
echo $WKD/'STAR_align/Other_results' >> $OUTMULTIQC/'results.paths.txt'
echo $WKD/'quantification' >> $OUTMULTIQC/'results.paths.txt'

##################
## 4. Build the final summary report
##################

echo ''
echo 'Run MULTIQC.............................. '`date`
echo ''

cd $OUTMULTIQC

singularity exec $IMAGES_PATH/$MULTIQC multiqc --file-list results.paths.txt --ignore '*_ReadsPerGene.out.tab' --exclude 'snippy' --title $PROJECT': Summary report transcriptome profiling (RNAseq data)' --filename $PROJECT'_Summary_Transcriptome_Profiling.html'

#######################
## 4. End Preprocessing
#######################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Processing Time: $DIFF seconds"
