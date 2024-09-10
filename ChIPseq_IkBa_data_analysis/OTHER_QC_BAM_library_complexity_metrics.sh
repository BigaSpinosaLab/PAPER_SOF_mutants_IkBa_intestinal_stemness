#!/bin/bash

#SBATCH --job-name=LC_Metrics
#SBATCH --partition=long
#SBATCH --cpus-per-task=2 
#SBATCH --mem=8G
#SBATCH --nodes=1  
#SBATCH --output=logs/LC_Metrics.out
#SBATCH --error=logs/LC_Metrics.err

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
##       QC on BAM files from ChIPseq experiment: 
##  Metrics: NRF, PBC1 and PBC2 as defined in ENCODE 
##  https://www.encodeproject.org/data-standards/terms/#library
################################################################################

# Use of Custom Python script from PEPATAC (ATAC-seq analysis)
# This custom script covers these three metrics from the obtained BAM file

BAMQC=$IMAGES_PATH/'pepatac/tools/bamQC.py'

###########################
## 1. Other relevant paths
###########################

# Folder where BAM files are stored: the ones with dups
DATA=$WKD'/Bowtie_align/BAM_Markdup'

# Folder where Library Complexity metrics are stored
OUT_BAMQC=$WKD'/Bowtie_align/QC/Library_Complexity'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# We need to call the PEPATAC images to load the proper Python modules

# Specify image/s name to be used (tool-related)
PEPATAC='pepatac_v0.10.3.sif' 

################################################################################
## 3. Execution in a loop
################################################################################

# REMARK: Execution in batch array mode was failing

# Create a running pepatac instance for all samples
singularity instance start -B $ROOTDIR:$ROOTDIR $IMAGES_PATH/$PEPATAC pepatac_instance

for FILENAME in $DATA/*_trimmed.sorted.unique.markdup.filtered_blacklisted.bam
do
     NAME=${FILENAME%.bam}
     SAMPLE=$(basename $NAME)

     # LC metrics execution using PEPATAC script
     
    singularity exec instance://pepatac_instance \
        $BAMQC -i $FILENAME \
        --cores 6 \
        -o $OUT_BAMQC/$SAMPLE.bamQC.tsv

done

# Close the running pepatac instance
singularity instance stop pepatac_instance


################################################################################
## 5. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'BAM library complexity metrics completed' 
echo "Processing Time: $DIFF seconds"
