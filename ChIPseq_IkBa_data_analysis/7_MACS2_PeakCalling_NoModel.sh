#!/bin/bash

#SBATCH --job-name=MACS2_peak_calling 
#SBATCH --partition=long
#SBATCH --cpus-per-task=2 
#SBATCH --mem=12G
#SBATCH --nodes=1  
#SBATCH --output=logs/MACS2.IkBa.nomodel.out
#SBATCH --error=logs/MACS2.IkBa.nomodel.err
#SBATCH --array=1-6%1

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

# SPECIFY if are you gonna called broad peaks? Comment the proper line
#BROAD="YES" # Broad peak calling
BROAD="NO" # No broad peak calling (narrow)

# SPECIFY the file name where the sample;input is included. Remember to create
# one txt per type of calling (i.e. broad or narrow) according to previous
# configuration. Include a Return in the last row file!
SAMPLESHEET=$WKD"/MACS2_peak_calling/Samples_Input_IkBa.txt"

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
##       MACS2 peak calling
################################################################################

# MACS2 peak caller requires aligned reads in BAM format (typically from Bowtie2)
# For ChIPseq: BAM should only include unique mappers and duplicates should be marked

# Link to MACS2 peak caller manual
# https://pypi.org/project/MACS2/

###########################
## 1. Other relevant paths
###########################

# Folder where input BAM files are available
DATA=$WKD'/Bowtie_align/BAM_Markdup'

# Folder where MACS2 output results will be stored: 
OUTPEAKS=$WKD'/MACS2_peak_calling/Other_results'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
MACS2='macs2_v2.2.7.1.sif '  #This image inludes MACS2 2.2.7.1

# Specify any particular tool parameters

# Effective genome size. MACS2 has precomputed values 
# Recommended to use these ones adapted to read length 
# https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
GSIZE=2701495711 # hg38 50bp read length 

# Adj p-val (q-val) to be used as threshold criteria: 5% by default
FDR=0.05

# Criteria (cut-off) for broad regions, if considered: 10% by default
BROAD_CUTOFF=0.1

# NOTE: MACS2 does not consider duplicates for peak calling
KEEPDUP="" # "--keep-dup all" if you wanna change this behaviour

################################################################################
## 3. Command file preparation: to execute batch array
################################################################################

while IFS=";" read -r sample input fragment; 
do
  # Sample name
  NAME=${sample%_trimmed.sorted.unique.markdup.filtered_blacklisted.bam}
  
  # Peak calling with MACS2:
  if [ $BROAD = "NO" ]; then
    echo "singularity exec $IMAGES_PATH/$MACS2 macs2 callpeak -B $KEEPDUP --nomodel --extsize $fragment -g $GSIZE -q $FDR -t $DATA/$sample -c $DATA/$input --outdir $OUTPEAKS -n $NAME"
  else
    echo "singularity exec $IMAGES_PATH/$MACS2 macs2 callpeak -B $KEEPDUP --nomodel --broad --broad-cutoff $BROAD_CUTOFF -g $GSIZE -q $FDR -t $DATA/$sample -c $DATA/$input --outdir $OUTPEAKS -n $NAME'_nomodel_relaxed_adjpval' "
  fi
  
done < $SAMPLESHEET > $WKD'/scripts/cmds/MACS2_peak_calling_samples.cmd'


################################################################################
## 4. MACS2 Peak calling
################################################################################

echo "-----------------------------------------"
echo "Starting MACS2 Peak Calling"
echo "-----------------------------------------"

DATE=$(date +%m-%d-%Y--%T)
echo "  Samples peak calling in array mode: $DATE"
echo " "

SEEDFILE=$WKD'/scripts/cmds/MACS2_peak_calling_samples.cmd'
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
eval $SEED

DATE=$(date +%m-%d-%Y--%T)
echo "  All samples peak-called: $DATE"


################################################################################
## 5. Move files different from the resulting BED files (narrow or broad Peak files)
################################################################################

if [ $BROAD = "NO" ]; then
	  mv $OUTPEAKS/*.narrowPeak $WKD'/MACS2_peak_calling/Peaks'
else
  mv $OUTPEAKS/*.broadPeak $WKD'/MACS2_peak_calling/Peaks'
fi

################################################################################
## 6. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'MACS2 peak calling completed' 
echo "Processing Time: $DIFF seconds"


