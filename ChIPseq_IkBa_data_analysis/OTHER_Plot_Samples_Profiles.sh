#!/bin/bash

#SBATCH --job-name=PlotProf
#SBATCH --partition=fast
# #SBATCH --cpus-per-task=2
# #SBATCH --mem=8G
#SBATCH --nodes=1  
#SBATCH --output=logs/Plot_Profiles.out
#SBATCH --error=logs/Plot_Profiles.err

#=========================
# User defined parameters: relevant paths and analysis type
#=========================

# SPECIFY Root directory in the cluster (usually /projects/cancer)
ROOTDIR="/projects/cancer"

# SPECIFY your project working directory. 
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
##       PlotProfile of different samples from deepTools
################################################################################

# First, create a count matrix based on the BigWig files (already normalized) for
# specific bed files and then plot the profile for those

# Required functions from deepTools: computeMatrix and plotProfile

###########################
## 1. Other relevant paths
###########################

# Folder where input BigWig files are available with normalized data
DATA=$WKD'/BigWig'

# Folder where results will be stored: 
OUT=$WKD'/Other/Profiles'

# Txt file where all required information is stored
INFO=$OUT"/Data_x_Profiles.txt"

# Samples from this project
bwfiles=$(ls $DATA/*R1_001.ext_bs10_smooth30_woCenter_wDups.bw)
labels=$(for f in $DATA/*R1_001.ext_bs10_smooth30_woCenter_wDups.bw; do basename $f | rev | cut -c75- | rev ;done)

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
DEEPTOOLS='deepTools_v3.5.1.simg'  #This image inludes deepTools v3.5.1

# Specify any particular tool parameters
# Number of processors
T=4

################################################################################
## 3. Compute matrices considering all bigwig files of interest
################################################################################

while IFS=";" read -r bed name; 
do
    # Compute Matrix 
    
    singularity exec $IMAGES_PATH/$DEEPTOOLS computeMatrix reference-point \
            --referencePoint center \
            -b 3000 -a 3000 \
            --regionsFileName $bed \
            --scoreFileName $bwfiles \
            --samplesLabel $labels \
            --binSize 50 \
            --averageTypeBins mean \
            -p max \
            --skipZeros \
            -out $OUT/$name'_Matrix_peaks.gz' \
            --outFileNameMatrix $OUT/$name'_Matrix_individual_values.tab'
            
      # Plot Profile
      
      singularity exec $IMAGES_PATH/$DEEPTOOLS plotProfile \
                --matrixFile $OUT/$name'_Matrix_peaks.gz' \
                -out $OUT/$name'_Peaks_profile.pdf' \
                --regionsLabel "" \
                --perGroup \
                --refPointLabel $name" center peaks"
                
done < $INFO
            
################################################################################
## 5. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'Plotting profiles completed' 
echo "Processing Time: $DIFF seconds"


