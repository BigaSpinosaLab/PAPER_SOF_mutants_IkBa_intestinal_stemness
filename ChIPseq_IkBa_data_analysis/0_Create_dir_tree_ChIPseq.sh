#!/bin/bash

################################################################################
##       Create Project directory tree for an standard ChIP-seq analysis with Bowtie2 and MACS2
################################################################################

# Project directory tree consists of (inside main project folder):
#   
#    .- raw_data: Includes raw data (FASTQ files). Subdirs: FASTQC and MULTIQC
#    .- trimmed_data: Includes trimmed data (from i.e. trimgalore). Subdirs: FASTQC and MULTIQC
#    .- Bowtie_align: Includes output from Bowtie alignment. 
#    .- MACS2_peak_calling: includes output from MACS2 (peak calling)
#    .- scripts: Includes all scripts used to analyze the data
#    .- BigWig: Includes generated bigwig files from BAMs

# IMPORTANT REMARK: It is assumed that an initial scripts folder has been previously where
# this script is run

#=========================
# User defined parameters
#=========================

# SPECIFY full path to project directory in the cluster 
PROJECTDIR="/projects/cancer/IKBa_SOF_mutants_LEspinosa/ChIPseq_Ikba_CRC_cell_lines"

mkdir -p $PROJECTDIR/raw_data/{FASTQC,MULTIQC}
mkdir -p $PROJECTDIR/trimmed_data/{FASTQC,MULTIQC}
mkdir -p $PROJECTDIR/Bowtie_align/{BAM,BAM_Markdup,BAM_NoDups,Other_results,QC/{PhantomPeakTools,Library_Complexity}}
mkdir -p $PROJECTDIR/MACS2_peak_calling/{Peaks,Other_results}
mkdir -p $PROJECTDIR/scripts/{cmds,logs}
mkdir -p $PROJECTDIR/BigWig
