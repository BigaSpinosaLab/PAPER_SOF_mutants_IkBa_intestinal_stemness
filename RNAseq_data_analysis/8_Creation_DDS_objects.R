################################################################################
##        Title : RNAseq data analysis - Induced ikB
##  Description : This script is for creating the raw and normalized (VST) 
##                DESeqDataSet objects required for downstream analyis (EA, DEA..)
##   Researcher : D.Alvarez-Villanueva
##         Date : 25th May 2022
################################################################################

# Set up paths ------------------------------------------------------------

################################################################################
## 1. Set up paths
################################################################################

.wd <- getwd()
.res <- file.path(.wd, "results")
.dat <- file.path(.wd,"data")
.RData <- file.path(.wd, "RData")

# Load packages -----------------------------------------------------------

################################################################################
## 2. Load packages
################################################################################

require("openxlsx")
require("DESeq2")

# Import counts and samples data ------------------------------------------

################################################################################
## 3. Import data: counts and samples information
################################################################################

# Import count data file
counts <- read.delim(file=file.path(.dat,"RNAseq_INDs_ikBa_gene_quantification.txt"), 
           header=TRUE, sep="\t", skip=1)
counts <- counts[,-which(colnames(counts) %in% c("Chr","Start","End","Strand","Length"))]

# Import metadata files
metadata <- read.table(file=file.path(.dat,"Targets.txt"), 
                       header=TRUE, sep="\t")

# metadata

# File.Name	Sample.id	Cell.type	Condition
# X.projects.cancer.RNAseq_INDs_ikBa_DANI.STAR_align.BAM.108CT1_04839AAD_TGACTTGG_Aligned.sortedByCoord.out.bam	108CT1	108	CT
# X.projects.cancer.RNAseq_INDs_ikBa_DANI.STAR_align.BAM.108CT2_04840AAD_TAAGATGG_Aligned.sortedByCoord.out.bam	108CT2	108	CT
# X.projects.cancer.RNAseq_INDs_ikBa_DANI.STAR_align.BAM.108D1_04841AAD_CCGCAACG_Aligned.sortedByCoord.out.bam	108D1	108	D
# X.projects.cancer.RNAseq_INDs_ikBa_DANI.STAR_align.BAM.108D2_04842AAD_CTGATAAG_Aligned.sortedByCoord.out.bam	108D2	108	D
# X.projects.cancer.RNAseq_INDs_ikBa_DANI.STAR_align.BAM.M3CT1_04843AAD_AATCCGTC_Aligned.sortedByCoord.out.bam	M3CT1	M3	CT
# X.projects.cancer.RNAseq_INDs_ikBa_DANI.STAR_align.BAM.M3CT2_04844AAD_AGAACCGC_Aligned.sortedByCoord.out.bam	M3CT2	M3	CT
# X.projects.cancer.RNAseq_INDs_ikBa_DANI.STAR_align.BAM.M3D1_04845AAD_GCAACGCC_Aligned.sortedByCoord.out.bam	M3D1	M3	D
# X.projects.cancer.RNAseq_INDs_ikBa_DANI.STAR_align.BAM.M3D2_04846AAD_CAACTACC_Aligned.sortedByCoord.out.bam	M3D2	M3	D
# X.projects.cancer.RNAseq_INDs_ikBa_DANI.STAR_align.BAM.WTCT1_04835AAD_TGGTTCTT_Aligned.sortedByCoord.out.bam	WTCT1	WT	CT
# X.projects.cancer.RNAseq_INDs_ikBa_DANI.STAR_align.BAM.WTCT2_04836AAD_TCGCATCT_Aligned.sortedByCoord.out.bam	WTCT2	WT	CT
# X.projects.cancer.RNAseq_INDs_ikBa_DANI.STAR_align.BAM.WTD1_04837AAD_GATGGTTG_Aligned.sortedByCoord.out.bam	WTD1	WT	D
# X.projects.cancer.RNAseq_INDs_ikBa_DANI.STAR_align.BAM.WTD2_04838AAD_AGAGGTTG_Aligned.sortedByCoord.out.bam	WTD2	WT	D


# Create DSeqDataSet object -----------------------------------------------

################################################################################
## 4. Create DSeqDataSet object
################################################################################

# 1. Transform to factors the main effects to be included in the GLM. Define correctly
#    the reference level of interest. An example below: 

# Define the reference level for the condition (CT)
metadata$Condition <- relevel(as.factor(metadata$Condition),ref="CT")


# 2. Generate the DESeqDataSet from HTseq counts
data.dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                           colData = metadata, 
                                           design =  ~ Cell.type + Cell.type:Condition, 
                                           tidy=TRUE)


dim(data.dds) 
# [1] 61552    12
colData(data.dds)
head(assay(data.dds,"counts"))

# 3. Remove low-count genes 

# Remove those with less than 10 counts accross all samples (as implemented in RNAseq_fetal_Ruiz)
data.dds.red <- data.dds[rowSums(counts(data.dds)) >=10]
nrow(data.dds.red)   
# [1] 22993

# Read Counts transformation ----------------------------------------------

################################################################################
## 6. Read Counts transformation: VST 
################################################################################

# VST transformation: less sensitive to high count outliers. 
#  Not blind to the experimental design: otherwise dispersion would be overestimated
vstData.dds<- vst(data.dds.red, blind=FALSE) 

# Store RData objects -----------------------------------------------------

################################################################################
## 7. Store RData objects
################################################################################

# Save DDS objects
save(data.dds.red, file = file.path(.RData,"Raw_DDS_INDs_ikB.RData"))
# Save normalized objects
save(vstData.dds, file= file.path(.RData,"VST_DDS_INDs_ikB.RData"))

# Expression matrices -----------------------------------------------------

################################################################################
## 8. Store csv files with expression matrices
################################################################################

hs <- createStyle(textDecoration = "BOLD")
l <- list("HUMAN_RAWcounts"= assay(data.dds.red),
          "HUMAN_NORMALIZEDcounts" = assay(vstData.dds))

write.xlsx(x = l, 
          file = file.path(.res,"RNAseq_INDs_iKB_expression_matrices.xlsx"), 
           rowNames=TRUE,
           headerStyle= hs)