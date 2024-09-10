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

# Load packages -----------------------------------------------------------

################################################################################
## 2. Load packages
################################################################################

library("openxlsx")
library("DESeq2")
require(pls)
require(reshape)
require(ggplot2)
require(plotly)
library("ComplexHeatmap")
require(vsn)
require(pvclust)


# Import VST data ---------------------------------------------------------

################################################################################
## 3. Import VST data elsewhere generated  + annotations
################################################################################

# VST data
load(file = "RData/VST_DDS_INDs_ikB.RData")

# Import targets file -----------------------------------------------------

################################################################################
## 4. Import targets file
################################################################################

# Targets
metadata <- read.table(file="data/Targets.txt", 
                       header=TRUE, sep="\t")

# Change the exprs matrix colnames to sample IDs
colnames(vstData.dds) <- metadata$Sample.id

# Hierarchical clustering -------------------------------------------------

################################################################################
## PLOT: Hierarchical clustering
################################################################################

# 1.Correlation coefficient as (dis)similarity distance 
distcorr.vst <- as.dist(1-cor(assay(vstData.dds), method="pearson"))

pdf(file = file.path(.fig,"HClust_VST_Inds_IkB.pdf"), width=6, height=6)
plot(hclust(distcorr.vst, 
           method = "complete"),
     labels=colnames(assay(vstData.dds)),
     main="Clustering based on Pearson correlation (Complete linkage)")
dev.off()

# PCA ---------------------------------------------------------------------

################################################################################
## PLOT: PCA
################################################################################

# 1. Using implemented DESeq2 function: quick view. 
plotPCA(vstData.dds,intgroup=c("Cell.type", "Condition"),) #, ntop=nrow(vstData.dds))

# 2. Manual computation with prcomp
mat <- t(assay(vstData.dds)) 
pca <-  prcomp(mat, scale=TRUE, center=TRUE)
PCs <- as.data.frame(pca$x)
PCs$Sample <- rownames(PCs)
PCs$Condition <- as.character(colData(vstData.dds)$Condition)
PCs$Cell.type <- as.character(colData(vstData.dds)$Cell.type)
Variance <- round(summary(pca)$importance[2,]*100, digits=1)

pdf(file = "PAPER_related/PCA_PC1PC2_iIKB_CTL_DOX_AllSamples.pdf", width=4, height=4)
ggplot(PCs, aes(PC1, PC2, color=Cell.type,label=Sample, shape=Condition)) +
  #geom_point(color="black",size=7,alpha=0.3) +
  geom_point(size=3,alpha=0.8) +
  #geom_text(aes(label=Sample),hjust=1, vjust=1,size=4, color="black")+
  xlab(paste0("PC1: ",Variance[1],"% variance")) +
  ylab(paste0("PC2: ",Variance[2],"% variance")) +
  scale_color_manual(values=c("cadetblue4","darkorange4","darkorchid4"))+
  theme_bw()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=14, hjust = 1),
        axis.text.y = element_text(size=14, hjust = 1)) +
  coord_fixed() +
  xlim(-100,150) +
  ylim(-100,100)
dev.off()


pdf(file = "PAPER_related/PCA_PC2PC3_iIKB_CTL_DOX_AllSamples.pdf", width=4, height=4)
ggplot(PCs, aes(PC2, PC3, color=Cell.type,label=Sample, shape=Condition)) +
  #geom_point(color="black",size=7,alpha=0.3) +
  geom_point(size=3,alpha=0.8) +
  #geom_text(aes(label=Sample),hjust=1, vjust=1,size=4, color="black")+
  xlab(paste0("PC2: ",Variance[2],"% variance")) +
  ylab(paste0("PC3: ",Variance[3],"% variance")) +
  scale_color_manual(values=c("cadetblue4","darkorange4","darkorchid4"))+
  theme_bw()+
  theme(legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=14, hjust = 1),
        axis.text.y = element_text(size=14, hjust = 1)) +
  coord_fixed() +
  xlim(-100,150) +
  ylim(-100,100)
dev.off()


# Boxplot -----------------------------------------------------------------

################################################################################
## PLOT:  Boxplot
################################################################################


toplot <- melt(assay(vstData.dds))
colnames(toplot) <- c("Gene","Sample","Transformed_Counts")

toplot$Condition <- rep(colData(vstData.dds)$Condition, 1, each=nrow(vstData.dds))
toplot$Cell.type <- rep(colData(vstData.dds)$Cell.type, 1, each=nrow(vstData.dds))

pdf(file = file.path(.fig,"Expression_Violin_plot_VST_INDs_iKB.pdf"),width=12, height=8)
ggplot(toplot, aes(x=Sample, y=Transformed_Counts, fill=Condition)) + 
  geom_violin() +
  facet_wrap("Cell.type",scale="free_x") +
  scale_fill_manual(values=c("blue","maroon2","forestgreen")) +
  stat_summary(fun=median, geom="point", size=2, color="grey")
dev.off()

