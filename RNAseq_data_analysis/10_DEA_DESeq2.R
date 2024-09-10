################################################################################
##        Title : RNAseq data analysis - Induced ikB
##  Description : This script is for conducting DEA from dds objects
##   Researcher : D.Alvarez
##         Date : 26th May 2022
################################################################################

# Set up paths ------------------------------------------------------------

################################################################################
## 1. Set up paths
################################################################################

.wd <- getwd()
.res <- file.path(.wd, "results")
.fig <- file.path(.res,"figures")
.dat <- file.path(.wd,"data")
.RData <- file.path(.wd,"RData")

# Load packages -----------------------------------------------------------

################################################################################
## 2. Load packages
################################################################################

library("DESeq2")
library("EnhancedVolcano")
require(ggpubr)
require(org.Hs.eg.db)

# Import dds objects -------------------------------------------------------

################################################################################
## 3. Import data: DESeq2 dataset objects
################################################################################

load(file.path(.RData,"Raw_DDS_INDs_ikB.RData"))

# Assign the object of interest to the next variable
data.dds <- data.dds.red
rm(data.dds.red)

# Import the VST norm values for further plotting
load(file.path(.RData, "VST_DDS_INDs_ikB.RData"))

# DEA ---------------------------------------------------------------------

################################################################################
## 4. DEA with DESeq2 and store results 
## Analyze the treatment effect on every cell type
################################################################################

# Let's relevel the group factor to get the comparison of interest: otherwise 
# apglm logFC shrinkage cannot be applicable

dea.dds <- DESeq(data.dds)
resultsNames(dea.dds)  #Identify the name of the comparison of interest

# [1] "Intercept"               "Cell.type_M3_vs_108"     "Cell.type_WT_vs_108"     "Cell.type108.ConditionD" "Cell.typeM3.ConditionD" 
# [6] "Cell.typeWT.ConditionD" 

# IMPORTANT: Execute complete_dea() function in the next section since it used
# here!!!!!

dea.108 <- complete.dea(DESeqSet = dea.dds,
                        comparison = "Cell.type108.ConditionD",
                        lim.ma.abs = 10,
                        fig.path = .fig)

# out of 22993 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 517, 2.2%
# LFC < 0 (down)     : 995, 4.3%
# outliers [1]       : 0, 0%
# low counts [2]     : 11145, 48%
# (mean count < 83)

dea.M3 <- complete.dea(DESeqSet = dea.dds,
                        comparison = "Cell.typeM3.ConditionD",
                        lim.ma.abs = 10,
                        fig.path = .fig)

# out of 22993 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3, 0.013%
# LFC < 0 (down)     : 17, 0.074%
# outliers [1]       : 0, 0%
# low counts [2]     : 21843, 95%
# (mean count < 6413)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

dea.WT <- complete.dea(DESeqSet = dea.dds,
                       comparison = "Cell.typeWT.ConditionD",
                       lim.ma.abs = 10,
                       fig.path = .fig)

# out of 22993 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 132, 0.57%
# LFC < 0 (down)     : 179, 0.78%
# outliers [1]       : 0, 0%
# low counts [2]     : 12928, 56%
# (mean count < 262)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


# Write results into the same excel file
final.res <- list("DEA_108_D_vs_CT" = dea.108,
                  "DEA_M3_D_vs_CT" = dea.M3,
                  "DEA_WT_D_vs_CT" = dea.WT)

# Save DDS objects
save(final.res, file = file.path(.RData,"DEA_RESULTS_INDs_ikB.RData"))

hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "#4F80BD")
write.xlsx(final.res,file = file.path(.res,"DEA_results_INDs_iKB_CellTypes_DvsCT_Comparisons.xlsx"), headerStyle = hs)


# Function to compute DEA and lfcShrink and add annotations ---------------

################################################################################
## 6. Function to compute DEA, lfcShrinkage, gather all information and add annotations
##    into the same data frame
################################################################################

# NOTE: Assumed to have ENSEMBL ids in the gene annotations. If not, function to
# be reviewed

complete.dea <- function(DESeqSet, #Obtain after computing DESeq2
                         comparison, 
                         lim.ma.abs= 10,  #Limit to consider in MA plot
                         fig.path= .fig)  #Path to store the MA plot figure
{
  # 1. Obtain results from DEA
  results_dea <-results(dea.dds,
                        name =  comparison,
                        alpha = 0.05)
  
  # 2. Print summary of results
  summary(results_dea)
  
  # 3. LogFC shrinkage
  res.shr <- lfcShrink(dds = dea.dds, 
                       coef=comparison,
                       res = results_dea, 
                       type = "apeglm")
  
  # 4. Store an MA-plot to see the effect
  pdf(file = file.path(.fig,paste("MAplot_",comparison,".pdf",sep="")), width=16, height=10)
  par(mfrow=c(1,2))    
  plotMA(results_dea,ylim=c(-lim.ma.abs,lim.ma.abs),main=paste("DEA ",comparison, sep=""))
  plotMA(res.shr, ylim=c(-lim.ma.abs,lim.ma.abs), main=paste("DEA ", comparison," - Shrinkage performed", sep=""))
  dev.off()

  # 5. Collect results in the same data frame: DEA resuls and logFC shrinkage
  DEA_complete <- cbind(as.data.frame(results_dea), 
                            "Shrunken_lFC" = as.data.frame(res.shr)$log2FoldChange,
                            "Shrunken_lFCSE" = as.data.frame(res.shr)$lfcSE)
  
  # Add a new column with the gene.annot
  DEA_complete$ENSEMBL <- rownames(DEA_complete)
  
  # 6. Additional annotation. REMARK: Assuming ENSEMBL ids as entries
    # SYMBOL!
  df  <- clusterProfiler::bitr(geneID = DEA_complete$ENSEMBL, fromType = "ENSEMBL", 
                                               toType = "SYMBOL",OrgDb = "org.Hs.eg.db",drop=FALSE)
    # Collapse the duplicates! 
  df <- df %>%
    dplyr::group_by(ENSEMBL) %>%
    dplyr::summarise(SYMBOL.col = paste(SYMBOL, collapse = ";"))
  
    # Put in the same order
  order <- match(DEA_complete$ENSEMBL,df$ENSEMBL)
  df <- df[order,]
    #Assign
  DEA_complete$SYMBOL <- df$SYMBOL.col
  
  # ENTREZ!
  df  <- clusterProfiler::bitr(geneID = DEA_complete$ENSEMBL, fromType = "ENSEMBL", 
                               toType = "ENTREZID",OrgDb = "org.Hs.eg.db",drop=FALSE)
  df <- df %>%
    dplyr::group_by(ENSEMBL) %>%
    dplyr::summarise(ENTREZ.col = paste(ENTREZID, collapse = ";"))
  
    # Put in the same order
  order <- match(DEA_complete$ENSEMBL,df$ENSEMBL)
  df <- df[order,]
    #Assign
  DEA_complete$ENTREZID <- df$ENTREZ.col

  # 7. Sort results by p-adj
  DEA_complete <- DEA_complete[order(DEA_complete$padj),]
  
  # 8. Return complete dataframe
  return(DEA_complete)
}
  


# Volcano plot ------------------------------------------------------------

################################################################################
## 7. Volcano plot
################################################################################

# For those enriched in DOX: darkorchid4 (WT), cadetblue4 (N108), darkorange4 (M3)
# For those enriched in CTL: mediumorchid (WT), cadetblue3 (N108), darkorange (M3)

for(i in 1:length(final.res))
{
  tab <- final.res[[i]]
  tab <- tab[-which(tab$SYMBOL == "NFKBIA"),]  #It worsens the plot
  tab <- tab[-which(tab$SYMBOL == "NA"),]
  
  # Create custom key-value shape for those outlier genes
  tab$Outlier <-ifelse(tab$padj < 1e-50, "Yes", "No")
  keyvals.shape <- ifelse(tab$Outlier == "Yes",17 , 19)
  names(keyvals.shape)[keyvals.shape == 19] <- 'adj pval > 1e-50'
  names(keyvals.shape)[keyvals.shape == 17] <- 'adj pval < 1e-50'
    # Limit to this adj pval max
  tab$padj[tab$padj < 1e-50] <- 1e-50
  
  # Create custom key-value colors: only coloured the DEGs
  tab$Sig <- "Non_sign"
  tab$Sig[which(tab$padj <0.05 & tab$Shrunken_lFC>0.2)]  <- "UP"
  tab$Sig[which(tab$padj <0.05 & tab$Shrunken_lFC< -0.2)]  <- "DOWN"
  
  if(length(grep("_WT_",names(final.res)[i]))>0)
  {
    colors <- c("darkorchid4", "mediumorchid")
  }else if(length(grep("_108_",names(final.res)[i]))>0)
  {
    colors <- c("cadetblue4", "cadetblue3")
  }else{
    colors <- c("darkorange4", "darkorange")
  }
  
  keyvals.colors <- ifelse(tab$Sig == "Non_sign","gray50" , ifelse(tab$Sig == "UP", colors[1], colors[2]))
  names(keyvals.colors)[keyvals.colors == colors[1]] <- 'UP'
  names(keyvals.colors)[keyvals.colors == colors[2]] <- 'DOWN'
  names(keyvals.colors)[keyvals.colors == "gray50"] <- 'Non_sign'
   
  p <- EnhancedVolcano(toptable = tab,
                  lab= tab$SYMBOL,
                  x ='Shrunken_lFC',
                  y= 'padj',
                  xlab = "Shrunken log2 FC",
                  ylab = "-Log10(adj pval)",
                  caption = paste0("adj pval < 0.05; | shrunken log2 FC | > 0.2"),
                  title = names(final.res)[i],
                  subtitle= "Differential Expression - INDs iKB Dataset",
                  xlim=c(-1.5,2),
                  ylim=c(0,50),
                  pCutoff = 0.05,
                  FCcutoff = 0.2,
                  legendLabels=c("Not sign.", "Shrunken lFC", "adj pval","adj pval & shrunken lFC"),
                  legendPosition = "right",
                  legendIconSize = 5.0,
                  legendLabSize = 10,
                  pointSize = 2,
                  labSize =3.5,
                  shapeCustom = keyvals.shape,
                  colCustom = keyvals.colors,
                  colAlpha = 0.75,
                  drawConnectors = TRUE,
                  widthConnectors = 0.2,
                  lengthConnectors = unit(0.005, "npc"),
                  typeConnectors = "open",
                  max.overlaps = 15,
                  border = "full",
                  gridlines.major = FALSE,
                  gridlines.minor =FALSE)
  
  pdf(file = paste0("PAPER_related/VolcanoPlot_", names(final.res)[i],".pdf"), width=10, height=8)
  print(p)
  dev.off()
}



# Venn Diagram  -----------------------------------------------------------

################################################################################
## Other: Venn Diagrams among the three comparisons
################################################################################

# Select "DEGs" based on a criteria over the adjusted p-value
pcriteria <- 0.01  #0.05
lfC_criteria = 0.0

degs.up <- lapply(final.res, function(res) {
  de <- res %>% dplyr::filter(padj < pcriteria, Shrunken_lFC > lfC_criteria, SYMBOL!="NA")
  de <- de[!is.na(de$SYMBOL),"SYMBOL"]
  })
names(degs.up) <- paste("UP", names(degs.up), sep="_")

degs.down <- lapply(final.res, function(res) {
  de <- res %>% dplyr::filter(padj < pcriteria, Shrunken_lFC < -(lfC_criteria), SYMBOL!="NA")
  de <- de[!is.na(de$SYMBOL),"SYMBOL"]
})
names(degs.down) <- paste("DOWN", names(degs.down), sep="_") 


#Upset plot: we need a list o ofnamed vectors

toupset <- c(degs.up[c(1,3)], degs.down[c(1,3)])
#names(toupset) #"UP_DEA_108_D_vs_CT"   "UP_DEA_WT_D_vs_CT"    "DOWN_DEA_108_D_vs_CT" "DOWN_DEA_WT_D_vs_CT" 

names(toupset) <- c("Up_Dox_N108","Up_Dox_WT","Down_Dox_N108","Down_Dox_WT") 

pdf(file = "PAPER_related/UPSET_DEGs_Up_Down_N108_WT_padj0.05_slFC_0.2.pdf", width=5, height=4)
upset(fromList(toupset),
      sets = c("Down_Dox_N108", "Down_Dox_WT", "Up_Dox_N108","Up_Dox_WT"),
      keep.order=TRUE,
      mainbar.y.label = "Unique/Common DEGs",
      sets.x.label = "DEGs number",
      point.size = 2,
      line.size = 0.5,
      text.scale = 1,
      set_size.show = TRUE,
      set_size.scale_max= 400,
      main.bar.color = c("cadetblue3", "cadetblue4","mediumorchid", "darkorchid4","gray60", "gray60","gray60"),
      sets.bar.color = c("cadetblue3", "mediumorchid","cadetblue4", "darkorchid4"),
      nsets = 4)
dev.off()

# In both cases (both filtering criteria), there are 4 genes in common between Up_Dox_N108
# and Donw_Dox_WT which seems weird. THose genes are:
intersect(toupset$Up_Dox_N108, toupset$Down_Dox_WT)
#"LAMC2" "ITGAV" "OPTN"  "NFKB2"

