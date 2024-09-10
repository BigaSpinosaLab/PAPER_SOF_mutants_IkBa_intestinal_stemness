################################################################################
##        Title : RNAseq data analysis - Induced ikB
##  Description : This script is for GSEA Pre ranked analysis for particular gene
##                signature: ISC signature Munoz et al. EMBOJ 2012
##   Researcher : D.Alvarez-Villanueva
##       Author : María Maqueda
##         Date : 27th July
################################################################################

# Set up paths ------------------------------------------------------------

################################################################################
## 1. Set up paths
################################################################################

# Load packages -----------------------------------------------------------

################################################################################
## 2. Load packages
################################################################################

require(fgsea)
require(org.Hs.eg.db)
require(openxlsx)
require(dplyr)
require(grid)  # For adding some annotations in the plot

# Import data -----------------------------------------------------------

################################################################################
## 3. Import data: complete list of genes 
################################################################################

# Object: final.res
load("~/Documents/Projects/Daniel_Alvarez/Ikba_mutants/RNAseq_INDs_IkB_Dani/RData/DEA_RESULTS_INDs_ikB.RData")

# Rank genes based on log2 FC shrunken and using gene symbols
ranked.genes <- lapply(final.res, function(case) {
  v <- case$Shrunken_lFC
  names(v) <- case$SYMBOL
  v <- sort(v, decreasing=TRUE)
  return(v)
})

# Genes with no Symbol annotation  out of the initial 22,993 genes
lapply(ranked.genes, function(ranking) length(which(names(ranking) == "NA")))
# [1] 5,898

# Remove those ranked genes without any name
ranked.genes.red <- lapply(ranked.genes, function(ranking) ranking[-which(names(ranking) == "NA")])
# 17,095

# Remove duplicated gene names. Check that these are exceptionally (too few values)
lapply(ranked.genes.red, function(ranking) length(which(duplicated(names(ranking)))))  # 32 out of 17,095
ranked.genes.red <- lapply(ranked.genes.red, function(ranking) ranking[-which(duplicated(names(ranking)))])

################################################################################
## 4. Import data: Gene sets to be used
################################################################################

#===============
# Muñoz signature: Muñoz et al. EMBO J 2012
#===============

ISC_Munoz <- readRDS("~/Documents/Projects/Daniel_Alvarez/Ikba_mutants/ISC_signatures_REF/ISC_Muñoz_genes.RDS")

################################################################################
## Compute Pre-Ranked GSEA for custom gene sets
################################################################################

GSEA.res <- lapply(ranked.genes.red, function(ranking)
  {
  set.seed(123)
  fgseaRes <- fgseaMultilevel(pathways = list("ISC_Munoz" = ISC_Merlos), 
                              stats = ranking,
                              scoreType = "std",
                              nPermSimple = 10000)
  # Fix leading edge column
  fgseaRes <- fgseaRes %>% 
    mutate(leadingEdge = sapply(leadingEdge, toString))
  
  return(fgseaRes)
})
names(GSEA.res) <- names(ranked.genes.red)

# LeadingEdge column: Transform entrez ids into symbols
# NOT REQUIRED since symbols provided
# fgseaRes[, leadingEdge2 := mapIdsList(
#       x=org.Hs.eg.db, 
#       keys=leadingEdge,
#       keytype="ENTREZID", 
#       column="SYMBOL")]
    
## Store results in an excel

tostore <- do.call(rbind, GSEA.res)
tostore$CONDITION <- names(ranked.genes.red)
tostore <- tostore[,c(9,1:8)]

hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "#4F80BD")
write.xlsx(x = tostore, 
           file= "PAPER_related/GSEA_ISC_Signatures_MERLOS-SUAREZ_N108_M3_WT.xlsx",
             headerStyle=hs)
  
################################################################################
## Enrichment plot with results
################################################################################

#ISC Munoz signature
param <- c("DEA_108_D_vs_CT", "cadetblue4",0.70)
param <- c("DEA_M3_D_vs_CT", "darkorange3",0.55)
param <- c("DEA_WT_D_vs_CT", "darkorchid4",0.65)


pdf(file = paste0("PAPER_related/GSEA_ISC_Munoz_",param[1],".pdf"), width=4, height=4)

plotEnrichment_adapted(pathway = ISC_Merlos, 
                       stats = ranked.genes.red[[param[1]]],
                       col_walking = param[2],
                       fgsea_res = tostore[tostore$CONDITION == param[1],],
                       ypos=param[3]) + 
  labs(title=paste0("ISC Signature (Muñoz et al.)\n",param[1])) 
dev.off()

################################################################################
## ANNEX. Enrichment plot adapted (aesthetics) from the original from fgsea package
################################################################################

require(data.table)  

plotEnrichmentData <- function(pathway, stats,
                                 gseaParam=1) {
    
    if (any(!is.finite(stats))){
      stop("Not all stats values are finite numbers")
    }
    
    rnk <- rank(-stats)
    ord <- order(rnk)
    
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
    
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    pathway <- unique(pathway)
    
    gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats = pathway,
                            returnAllExtremes = TRUE)
    
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.table(rank=c(0, xs, n + 1), ES=c(0, ys, 0))
    ticks <- data.table(rank=pathway, stat=statsAdj[pathway])
    stats <- data.table(rank=seq_along(stats), stat=statsAdj)
    
    res <- list(
      curve=toPlot,
      ticks=ticks,
      stats=stats,
      posES=max(tops),
      negES=min(bottoms),
      spreadES=max(tops)-min(bottoms),
      maxAbsStat=max(abs(statsAdj)))
  }
  

plotEnrichment_adapted <- function(pathway, 
                                   stats,
                                   fgsea_res,
                                   col_walking,
                                   ypos,
                             gseaParam=1,
                             ticksSize=0.2) {
    
    pd <- plotEnrichmentData(
      pathway = pathway,
      stats = stats,
      gseaParam = gseaParam)
    
    # Create the text to be added in the plot regarding NES and padj
    nes = round(fgsea_res$NES, digits=1)
    padj = formatC(fgsea_res$pval, digits=1, format="e")  # Let's take the pval since we have only tested two gene sets
    # It does not make sense to make FDR in this scenario
    
    grob <- grobTree(textGrob(paste("NES=", nes,"\npval=",padj,sep=""), x=0.1,  y=0.15, hjust=0,
                              gp=gpar(col="gray43", 
                                      fontsize=14, 
                                      fontface="bold")))

    # Create the text to include KO and HE-WT labels
    grob2 <- grobTree(textGrob("DOX", x=0.1,  y=ypos, hjust=0,
                              gp=gpar(col="darkred", 
                                      fontsize=15, 
                                      fontface="bold")))
    
    grob3 <- grobTree(textGrob("CTL", x=0.8,  y=ypos, hjust=0,
                               gp=gpar(col="khaki4", 
                                       fontsize=15, 
                                       fontface="bold")))
    
    
    with(pd,
         ggplot(data=curve) +
           geom_line(aes(x=rank, y=ES), color=col_walking) +  
           geom_segment(data=ticks,
                        mapping=aes(x=rank, y=-spreadES/16,
                                    xend=rank, yend=spreadES/16),
                        linewidth=ticksSize) +
           annotation_custom(grob) +
           annotation_custom(grob2) +
           annotation_custom(grob3) +
           geom_hline(yintercept=posES, colour="red", linetype="dashed") +
           geom_hline(yintercept=negES, colour="red", linetype="dashed") +
           geom_hline(yintercept=0, colour="black") +
           theme(
             #panel.border=element_rect(color="black"),
             panel.background = element_blank(),
             panel.grid.major =element_blank()
             #panel.grid.major=element_line(color="grey92")
           ) +
           labs(x="Rank", y="Enrichment Score"))
  }
  

