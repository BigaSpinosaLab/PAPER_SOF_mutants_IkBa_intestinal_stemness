################################################################################
##        Title : RNAseq data analysis - mESCs and EBs (48h and 96h) ikba KO
##  Description : This script is for conducting ORA against ChEA
##   Researcher : Luis Galán/D.Álvarez
##         Date : 23rd July 2022
################################################################################

# ONLY ENCODE and ChEA Consensus TFs from ChIP-X is used
# https://maayanlab.cloud/Enrichr/#libraries

# With ChEA 2022 we do not get the same results

# Set up paths ------------------------------------------------------------

################################################################################
## 1. Set up paths
################################################################################


# Load packages -----------------------------------------------------------

################################################################################
## 2. Load packages
################################################################################

require(clusterProfiler)
require(org.Hs.eg.db)
require(openxlsx)
require(enrichplot)
require(GSA) # required for reading the .gmt file (ChEA database)

# Import data -------------------------------------------------------------

################################################################################
## 3. Import data: DEGs
################################################################################

# All results  (final.res object)
load("~/Documents/Projects/Daniel_Alvarez/Ikba_mutants/RNAseq_INDs_IkB_Dani/RData/DEA_RESULTS_INDs_ikB.RData")

# Select DEGs in N108 - Down regulated
pcriteria <- 0.05  #0.05
lfC_criteria = 0.2

degs.down <- lapply(final.res, function(res) {
  de <- res %>% dplyr::filter(padj < pcriteria, Shrunken_lFC < -(lfC_criteria), SYMBOL!="NA")
  
  # OPTIONAL: order them by log FC and select the top 50
  #de <- de[order(de$padj, decreasing=FALSE),]
  #de <- de[1:100,]
  
  # Get only the symbol
  symbol <- de[!is.na(de$SYMBOL),"SYMBOL"]
  return(symbol)
})
names(degs.down) <- paste("DOWN", names(degs.down), sep="_") 


# 
# Downloaded from: https://maayanlab.cloud/Enrichr/#libraries
chea.encode <- GSA::GSA.read.gmt(filename = "data/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.txt")

################################################################################
## 4. Prepare ChEA database for the functional analysis
################################################################################

chea.encode <- chea.encode[-3]   # There're no descriptions to deal with 

# In all gene sets there is an empty element. Let's remove it
chea.encode$genesets <- lapply(chea.encode$genesets, function(gs) gs[-length(gs)])

# For enriching ChEA: a df must be prepared following instructions in
# http://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/

lengths <- (unlist(lapply(chea.encode$genesets, function(gs) length(gs))))

chea.encode.ora <- data.frame(
  "tf.id" = unlist(sapply(seq(1,length(chea.encode$geneset.names)), function(r) 
    rep(chea.encode$geneset.names[r], each = lengths[r]))),
  "genes" = unlist(chea.encode$genesets))


################################################################################
## 5. ChEA - ENCODE  overrepresentation analysis
################################################################################

ora.ChEA <- lapply(degs.down[c(1,3)], function(g){
  enricher(gene   = g,
           pvalueCutoff = 1, # To explore a longer list
           qvalueCutoff = 1, # To explore a longer list
           pAdjustMethod = "BH",
           minGSSize = 10,
           maxGSSize = 5000, # Largest GS is 4643
           TERM2GENE = chea.encode.ora)
})

################################################################################
## 6. ChEA - ENCODE  ORA plot
################################################################################

# Only significant TFs are shown (adjusted pval < 0.05)


results = ora.ChEA

# Dotplot (from clusterProfiler) of the subset scenario
plots <- lapply(names(results), function(scenario){
  # Results are already sorted per adjusted pvalue
  long = length(which(as.data.frame(results[[scenario]])$p.adjust < 0.05))
  dotplot(results[[scenario]], showCategory = long, #ifelse(long >10, 10, long),
          title=scenario) +
    scale_colour_gradient(limits=c(0, 0.05), low="red") +
    theme_bw(base_size = 10) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.y = element_text(size=10, face="bold"),
          axis.text.x = element_text(size=10, face="bold")) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 24)) +
    scale_x_continuous(expand = expansion(mult=0, add=0.025)) 
})

pdf(file = "PAPER_related/ORA_ChEA_ENCODE_Down_DEGs_108.pdf",
    width=4, height=4)
print(plots[[1]])
dev.off()

pdf(file = "PAPER_related/ORA_ChEA_ENCODE_Down_DEGs_WT.pdf",
    width=4, height=4)
print(plots[[2]])
dev.off()



################################################################################
## FINAL. Store all results in an excel file and RDS file
# We will store an excel file per timepoint
################################################################################

hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "#4F80BD")

write.xlsx(x = ora.ChEA, 
           file= "results/ORA_results_ChEA_ENCODE_N108_WT.xlsx", 
           headerStyle=hs)





