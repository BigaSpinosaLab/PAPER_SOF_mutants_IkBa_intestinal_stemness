################################################################################
##        Title : RNAseq data analysis - Induced ikB
##  Description : This script is for GSEA Pre ranked analysis
##   Researcher : D. Alvarez- Villanueva
##         Date : 25th May 2022
################################################################################

# Set up paths ------------------------------------------------------------

################################################################################
## 1. Set up paths
################################################################################

.wd <- getwd()
.res <- file.path(.wd, "results")
.dat <- file.path(.wd,"data")
.tab <- file.path(.res,"tables")
.fig <- file.path(.res,"figures")
.RData <- file.path(.wd,"RData")

# Load packages -----------------------------------------------------------

################################################################################
## 2. Load packages
################################################################################

require(fgsea)
require(gage)  # To be used for retrieving the KEGG and GO gene sets
require(org.Hs.eg.db)
require(openxlsx)
require(msigdbr)

# Import data -----------------------------------------------------------

################################################################################
## 3. Import data: complete list of genes 
################################################################################

# Import the complete lists of genes (three comparisons)
load(file = file.path(.RData,"DEA_RESULTS_INDs_ikB.RData"))

# Rank genes based on log2 FC shrunken and using Entrez ID names

ranked.genes <- lapply(final.res, function(genes.list){
  rank <- genes.list$Shrunken_lFC
  names(rank) <- genes.list$ENTREZID
  rank <- sort(rank,decreasing=TRUE)
})

# Genes with no Entrez annotation  out of the initial 22993
# > lapply(ranked.genes, function(r) length(which(names(r) %in% "NA")))
# $DEA_108_D_vs_CT
# [1] 5898
# 
# $DEA_M3_D_vs_CT
# [1] 5898
# 
# $DEA_WT_D_vs_CT
# [1] 5898

# Remove those ranked genes without any name
ranked.genes.red <- lapply(ranked.genes, function(r) r[-which(names(r) %in% "NA")])
  

################################################################################
## 4. Gene sets to be used
################################################################################

# MSigdB HallMark Gene sets
# Hallmark collection (others can be imported also from MSigdb)
# Check them: https://www.gsea-msigdb.org/gsea/msigdb/
msigdb_hallmark <- msigdbr(species = "Homo sapiens", category = "H") 
# Put them in the required format
msigdbr_hallmark.gs = split(x = msigdb_hallmark$entrez_gene, f = msigdb_hallmark$gs_name)

# Store gene sets from different databases into a list
databases.gs <- list("Hallmark" = msigdbr_hallmark.gs)

################################################################################
## 5. Compute Pre-Ranked GSEA for KEGG gene sets
################################################################################

# Define the list of genes to be analysed and run the rest of the code
#ranking <- ranked.genes.red[["DEA_108_D_vs_CT"]]
#ranking <- ranked.genes.red[["DEA_M3_D_vs_CT"]]
ranking <- ranked.genes.red[["DEA_WT_D_vs_CT"]]

# Compute fgsea for each database
# About score type check: https://github.com/ctlab/fgsea/issues/87

fgseaRes.all <- lapply(databases.gs, function(database){
  set.seed(123)
  fgseaRes <- fgseaMultilevel(pathways = database, 
                              stats = ranking,
                              scoreType = "std",
                              minSize=10,
                              maxSize=500)
  
  # LeadingEdge column: Transform entrez ids into symbols
  fgseaRes[, leadingEdge2 := mapIdsList(
    x=org.Hs.eg.db, 
    keys=leadingEdge,
    keytype="ENTREZID", 
    column="SYMBOL")]
  
  # Return results
  return(fgseaRes)
})

# Example of plotting just one pathway
plotEnrichment(databases.gs$Hallmark[["HALLMARK_MYC_TARGETS_V1"]],
               ranking) + labs(title="HALLMARK_MYC_TARGETS_V1")

################################################################################
## 6. Store results in an excel and RData file
################################################################################

fgseaRes.all <- lapply(fgseaRes.all, function(database){
  database <- database %>% 
    mutate(leadingEdge = sapply(leadingEdge, toString),
           leadingEdge2 = sapply(leadingEdge2, toString))
})
 
hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "#4F80BD")

write.xlsx(x = fgseaRes.all, file= file.path(.res,
                                             "GSEA_results_WT_D_vs_CT.xlsx"), # Remember to specify the proper name!!!
           headerStyle=hs)

# Store results in a RData
save(fgseaRes.all, file = file.path(.RData, "GSEA_Results_WT_D_vs_CT.RData")) # Remember to specify the proper name!!!
