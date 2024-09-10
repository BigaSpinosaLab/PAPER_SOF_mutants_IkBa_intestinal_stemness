################################################################################
##        Title : IkBa SOF mutants - ChIPseq IkBa FLAG over CRC cell lines 
##                (HCT and HT29)
##  Description : This script is for annotating genes to identified ChIPseq peaks
##                and create a Consensus Peakset
##   Researcher : Dani Álvarez 
##       Author : María Maqueda
##         Date : 20th March 2024
################################################################################

# Objective: 
# a) Annotate genes to identified peaks in HCT and HT29 cell lines
# b) Create a consensus peak set per cell line
# c) Assess the peaks overlapping between both cell lines

# NOTE: Peaks were called with MACS2 (adj pval < 0.05 applied) since IkBa pattern
# is assumed to have a narrow profile

################################################################################
## 0. Load packages
################################################################################

require(ChIPpeakAnno)
library(rtracklayer)
require(openxlsx)
require(dplyr)

# Required for annotation
require(ChIPseeker)  
require(GenomicFeatures)  # For creating a TxDb object from gtf file
require(org.Hs.eg.db)

################################################################################
## 1. Import files with called peaks
################################################################################

results.path <- "/Volumes/grcmc/ABigas_Bioinfo/PROJECTS/Ikba_mutants_DANI/ChIP_IkBa_Flag_CRC_cell_lines/Peaks/source"

files <- list("HCT" = list.files(path = results.path,full.names = TRUE, pattern = "HCT"),
              "HT" = list.files(path = results.path,full.names = TRUE, pattern = "HT"))

# Peak files
# ---------------

peak.files <- lapply(files, function(f) {
  lapply(f, function(ff) {
    read.table(file = ff, header=FALSE,comment.char = "")})
  })
  
names(peak.files$HCT) <- c("HCT_1", "HCT_2", "HCT_3")
names(peak.files$HT) <- c("HT_1", "HT_2", "HT_3")

peak.files <- lapply(peak.files, function(cell_line)
{
  lapply(cell_line, function(set) {
    colnames(set) <- c("Chrom", "Start", "End",  "Peak_name","Score", "Strand",
                       "FC.summit", "log10pval", "log10qval", "Relative.summit.pos")
    return(set)
  }) })

################################################################################
## 2. Number of peaks per sample
## Number of peaks FDR adj pval 5% (this is the default criteria in previous files)
################################################################################

# Number of peaks
lapply(peak.files, function(cl) {lapply(cl, nrow)})

# $HCT
# $HCT$HCT_1
# [1] 4349
# 
# $HCT$HCT_2
# [1] 4973
# 
# $HCT$HCT_3
# [1] 4639
# 
# 
# $HT
# $HT$HT_1
# [1] 1178
# 
# $HT$HT_2
# [1] 712
# 
# $HT$HT_3
# [1] 1063


# (Optional) Make a more restrictive selection for downstream 
# threshold = 5
# lapply(peak.files, function(l) {
#   lapply(l, function(s) {
#     a <- s %>% filter(FC.summit > threshold)
#     nrow(a)
#   })
# })

# $HCT$HCT_1
# [1] 2327
# 
# $HCT$HCT_2
# [1] 2649
# 
# $HCT$HCT_3
# [1] 2472
# 
# 
# $HT
# $HT$HT_1
# [1] 800
# 
# $HT$HT_2
# [1] 545
# 
# $HT$HT_3
# [1] 671


################################################################################
## 3. Create GRanges objects from those peak files
################################################################################

# FIRST, peaks list must be converted GRanges objects
p.gr <- lapply(peak.files, function(cl)
{
  lapply(cl, function(res) {
    makeGRangesFromDataFrame(df = res,
                             keep.extra.columns=TRUE,
                             seqnames.field="Chrom",
                             start.field="Start",
                             end.field="End",
                             strand.field = "Strand",
                             starts.in.df.are.0based=TRUE)
    
  })
})
  

p.gr <- lapply(p.gr, function(pp) GRangesList(pp))

# Let's add names to the peaks properly
# HCT cell line
samples <- names(p.gr$HCT)
p.gr$HCT <-GRangesList(lapply(samples, function(regions) 
  {
  names(p.gr$HCT[[regions]]) <- p.gr$HCT[[regions]]$Peak_name
  return(p.gr$HCT[[regions]])
}))
names(p.gr$HCT) <- samples

# HT cell line
samples <- names(p.gr$HT)
p.gr$HT <-GRangesList(lapply(samples, function(regions) 
{
  names(p.gr$HT[[regions]]) <- p.gr$HT[[regions]]$Peak_name
  return(p.gr$HT[[regions]])
}))
names(p.gr$HT) <- samples

################################################################################
## 4. Find overlapping peaks (by at least 1bp). These will be merged
################################################################################

# Overlapping peaks
overlapping <- findOverlapsOfPeaks(p.gr$HT,
                                    ignore.strand = TRUE, 
                                    connectedPeaks = "keepAll", 
                                    maxgap = -1)  # Maxgap default value: 1base in common

# Results: venn counts
# overlapping$venn_cnt

write.table(x = as.table(overlapping$venn_cnt), 
            file = "results/Consensus_PeakSet/other/Venn_Counts_HT29_replicates.txt",
            row.names = FALSE, quote = FALSE)

# Store the overlapping as as BED file 
#rtracklayer::export.bed(overlapping$mergedPeaks,con='results/Consensus_PeakSet/BED/ChIP_IkBa_ConsensusPeakset_HCT_all_chr.bed')
rtracklayer::export.bed(overlapping$mergedPeaks,con='results/Consensus_PeakSet/BED/ChIP_IkBa_ConsensusPeakset_HT29_all_chr.bed')


#==================================================
# Final HCT/HT29 consensus without the non-canonical chromosomes
#consensus.HCT <- keepSeqlevels(overlapping$mergedPeaks, c(1:22, "X"), pruning.mode="coarse")  # Optional: remove non-canonical chrom (few cases)
consensus.HT29 <- keepSeqlevels(overlapping$mergedPeaks, c(1:22, "X"), pruning.mode="coarse")  # Optional: remove non-canonical chrom (few cases)

# Include in the overlapping peaks, the maximum FC.summit among the replicates (with overlapped peaks)
info <- overlapping$peaksInMergedPeaks

for(i in 1:length(consensus.HT29$peakNames))
{
  peaks_reps <- consensus.HT29$peakNames[[i]]
  fc.max <-  max(as.data.frame(info[which(names(info) %in% peaks_reps),])$FC.summit)
  consensus.HT29$FC.summit.max[i] = fc.max
}

rtracklayer::export.bed(consensus.HT29,con='results/Consensus_PeakSet/BED/ChIP_IkBa_ConsensusPeakset_HT29_canonical_chr.bed')

################################################################################
## 5. Create txdb for Annotate the consensus PeakSet
################################################################################

# Create the txdb according to the same used during alignment
metadata <- data.frame(name="Resource URL",
                       value=paste0("ftp://ftp.ensemblgenomes.org/pub/","release-102/gtf/mus_musculus/"))

txdb <- makeTxDbFromGFF(file = "/Volumes/projectscomput/cancer/db_files/Genomes/Ensembl/human/hg38/release-106/Homo_sapiens.GRCh38.106.gtf",
                        format="gtf",
                        dataSource="Ensembl_FTP_repository",
                        organism="Homo Sapiens",
                        taxonomyId=9606,
                        circ_seqs=NULL,
                        chrominfo=NULL,
                        miRBaseBuild=NA,
                        metadata=metadata,
                        dbxrefTag="gene_id")


################################################################################
## 6. Annotate the consensus PeakSet and the individual peaks
################################################################################

# Consensus
Consensus <- annotatePeak(peak = consensus.HT29, 
                 tssRegion = c(-5000,100), 
                 TxDb = txdb,
                 genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                               "Downstream", "Intergenic"),  # Default genomic Annotation
                 annoDb = "org.Hs.eg.db")

# Alternative, if some categories want to be aggregated (i.e. all Promoter cases)

stat <- Consensus@annoStat

stat.summary <- data.frame("Feature"= c("Promoter", "Other", "Intron", "Distal Intergenic"),
                           "Frequency" = c(sum(stat$Frequency[grep("Promoter", stat$Feature)]),
                                           sum(stat$Frequency[grep("UTR|Downsream|Exon", stat$Feature)]),
                                           sum(stat$Frequency[grep("Intron", stat$Feature)]),
                                           sum(stat$Frequency[grep("Distal", stat$Feature)]))
                           )

pdf(file = "results/Consensus_PeakSet/ChIP_IkBa_Annotation_Distribution_HT29_Consensus_SUMMARY.pdf", width=6, heigh =4)
pie(stat.summary$Frequency, 
    labels = paste0(stat.summary$Feature,"\n", round(stat.summary$Frequency,digits=1), "%"),
    col = c("lightblue", "grey50","orchid", "gold2"),
    main="IkBa consensus peaks - HCT cell line")
dev.off()

# Add annotations 
peakAnno <- lapply(p.gr$HT, function(set){
  p <-annotatePeak(peak = set, 
               tssRegion = c(-5000,100), 
               TxDb = txdb,
               genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                                         "Downstream", "Intergenic"),  # Default genomic Annotation
              annoDb = "org.Hs.eg.db")
  #p <- cbind(names(set), as.data.frame(p))
  return(as.data.frame(p))
})


# Store annotations: Consensus + Individual peaks
tostore <- peakAnno
tostore$Consensus <- as.data.frame(Consensus)
                
hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "#4F80BD")
write.xlsx(x= tostore, file = file.path("results/Consensus_PeakSet/ChIP_IkBa_Consensus_PeakSet_HT29.xlsx"), headerStyle=hs)

################################################################################
## 6. ANNEX: Consensus of consensus
################################################################################

overlap_of_overlap <- findOverlapsOfPeaks(consensus.HT29, consensus.HCT,
                                          ignore.strand = TRUE, 
                                          connectedPeaks = "keepAll", 
                                          maxgap = -1)  # Maxgap default value: 1base in common

# Most of the peaks are in common: 524 out of the 610 from HT29. Let's store the HT29 exclusive

exclusive.HT29 <- annotatePeak(peak = overlap_of_overlap$uniquePeaks[grep("consensus.HT29", names(overlap_of_overlap$uniquePeaks)),], 
                          tssRegion = c(-5000,100), 
                          TxDb = txdb,
                          genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                                        "Downstream", "Intergenic"),  # Default genomic Annotation
                          annoDb = "org.Hs.eg.db")
write.xlsx(x= exclusive.HT29, file = file.path("results/Consensus_PeakSet/ChIP_IkBa_Exclusive_HT29_against_HCT.xlsx"), headerStyle=hs)


rtracklayer::export.bed(overlap_of_overlap$uniquePeaks[grep("consensus.HT29", names(overlap_of_overlap$uniquePeaks)),],
                        con='results/Consensus_PeakSet/BED/ChIP_IkBa_Exclusive_HT29_against_HCT.bed')
