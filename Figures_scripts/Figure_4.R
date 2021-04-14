library(GenomicRanges)
library(RepViz)
library(ChIPpeakAnno)

###################################################
##
## Figure 1 PRKCI chr3:169940220-170023770 H3K27Ac
##
###################################################

region<- GRanges("chr14:25102000-25105000")
jpeg("Supplementary_3C_chr14:25102000-25105000.jpeg", width = 1400, height = 1000)
RepViz(region = region,
           genome = "hg19",
           BAM = "BAM_input.csv",
           BED = "BED_input.csv",
           avgTrack = T,
           geneTrack = T,
           verbose=T,
       col = c("orange","red","green","violet","blue","black","black"),cex = 2)
dev.off()


region<- GRanges("chr1:25390500-25393650")
jpeg("supplementary_3B_chr1:25390500-25393650.jpeg", width = 1400, height = 1000)
RepViz(region = region,
       genome = "hg19",
       BAM = "BAM_input.csv",
       BED = "BED_input.csv",
       avgTrack = T,
       geneTrack = T,
       verbose=T,
       col = c("orange","red","green","violet","blue","black","black"),cex = 2)
dev.off()


ROTS <- read.table("ROTS_l1pf.bed",sep="\t",header=F)
DB <- read.table("DiffBind.bed",sep="\t",header=F)
#DB_edger <- read.table("../Software_comparison/DBsitesH3K36me3_edger.bed",sep="\t",header=F)
DR  <- read.table("diffReps.bed",sep="\t",header=F)
PePr <- read.table("PePr.bed",sep="\t",header=F)
THOR <- read.table("THOR.bed")


THOR <- THOR[order(THOR$V7, decreasing = F),]
THOR$rank <- 1:dim(THOR)[1]
PePr <- PePr[order(PePr$V8, decreasing = F),]
PePr$rank <- 1:dim(PePr)[1]

ROTS <- ROTS[which(ROTS$V8 < 0.05),]
DB <- DB[which(DB$V7 < 0.05),]
DR <- DR[which(DR$V7 < 0.05),]
THOR <- THOR[which(THOR$V7 < 0.05),]
PePr <- PePr[which(PePr$V8 < 0.1),]
#################################################################################
##
##
##        Main
##
##
#################################################################################

## Genome





ROTS <- ROTS[1:2000,c(1:3,9)]

DB <- DB[1:2000,c(1:3,8)]

DR <- DR[1:2000,c(1:3,8)]

PePr <- PePr[1:2000,c(1:3,9)]

THOR <- THOR[1:2000,c(1:3,8)]




colnames(ROTS) <- c("space","start","end","fc")
colnames(DB) <- c("space","start","end","fc")
colnames(DR) <- c("space","start","end","fc")
colnames(PePr) <- c("space","start","end","fc")
colnames(THOR) <- c("space","start","end","fc")


## Conversion of the Data frames to GRanges objects

ROTS <- toGRanges(ROTS)
DB <- toGRanges(DB)
DR <- toGRanges(DR)
PePr <- toGRanges(PePr)
THOR <- toGRanges(THOR)