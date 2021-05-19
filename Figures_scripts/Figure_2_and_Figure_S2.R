library(GenomicRanges)
library(IRanges)
library(corrplot)
library(ChIPpeakAnno)

## YF_ATAC files
YF_ROTS <- read.table("source/YF_ATAC/ROTS_YF_narrow_pooledpeaks.bed",sep=" ",header=T)
YF_DB <- read.table("source/YF_ATAC/DBsites_YF_deseq2_narrow_pooledpeaks.bed",sep="\t",header=F)
YF_DR  <- read.table("source/YF_ATAC/diffReps.bed",sep="\t",header=F)
YF_PePr <- read.table("source/YF_ATAC/PePr2.bed",sep="\t",header=F)
YF_THOR <- read.table("source/YF_ATAC/THOR.bed")
YF_MAnorm <- read.table("source/YF_ATAC/MAnorm2_YF_pooledpeaks.bed")

## IFN_ATAC files
IFN_ROTS <- read.table("source/IFN_ATAC/ROTS_IFN_narrow_pooledpeaks.bed",sep=" ",header=T)
IFN_DB <- read.table("source/IFN_ATAC/DBsites_IFN_deseq2_narrow_pooledpeaks.bed",sep="\t",header=F)
IFN_DR  <- read.table("source/IFN_ATAC/diffReps.bed",sep="\t",header=F)
IFN_PePr <- read.table("source/IFN_ATAC/PePr.bed",sep="\t",header=F)
IFN_THOR <- read.table("source/IFN_ATAC/THOR.bed")
IFN_MAnorm <- read.table("source/IFN_ATAC/MAnorm2_IFN_pooledpeaks.bed")

## RA H3K4me3 files
K4_ROTS <- read.table("source/RA_H3K4me3/ROTS_K4_narrow_pooledpeaks.bed",sep=" ",header=T)
K4_DB <- read.table("source/RA_H3K4me3/DBsites_H3K4me3_deseq2_narrow_pooledpeaks.bed",sep="\t",header=F)
K4_DR  <- read.table("source/RA_H3K4me3/diffReps.bed",sep="\t",header=F)
K4_PePr <- read.table("source/RA_H3K4me3/PePr.bed",sep="\t",header=F)
K4_THOR <- read.table("source/RA_H3K4me3/THOR.bed")
K4_MAnorm <- read.table("source/RA_H3K4me3/MAnorm2_K4_pooledpeaks.bed")

## RA H3K36me3 files
K36_ROTS <- read.table("source/RA_H3K36me3/ROTS_K36_pooledpeaks.bed",sep=" ",header=T)
K36_DB <- read.table("source/RA_H3K36me3/DBsitesH3K36me3_deseq2_pooledpeaks.bed",sep="\t",header=F)
K36_DR  <- read.table("source/RA_H3K36me3/diffReps.bed",sep="\t",header=F)
K36_PePr <- read.table("source/RA_H3K36me3/PePr.bed",sep="\t",header=F)
K36_THOR <- read.table("source/RA_H3K36me3/THOR.bed")
K36_MAnorm <- read.table("source/RA_H3K36me3/MAnorm2_K36_pooledpeaks.bed")

# switch 1 for 2000 first switch 2 for pval < 0.05
main_function <- function(DB, ROTS, MAnorm, DR, PePr, THOR, name="generic_name",switch){
  ## Sorting peaks and ranking them if needed
  THOR <- THOR[order(THOR$V7, decreasing = F),]
  THOR$rank <- 1:dim(THOR)[1]
  PePr <- PePr[order(PePr$V8, decreasing = F),]
  PePr$rank <- 1:dim(PePr)[1]
  DR <- DR[order(DR$V8, decreasing = F),]
  DR$rank <- 1:dim(DR)[1]
  MAnorm <- MAnorm[order(MAnorm$V7, decreasing = F),]
  MAnorm$rank <- 1:dim(MAnorm)[1]
  
  ## Choose according to purpose 2000 first or < 0.05
  ## Select the number of peaks to include in the figure
  if(switch == 1){
    ROTS <- ROTS[1:2000,c(2:4,12)]
    DB <- DB[1:2000,c(1:3,6)]
    MAnorm <- MAnorm[1:2000,c(1:3,9)]
    DR <- DR[1:2000,c(1:3,8)]
    PePr <- PePr[1:2000,c(1:3,9)]
    THOR <- THOR[1:2000,c(1:3,8)]
  }
  if(switch == 2){
    ROTS <- ROTS[which(ROTS$FDR < 0.05),c(2:4,12)]
    DB <- DB[which(DB$V5 < 0.05),c(1:3,6)]
    MAnorm <- MAnorm[which(MAnorm$V7 < 0.05),c(1:3,9)]
    DR <- DR[which(DR$V7 < 0.05),c(1:3,8)]
    PePr <- PePr[which(PePr$V8 < 0.05),c(1:3,9)]
    THOR <- THOR[which(THOR$V7 < 0.05),c(1:3,8)]
  }
  

  
  ## Rename the columns of the different data frames to be
  ## able to convert them to GRanges
  colnames(ROTS) <- c("space","start","end","fc")
  colnames(DB) <- c("space","start","end","fc")
  colnames(MAnorm) <- c("space","start","end","fc")
  colnames(DR) <- c("space","start","end","fc")
  colnames(PePr) <- c("space","start","end","fc")
  colnames(THOR) <- c("space","start","end","fc")
  
  
  ## Conversion of the Data frames to GRanges objects with ChIPpeakAnno
  ROTS <- toGRanges(ROTS)
  DB <- toGRanges(DB)
  MAnorm <- toGRanges(MAnorm)
  DR <- toGRanges(DR)
  PePr <- toGRanges(PePr)
  THOR <- toGRanges(THOR)
  
  #create the return matrix
  my_matrix <- matrix(0,6,6)
  if(switch == 1){
    my_matrix[lower.tri(my_matrix)] <- c(min(sum(countOverlaps(ROTS,DB)>0),sum(countOverlaps(DB,ROTS)>0))/2000,
                                         min(sum(countOverlaps(ROTS,MAnorm)>0),sum(countOverlaps(MAnorm,ROTS)>0))/2000,
                                         min(sum(countOverlaps(ROTS,DR)>0),sum(countOverlaps(DR,ROTS)>0))/2000,
                                         min(sum(countOverlaps(ROTS,PePr)>0),sum(countOverlaps(PePr,ROTS)>0))/2000,
                                         min(sum(countOverlaps(ROTS,THOR)>0),sum(countOverlaps(THOR,ROTS)>0))/2000,
                                         min(sum(countOverlaps(DB,MAnorm)>0),sum(countOverlaps(MAnorm,DB)>0))/2000,
                                         min(sum(countOverlaps(DB,DR)>0),sum(countOverlaps(DR,DB)>0))/2000,
                                         min(sum(countOverlaps(DB,PePr)>0),sum(countOverlaps(PePr,DB)>0))/2000,
                                         min(sum(countOverlaps(DB,THOR)>0),sum(countOverlaps(THOR,DB)>0))/2000,
                                         min(sum(countOverlaps(MAnorm,DR)>0),sum(countOverlaps(DR,MAnorm)>0))/2000,
                                         min(sum(countOverlaps(MAnorm,PePr)>0),sum(countOverlaps(PePr,MAnorm)>0))/2000,
                                         min(sum(countOverlaps(MAnorm,THOR)>0),sum(countOverlaps(THOR,MAnorm)>0))/2000,
                                         min(sum(countOverlaps(DR,PePr)>0),sum(countOverlaps(PePr,DR)>0))/2000,
                                         min(sum(countOverlaps(DR,THOR)>0),sum(countOverlaps(THOR,DR)>0))/2000,
                                         min(sum(countOverlaps(PePr,THOR)>0),sum(countOverlaps(THOR,PePr)>0))/2000)
    
    my_matrix[upper.tri(my_matrix)] <- c(min(sum(countOverlaps(ROTS,DB)>0),sum(countOverlaps(DB,ROTS)>0))/2000,
                                         min(sum(countOverlaps(ROTS,MAnorm)>0),sum(countOverlaps(MAnorm,ROTS)>0))/2000,
                                         min(sum(countOverlaps(DB,MAnorm)>0),sum(countOverlaps(MAnorm,DB)>0))/2000,
                                         min(sum(countOverlaps(ROTS,DR)>0),sum(countOverlaps(DR,ROTS)>0))/2000,
                                         min(sum(countOverlaps(DB,DR)>0),sum(countOverlaps(DR,DB)>0))/2000,
                                         min(sum(countOverlaps(MAnorm,DR)>0),sum(countOverlaps(DR,MAnorm)>0))/2000,
                                         min(sum(countOverlaps(ROTS,PePr)>0),sum(countOverlaps(PePr,ROTS)>0))/2000,
                                         min(sum(countOverlaps(DB,PePr)>0),sum(countOverlaps(PePr,DB)>0))/2000,
                                         min(sum(countOverlaps(MAnorm,PePr)>0),sum(countOverlaps(PePr,MAnorm)>0))/2000,
                                         min(sum(countOverlaps(DR,PePr)>0),sum(countOverlaps(PePr,DR)>0))/2000,
                                         min(sum(countOverlaps(ROTS,THOR)>0),sum(countOverlaps(THOR,ROTS)>0))/2000,
                                         min(sum(countOverlaps(DB,THOR)>0),sum(countOverlaps(THOR,DB)>0))/2000,
                                         min(sum(countOverlaps(MAnorm,THOR)>0),sum(countOverlaps(THOR,MAnorm)>0))/2000,
                                         min(sum(countOverlaps(DR,THOR)>0),sum(countOverlaps(THOR,DR)>0))/2000,
                                         min(sum(countOverlaps(PePr,THOR)>0),sum(countOverlaps(THOR,PePr)>0))/2000)
  }
  if(switch == 2){
  my_matrix[lower.tri(my_matrix)] <- c(sum(countOverlaps(ROTS,DB)>0)/length(ROTS),
                                       sum(countOverlaps(ROTS,MAnorm)>0)/length(ROTS),
                                       sum(countOverlaps(ROTS,DR)>0)/length(ROTS),
                                       sum(countOverlaps(ROTS,PePr)>0)/length(ROTS),
                                       sum(countOverlaps(ROTS,THOR)>0)/length(ROTS),
                                       sum(countOverlaps(DB,MAnorm)>0)/length(DB),
                                       sum(countOverlaps(DB,DR)>0)/length(DB),
                                       sum(countOverlaps(DB,PePr)>0)/length(DB),
                                       sum(countOverlaps(DB,THOR)>0)/length(DB),
                                       sum(countOverlaps(MAnorm,DR)>0)/length(MAnorm),
                                       sum(countOverlaps(MAnorm,PePr)>0)/length(MAnorm),
                                       sum(countOverlaps(MAnorm,THOR)>0)/length(MAnorm),
                                       sum(countOverlaps(DR,PePr)>0)/length(DR),
                                       sum(countOverlaps(DR,THOR)>0)/length(DR),
                                       sum(countOverlaps(PePr,THOR)>0)/length(PePr))
  
  my_matrix[upper.tri(my_matrix)] <- c(sum(countOverlaps(DB,ROTS)>0)/length(DB),
                                       sum(countOverlaps(MAnorm,ROTS)>0)/length(MAnorm),
                                       sum(countOverlaps(MAnorm,DB)>0)/length(MAnorm),
                                       sum(countOverlaps(DR,ROTS)>0)/length(DR),
                                       sum(countOverlaps(DR,DB)>0)/length(DR),
                                       sum(countOverlaps(DR,MAnorm)>0)/length(DR),
                                       sum(countOverlaps(PePr,ROTS)>0)/length(PePr),
                                       sum(countOverlaps(PePr,DB)>0)/length(PePr),
                                       sum(countOverlaps(PePr,MAnorm)>0)/length(PePr),
                                       sum(countOverlaps(PePr,DR)>0)/length(PePr),
                                       sum(countOverlaps(THOR,ROTS)>0)/length(THOR),
                                       sum(countOverlaps(THOR,DB)>0)/length(THOR),
                                       sum(countOverlaps(THOR,MAnorm)>0)/length(THOR),
                                       sum(countOverlaps(THOR,DR)>0)/length(THOR),
                                       sum(countOverlaps(THOR,PePr)>0)/length(THOR))
  }
  
  colnames(my_matrix) <- c("ROTS","DiffBind","MAnorm","diffReps","PePr","THOR")
  rownames(my_matrix) <- c("ROTS","DiffBind","MAnorm","diffReps","PePr","THOR")
  return(my_matrix)
  

}

#Call the main function
YF_2000 <- main_function(YF_DB,YF_ROTS,YF_MAnorm,YF_DR,YF_PePr,YF_THOR,"YF_",switch =1)
IFN_2000 <- main_function(IFN_DB,IFN_ROTS,IFN_MAnorm,IFN_DR,IFN_PePr,IFN_THOR,"IFN_",switch =1)
K4_2000 <- main_function(K4_DB,K4_ROTS,K4_MAnorm,K4_DR,K4_PePr,K4_THOR,"K4_",switch =1)
K36_2000 <- main_function(K36_DB,K36_ROTS,K36_MAnorm,K36_DR,K36_PePr,K36_THOR,"K36_",switch =1)

YF_05 <- main_function(YF_DB,YF_ROTS,YF_MAnorm,YF_DR,YF_PePr,YF_THOR,"YF_",switch =2)
IFN_05 <- main_function(IFN_DB,IFN_ROTS,IFN_MAnorm,IFN_DR,IFN_PePr,IFN_THOR,"IFN_",switch =2)
K4_05 <- main_function(K4_DB,K4_ROTS,K4_MAnorm,K4_DR,K4_PePr,K4_THOR,"K4_",switch =2)
K36_05 <- main_function(K36_DB,K36_ROTS,K36_MAnorm,K36_DR,K36_PePr,K36_THOR,"K36_",switch =2)


col4 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","#7FFF7F", "cyan", "#007FFF", "blue","#00007F")) 

#############################################################
###                   figure 2                            ###
#############################################################

pdf("correlogramme_2000_IFN_MAnorm_pooled.pdf",width=12,height = 12)
par(cex = 2)
corrplot(IFN_2000,p.mat=IFN_2000,cl.lim=c(0,1), col = col4(500),sig.level = -1,insig = "p-value",tl.col ='black',na.label = 0)
dev.off()

pdf("correlogramme_2000_YF_MAnorm_pooled.pdf",width=12,height = 12)
par(cex = 2)
corrplot(YF_2000,p.mat=YF_2000,cl.lim=c(0,1), col = col4(500),sig.level = -1,insig = "p-value",tl.col ='black',na.label = 0)
dev.off()

pdf("correlogramme_2000_K4_MAnorm_pooled.pdf",width=12,height = 12)
par(cex = 2)
corrplot(K4_2000,p.mat=K4_2000,cl.lim=c(0,1), col = col4(500),sig.level = -1,insig = "p-value",tl.col ='black',na.label = 0)
dev.off()

pdf("correlogramme_2000_K36_MAnorm_pooled.pdf",width=12,height = 12)
par(cex = 2)
corrplot(K36_2000,p.mat=K36_2000,cl.lim=c(0,1), col = col4(500),sig.level = -1,insig = "p-value",tl.col ='black',na.label = 0)
dev.off()

#############################################################
###            Supplementary Figure 2                     ###
#############################################################

pdf("correlogramme_005_IFN_MAnorm_pooled.pdf",width=12,height = 12)
par(cex = 2)
corrplot(IFN_05,p.mat=IFN_05,cl.lim=c(0,1), col = col4(500),sig.level = -1,insig = "p-value",tl.col ='black',na.label = 0)
dev.off()

pdf("correlogramme_005_YF_MAnorm_pooled.pdf",width=12,height = 12)
par(cex = 2)
corrplot(YF_05,p.mat=YF_05,cl.lim=c(0,1), col = col4(500),sig.level = -1,insig = "p-value",tl.col ='black',na.label = 0)
dev.off()

pdf("correlogramme_005_K4_MAnorm_pooled.pdf",width=12,height = 12)
par(cex = 2)
corrplot(K4_05,p.mat=K4_05,cl.lim=c(0,1), col = col4(500),sig.level = -1,insig = "p-value",tl.col ='black',na.label = 0)
dev.off()

pdf("correlogramme_005_K36_MAnorm_pooled.pdf",width=12,height = 12)
par(cex = 2)
corrplot(K36_05,p.mat=K36_2000,cl.lim=c(0,1), col = col4(500),sig.level = -1,insig = "p-value",tl.col ='black',na.label = 0)
dev.off()


