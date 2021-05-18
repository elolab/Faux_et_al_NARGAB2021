# # boxplot of the width of the different peaksfor all the datasets
# 
# 
# # loading YF data
# 
 ROTS <- read.table("YF_ATAC/ROTS_YF_narrow_pooledpeaks.bed",sep=" ",header=T)
 DB <- read.table("YF_ATAC/DBsites_YF_deseq2_narrow_pooledpeaks.bed",sep="\t",header=F)
 MAnorm <- read.table("YF_ATAC/MAnorm2_YF_pooledpeaks.bed",sep="\t",header=F)
 DR  <- read.table("YF_ATAC/diffReps.bed",sep="\t",header=F)
 PePr <- read.table("YF_ATAC/PePr2.bed",sep="\t",header=F)
 THOR <- read.table("YF_ATAC/THOR.bed")

# Sort the output when needed
 THOR <- THOR[order(THOR$V7, decreasing = F),]
 THOR$rank <- 1:dim(THOR)[1]
 PePr <- PePr[order(PePr$V8, decreasing = F),]
 PePr$rank <- 1:dim(PePr)[1]

# Take the 2000 first peaks for each software
 ROTS <- ROTS[1:2000,c(2:4,12)]
 DB <- DB[1:2000,c(1:3,6)]
 MAnorm <- MAnorm[1:2000,c(1:3,9)]
 DR <- DR[1:2000,c(1:3,8)]
 PePr <- PePr[1:2000,c(1:3,9)]
 THOR <- THOR[1:2000,c(1:3,8)]

# Calculate the width of the peaks
 ROTS[,5] <- ROTS[,3]-ROTS[,2]
 DB[,5] <- DB[,3]-DB[,2]
 MAnorm[,5] <- MAnorm[,3]-MAnorm[,2]
 DR[,5] <- DR[,3]-DR[,2]
 THOR[,5] <- THOR[,3]-THOR[,2]
 PePr[,5] <- PePr[,3]-PePr[,2]

# Create an object
 ROTS_YF <- ROTS
 DB_YF <- DB
 MAnorm_YF <- MAnorm
 DR_YF <- DR
 THOR_YF <- THOR
 PePr_YF <- PePr
# 
# # loading IFN data
# 
 ROTS <- read.table("IFN_ATAC/ROTS_IFN_narrow_pooledpeaks.bed",sep=" ",header=T)
 DB <- read.table("IFN_ATAC/DBsites_IFN_deseq2_narrow_pooledpeaks.bed",sep="\t",header=F)
 MAnorm <- read.table("IFN_ATAC/MAnorm2_IFN_pooledpeaks.bed",sep="\t",header=F)
 DR  <- read.table("IFN_ATAC-seq/diffReps.bed",sep="\t",header=F)
 PePr <- read.table("IFN_ATAC-seq/PePr.bed",sep="\t",header=F)
 THOR <- read.table("IFN_ATAC-seq/THOR.bed")
# 
# Sort the output when needed
 THOR <- THOR[order(THOR$V7, decreasing = F),]
 THOR$rank <- 1:dim(THOR)[1]
 PePr <- PePr[order(PePr$V8, decreasing = F),]
 PePr$rank <- 1:dim(PePr)[1]

# Take the 2000 first peaks for each software
 ROTS <- ROTS[1:2000,c(2:4,12)]
 DB <- DB[1:2000,c(1:3,6)]
 MAnorm <- MAnorm[1:2000,c(1:3,9)]
 DR <- DR[1:2000,c(1:3,8)]
 PePr <- PePr[1:2000,c(1:3,9)]
 THOR <- THOR[1:2000,c(1:3,8)]

# Calculate the width of the peaks
 ROTS[,5] <- ROTS[,3]-ROTS[,2]
 DB[,5] <- DB[,3]-DB[,2]
 MAnorm[,5] <- MAnorm[,3]-MAnorm[,2]
 DR[,5] <- DR[,3]-DR[,2]
 THOR[,5] <- THOR[,3]-THOR[,2]
 PePr[,5] <- PePr[,3]-PePr[,2]

# Create an object
 ROTS_IFN <- ROTS
 DB_IFN <- DB
 MAnorm_IFN <- MAnorm
 DR_IFN <- DR
 THOR_IFN <- THOR
 PePr_IFN <- PePr

# loading H3K4me3 data
 ROTS <- read.table("RA_H3K4me3/ROTS_K4_narrow_pooledpeaks.bed",sep=" ",header=T)
 DB <- read.table("RA_H3K4me3/DBsites_H3K4me3_deseq2_narrow_pooledpeaks.bed",sep="\t",header=F)
 MAnorm <- read.table("RA_H3K4me3/MAnorm2_K4_pooledpeaks.bed",sep="\t",header=F)
 DR  <- read.table("RA_H3K4me3/diffReps.bed",sep="\t",header=F)
 PePr <- read.table("RA_H3K4me3/PePr.bed",sep="\t",header=F)
 THOR <- read.table("RA_H3K4me3/THOR.bed")
 
 # Sort the output when needed
 THOR <- THOR[order(THOR$V7, decreasing = F),]
 THOR$rank <- 1:dim(THOR)[1]
 PePr <- PePr[order(PePr$V8, decreasing = F),]
 PePr$rank <- 1:dim(PePr)[1]


# Take the 2000 first peaks for each software
 ROTS <- ROTS[1:2000,c(2:4,12)]
 DB <- DB[1:2000,c(1:3,6)]
 DMAnorm <- MAnorm[1:2000,c(1:3,9)]
 DR <- DR[1:2000,c(1:3,8)]
 PePr <- PePr[1:2000,c(1:3,9)]
 THOR <- THOR[1:2000,c(1:3,8)]

# Calculate the width of the peaks
 ROTS[,5] <- ROTS[,3]-ROTS[,2]
 DB[,5] <- DB[,3]-DB[,2]
 MAnorm[,5] <- MAnorm[,3]-MAnorm[,2]
 DR[,5] <- DR[,3]-DR[,2]
 THOR[,5] <- THOR[,3]-THOR[,2]

# Create an object
 ROTS_K4 <- ROTS
 DB_K4 <- DB
 MAnorm_K4 <- MAnorm
 DR_K4 <- DR
 THOR_K4 <- THOR
 PePr_K4 <- PePr

# loading H3K36me3 data
 ROTS <- read.table("RA_H3K36me3/ROTS_K36_pooledpeaks.bed",sep=" ",header=T)
 DB <- read.table("RA_H3K36me3/DBsitesH3K36me3_deseq2_pooledpeaks.bed",sep="\t",header=F)
 MAnorm <- read.table("RA_H3K36me3/MAnorm2_K36_pooledpeaks.bed",sep="\t",header=F)
 DR  <- read.table("RA_H3K36me3/diffReps.bed",sep="\t",header=F)
 PePr <- read.table("RA_H3K36me3/PePr.bed",sep="\t",header=F)
 THOR <- read.table("RA_H3K36me3/THOR.bed")
 
 # Sort the output when needed
 THOR <- THOR[order(THOR$V7, decreasing = F),]
 THOR$rank <- 1:dim(THOR)[1]
 PePr <- PePr[order(PePr$V8, decreasing = F),]
 PePr$rank <- 1:dim(PePr)[1]
 
 # Take the 2000 first peaks for each software
 ROTS <- ROTS[1:2000,c(2:4,12)]
 DB <- DB[1:2000,c(1:3,6)]
 DB <- DB[1:2000,c(1:3,9)]
 DR <- DR[1:2000,c(1:3,8)]
 PePr <- PePr[1:2000,c(1:3,9)]
 THOR <- THOR[1:2000,c(1:3,8)]
 
 # Calculate the width of the peaks
 ROTS[,5] <- ROTS[,3]-ROTS[,2]
 DB[,5] <- DB[,3]-DB[,2]
 MAnorm[,5] <- MAnorm[,3]-MAnorm[,2]
 DR[,5] <- DR[,3]-DR[,2]
 THOR[,5] <- THOR[,3]-THOR[,2]
 PePr[,5] <- PePr[,3]-PePr[,2]
 
 # Create an object
 ROTS_K36 <- ROTS
 DB_K36 <- DB
 MAnorm_K36 <- MAnorm
 DR_K36 <- DR
 THOR_K36 <- THOR
 PePr_K36 <- PePr
 

 #boxplot call

par(mar=c(12,11,3,3))
boxplot(ROTS_YF[,5], DB_YF[,5],MAnorm_YF[,5], DR_YF$V5, PePr_YF$V5, THOR_YF$V5,
        ROTS_IFN[,5],DB_IFN[,5],MAnorm_IFN[,5], DR_IFN$V5, PePr_IFN$V5, THOR_IFN$V5,
        ROTS_K4[,5], DB_K4[,5],MAnorm_K4[,5], DR_K4$V5, PePr_K4$V5, THOR_K4$V5,
        ROTS_K36[,5], DB_K36[,5],MAnorm_K36[,5],DR_K36$V5, PePr_K36$V5, THOR_K36$V5,
        main = " ",
        at = c(1,2,3,4,5,6,8,9,10,11,12,13,15,16,17,18,19,20,22,23,24,25,26,27),
        names = c(rep(c("ROTS", "DiffBind","MAnorm2", "diffReps","PePr","THOR"),4)),
        las = 2,
        col = c(rep(c("grey","lightblue","violetred","lightgreen","tomato","orange"),5)),
        border = c(rep(c("black","blue","violet","green","red","Chocolate1"),5)),
        notch = TRUE,
        ylim = c(0,10000),
        ylab = "",
        cex.lab=3,
        cex.axis=3
)



