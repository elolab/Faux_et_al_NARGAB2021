png("YF_density_pp.png", width = 1200,height = 1200)

#Load the read counts, the counts are splited in 100 bins
load("RData/YF_ATAC/YF_DR.RData")
YF_DR <- DR
load("RData/YF_ATAC/YF_DB_pp.RData")
YF_DB <- DB
load("RData/YF_ATAC/YF_PePr.RData")
YF_PePr <- PePr
load("RData/YF_ATAC/YF_THOR.RData")
YF_THOR <- THOR
load("RData/YF_ATAC/YF_ROTS_pp.RData")
YF_ROTS <- ROTS
load("RData/YF_ATAC/YF_MAnorm_pp.RData")
YF_MAnorm <- MAnorm

peaks<- YF_ROTS
soft <- "ROTS_pp"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/YF_ATAC/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

plot(x=seq(-2000, 1960, length.out=100),
     y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
     ty='l',
     col='black',lwd=4,
     ylab=' ',
     cex.axis=4,
     xaxt="n",
     xlab=' ')

peaks<- YF_MAnorm
soft <- "MAnorm_pp"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/YF_ATAC/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='violet',lwd=4)

peaks<- YF_DR
soft <- "DR"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/YF_ATAC/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='green',lwd=4)

peaks<- YF_THOR
soft <- "THOR"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/YF_ATAC/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='orange',lwd=4)




peaks<- YF_DB
soft <- "DB_pp"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/YF_ATAC/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='blue',lwd=4)


peaks<- YF_PePr
soft <- "PePr"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/YF_ATAC/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='red',lwd=4)
dev.off()



png("K36_density_pp.png", width = 1200,height = 1200)

load("RData/RA_H3K36me3/K36_DR.RData")
K36_DR <- DR
  load("RData/RA_H3K36me3/K36_DB_pp.RData")
K36_DB <- DB
load("RData/RA_H3K36me3/K36_PePr.RData")
K36_PePr <- PePr
load("RData/RA_H3K36me3/K36_THOR.RData")
K36_THOR <- THOR
load("RData/RA_H3K36me3/K36_ROTS_pp.RData")
K36_ROTS <- ROTS
load("RData/RA_H3K36me3/K36_MAnorm_pp.RData")
K36_MAnorm <- MAnorm

peaks<- K36_ROTS
soft <- "ROTS_pp"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/RA_H3K36me3/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

plot(x=seq(-2000, 1960, length.out=100),
     y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
     ty='l',
     col='black',lwd=4,ylim=c(0,10),
     ylab=' ',
     cex.axis=4,
     xaxt="n",
     xlab=' ')


peaks<- K36_MAnorm
soft <- "MAnorm_pp"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/RA_H3K36me3/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='violet',lwd=4)

peaks<- K36_DR
soft <- "DR"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/RA_H3K36me3/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='green',lwd=4)

peaks<- K36_THOR
soft <- "THOR"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/RA_H3K36me3/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='orange',lwd=4)


peaks<- K36_DB
soft <- "DB_pp"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/RA_H3K36me3/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='blue',lwd=4)


peaks<- K36_PePr
soft <- "PePr"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/RA_H3K36me3/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='red',lwd=4)

dev.off()


png("K4_density_pp.png", width = 1200,height = 1200)

load("RData/RA_H3K4me3/K4_DR.RData")
K4_DR <- DR
load("RData/RA_H3K4me3/K4_DB_pp.RData")
K4_DB <- DB
load("RData/RA_H3K4me3/K4_PePr.RData")
K4_PePr <- PePr
load("RData/RA_H3K4me3/K4_THOR.RData")
K4_THOR <- THOR
load("RData/RA_H3K4me3/K4_ROTS_pp.RData")
K4_ROTS <- ROTS
load("RData/RA_H3K4me3/K4_MAnorm_pp.RData")
K4_MAnorm <- MAnorm

peaks<- K4_THOR
soft <- "THOR"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/RA_H3K4me3/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

plot(x=seq(-2000, 1960, length.out=100),
     y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
     ty='l',
     col='orange',lwd=4,
     ylab=' ',
     cex.axis=4,
     xaxt="n",
     xlab=' ')


peaks<- K4_ROTS
soft <- "ROTS_pp"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/RA_H3K4me3/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='black',lwd=4)

peaks<- K4_MAnorm
soft <- "MAnorm_pp"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/RA_H3K4me3/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='violet',lwd=4)

peaks<- K4_DR
soft <- "DR"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/RA_H3K4me3/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='green',lwd=4)


peaks<- K4_DB
soft <- "DB_pp"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/RA_H3K4me3/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='blue',lwd=4)


peaks<- K4_PePr
soft <- "PePr"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/RA_H3K4me3/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 1955,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='red',lwd=4)

dev.off()


png("IFN_density_pp.png", width = 1200,height = 1200)

load("RData/IFN_ATAC/IFN_DR.RData")
IFN_DR <- DR
load("RData/IFN_ATAC/IFN_DB_pp.RData")
IFN_DB <- DB
load("RData/IFN_ATAC/IFN_PePr.RData")
IFN_PePr <- PePr
load("RData/IFN_ATAC/IFN_THOR.RData")
IFN_THOR <- THOR
load("RData/IFN_ATAC/IFN_ROTS_pp.RData")
IFN_ROTS <- ROTS
load("RData/IFN_ATAC/IFN_MAnorm_pp.RData")
IFN_MAnorm <- MAnorm

peaks<- IFN_THOR
soft <- "THOR"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/IFN_ATAC/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

plot(x=seq(-2000, 1960, length.out=100),
     y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
     ty='l',
     col='orange',lwd=4,
     ylab=' ',
     cex.axis=4,
     xaxt="n",
     xlab=' ')


peaks<- IFN_ROTS
soft <- "ROTS_pp"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/IFN_ATAC/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='black',lwd=4)

peaks<- IFN_MAnorm
soft <- "MAnorm_pp"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/IFN_ATAC/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='violet',lwd=4)
#legend("topleft", legend = c("ROTS", "DiffBind","PePr","diffReps","THOR"),
#       col = c("orange","blue","red","green","black"),
#       lty = 1, lwd = 2, box.col = "white", bg = "white", text.font  = 2)

peaks<- IFN_DR
soft <- "DR"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/IFN_ATAC/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='green',lwd=4)


peaks<- IFN_DB
soft <- "DB_pp"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/IFN_ATAC/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='blue',lwd=4)


peaks<- IFN_PePr
soft <- "PePr"
nbtiles<-100
peaks <- cbind(names(peaks),data.frame(peaks))
peaks$order <- c(1:dim(peaks)[1])
peaksFCup <- peaks[which(peaks$fc >= 0),]
peaksFCdn <- peaks[which(peaks$fc < 0),]
myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])
load(paste0("RData/IFN_ATAC/matrices_",soft,"_4000bp.RData"))
dat <- matrix(data=c(10),nrow = 2000,ncol = 100)
sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

lines(x=seq(-2000, 1960, length.out=100),
      y=colMeans(rbind(sumlist[[1]],sumlist[[2]])),
      ty='l',
      col='red',lwd=4)

dev.off()






