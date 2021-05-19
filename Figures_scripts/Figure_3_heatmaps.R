
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(chipseq)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ChIPpeakAnno)


plotHeatmap <- function(mymatrix,nbtiles,nb,saturation,colors,main,myorder){

  #create some saturation after 50 reads
  mymatrix[mymatrix>saturation]<-saturation
  mymatrix[mymatrix< -saturation] <- -saturation

  ##this is here as a color scale purpose, the actual tiles are value zero and will be zero on the plot
  mymatrix[2000,100] <- -50
  mymatrix[2000,1] <- 50
  ##########################
  image(x=seq(-2000, 1960, length.out=nbtiles),
        y=1:nrow(mymatrix),
        z=t(mymatrix[as.integer(myorder$order),]),
        #z=t(mymatrix[order(rowSums(mymatrix)),]),
        col=colors,
        xlab=' ',
        ylab=' ', lwd=2,cex.lab=3,cex.axis=3,
        main=" ",xaxt="n",yaxt="n")#
  box(col='black', lwd=2)
  at1 <- seq(-2000, 2000, 1000)
  axis(side =1, at1, labels = F)
  abline(h=nb+0.5, lwd=1, col='gray')
}


IFNplotheatmap <- function(soft,peaks,dataset){

  #load the peaks and organize them by width and fold change
  peaks <- cbind(names(peaks),data.frame(peaks))
  nbtiles<-100
  peaks$order <- c(1:dim(peaks)[1])
  peaksFCup <- peaks[which(peaks$fc > 0),]
  peaksFCdn <- peaks[which(peaks$fc < 0),]
  myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])

  #load the read counts
  load(paste0("RData/",dataset,"/matrices_",soft,"_4000bp.RData"))


  png(paste0(dataset,"_",soft,"_4000bp.png"),height =1200, width = 1200)

  #take the average read counts over all the replicates
  dat <- matrix(data=c(8),nrow = 2000,ncol = 100)
  sumlist <- list(Reduce('+',matrices[1:5])/dat,Reduce('+',matrices[6:10])/dat)

  #create the color palette
  colors = colorRampPalette(c('blue','white','red'))(101)

  #Subtract the values between condition 1 and condition 2
  test <- sumlist[[2]]-sumlist[[1]]

  #call the ploting function
  plotHeatmap(test,nbtiles,dim(peaksFCup)[1],50,colors,main=dataset,myorder)
  
  dev.off()
}

RAplotheatmap <- function(soft,peaks,dataset){
  #load the peaks and organize them by width and fold change
  peaks <- cbind(names(peaks),data.frame(peaks))
  peaks <- peaks[1:2000,]
  nbtiles<-100
  peaks$order <- c(1:dim(peaks)[1])
  peaksFCup <- peaks[which(peaks$fc > 0),]
  peaksFCdn <- peaks[which(peaks$fc < 0),]
  myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])

  #load the read counts
  load(paste0("RData/",dataset,"/matrices_",soft,"_4000bp.RData"))
  png(paste0(dataset,"_",soft,"_4000bp.png"),height =1200, width = 1200)

  #take the average read counts over all the replicates
  dat <- matrix(data=c(8),nrow = 2000,ncol = 100)
  sumlist <- list(Reduce('+',matrices[1:10])/dat,Reduce('+',matrices[11:20])/dat)

  #create the color palette
  colors = colorRampPalette(c('blue','white','red'))(101)

  #Subtract the values between condition 1 and condition 2
  test <- sumlist[[2]]-sumlist[[1]]

  #call the ploting function
  plotHeatmap(test,nbtiles,dim(peaksFCup)[1],50,colors,main=dataset,myorder)
  
  dev.off()
}

YFplotheatmap <- function(soft,peaks,dataset){
  #load the peaks and organize them by width and fold change
  peaks <- cbind(names(peaks),data.frame(peaks))
  nbtiles<-100
  peaks$order <- c(1:dim(peaks)[1])
  peaksFCup <- peaks[which(peaks$fc > 0),]
  peaksFCdn <- peaks[which(peaks$fc < 0),]
  myorder <- rbind(peaksFCup[order(peaksFCup$width),],peaksFCdn[order(peaksFCdn$width),])

  #load the read counts
  load(paste0("RData/",dataset,"/matrices_",soft,"_4000bp.RData"))
  png(paste0(dataset,"_",soft,"_4000bp.png"),height =1200, width = 1200)

  #take the average read counts over all the replicates
  dat <- matrix(data=c(8),nrow = 2000,ncol = 100)
  sumlist <- list(Reduce('+',matrices[1:8])/dat,Reduce('+',matrices[9:16])/dat)

  #create the color palette
  colors = colorRampPalette(c('blue','white','red'))(101)

  #Subtract the values between condition 1 and condition 2
  test <- sumlist[[2]]-sumlist[[1]]

  #call the ploting function
  plotHeatmap(test,nbtiles,dim(peaksFCup)[1],50,colors,main=dataset,myorder)
  
  dev.off()
}

###main
load("RData/IFN_ATAC/IFN_DR.RData")
IFN_DR <- DR
load("RData/IFN_ATAC/IFN_PePr.RData")
IFN_PePr <- PePr
load("RData/IFN_ATAC/IFN_THOR.RData")
IFN_THOR <- THOR
load("RData/IFN_ATAC/IFN_ROTS_pp.RData")
IFN_ROTS <- ROTS
load("RData/IFN_ATAC/IFN_MAnorm_pp.RData")
IFN_MAnorm <- MAnorm
load("RData/IFN_ATAC/IFN_DB_pp.RData")
IFN_DB <- DB

IFNplotheatmap("DR",IFN_DR,"IFN_ATAC")
IFNplotheatmap("PePr",IFN_PePr,"IFN_ATAC")
IFNplotheatmap("THOR",IFN_THOR,"IFN_ATAC")
IFNplotheatmap("DB_pp",IFN_DB,"IFN_ATAC")
IFNplotheatmap("ROTS_pp",IFN_ROTS,"IFN_ATAC")
IFNplotheatmap("MAnorm_pp",IFN_MAnorm,"IFN_ATAC")

load("RData/YF_ATAC/YF_DR.RData")
YF_DR <- DR
load("RData/YF_ATAC/YF_PePr.RData")
YF_PePr <- PePr
load("RData/YF_ATAC/YF_THOR.RData")
YF_THOR <- THOR
load("RData/YF_ATAC/YF_ROTS_pp.RData")
YF_ROTS <- ROTS
load("RData/YF_ATAC/YF_MAnorm_pp.RData")
YF_MAnorm <- MAnorm
load("RData/YF_ATAC/YF_DB_pp.RData")
YF_DB <- DB

YFplotheatmap("DR",YF_DR,"YF_ATAC")
YFplotheatmap("PePr",YF_PePr,"YF_ATAC")
YFplotheatmap("THOR",YF_THOR,"YF_ATAC")
YFplotheatmap("DB_pp",YF_DB,"YF_ATAC")
YFplotheatmap("ROTS_pp",YF_ROTS,"YF_ATAC")
YFplotheatmap("MAnorm_pp",YF_MAnorm,"YF_ATAC")

load("RData/RA_H3K4me3/K4_DR.RData")
K4_DR <- DR
load("RData/RA_H3K4me3/K4_PePr.RData")
K4_PePr <- PePr
load("RData/RA_H3K4me3/K4_THOR.RData")
K4_THOR <- THOR
load("RData/RA_H3K4me3/K4_ROTS_pp.RData")
K4_ROTS <- ROTS
load("RData/RA_H3K4me3/K4_DB_pp.RData")
K4_DB <- DB
load("RData/RA_H3K4me3/K4_MAnorm_pp.RData")
K4_MAnorm <- MAnorm

RAplotheatmap("DR",K4_DR,"RA_H3K4me3")
RAplotheatmap("PePr",K4_PePr,"RA_H3K4me3")
RAplotheatmap("THOR",K4_THOR,"RA_H3K4me3")
RAplotheatmap("ROTS_pp",K4_ROTS,"RA_H3K4me3")
RAplotheatmap("DB_pp",K4_DB,"RA_H3K4me3")
RAplotheatmap("MAnorm_pp",K4_MAnorm,"RA_H3K4me3")

load("RData/RA_H3K36me3/K36_DR.RData")
K36_DR <- DR
load("RData/RA_H3K36me3/K36_PePr.RData")
K36_PePr <- PePr
load("RData/RA_H3K36me3/K36_THOR.RData")
K36_THOR <- THOR
load("RData/RA_H3K36me3/K36_ROTS_pp.RData")
K36_ROTS <- ROTS
load("RData/RA_H3K36me3/K36_DB_pp.RData")
K36_DB <- DB
load("RData/RA_H3K36me3/K36_MAnorm_pp.RData")
K36_MAnorm <- MAnorm

RAplotheatmap("DR",K36_DR,"RA_H3K36me3")
RAplotheatmap("PePr",K36_PePr,"RA_H3K36me3")
RAplotheatmap("THOR",K36_THOR,"RA_H3K36me3")
RAplotheatmap("ROTS_pp",K36_ROTS,"RA_H3K36me3")
RAplotheatmap("DB_pp",K36_DB,"RA_H3K36me3")
RAplotheatmap("MAnorm_pp",K36_MAnorm,"RA_H3K36me3")

colors = colorRampPalette(c('blue','white','red'))(101)
image(seq(-50, 50, length.out=100), 1,
      matrix(seq(-50, 50, length.out=100),100,1),
      col = colors,
      xlab='Difference in read number', ylab='',cex.axis=2,cex.lab=2,
      main='', yaxt='n',
      lwd=3, axes=TRUE)
