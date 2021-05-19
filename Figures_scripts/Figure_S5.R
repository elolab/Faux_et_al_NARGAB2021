super_scatter_toplist <- function(TABLE, xlim=NULL, ylim=NULL, size=1, main="No title", method="pearson"){
  
  plot(deduplicate(TABLE)$fc[1:size],
       deduplicate(TABLE)$logFC[1:size],
       xlim=xlim,
       ylim=ylim,
       main=main,
       ylab="Differential RNA-seq logFC",
       xlab="Differential ChIP-seq logFC", cex.lab=2,cex.main=3,cex=2,pch=20)
  
  
  abline(h = 0, col="red", lty=2)
  abline(v = 0, col="red", lty=2)
  
  message(main,", ",size," = ",round(
    cor(
      deduplicate(TABLE)$fc[1:size],
      deduplicate(TABLE)$logFC[1:size],
      method = method
    ),
    digits = 3
  ))
}


###############################################################################
## make an average of the Fold change for the multiple peakse annotated to the
## same gene
###############################################################################
## Merging the data from ChIP-seq and RNA-seq experiment and then ordering
## based on FDR (in ROTS som FDR values can be equal to 0, to deal with the
## ties I generate a rank first and then order them based on that (the
## p-vals can have ties in ROTS also))

ROTS <- merge(ROTS,toptable, by.x=option,by.y="ID")
ROTS <- ROTS[order(ROTS$rank, decreasing = F),]
ROTS_agg <- aggregate(x=ROTS[,c("fc")],by = list(ROTS$gene_name,ROTS$logFC,ROTS$rank), FUN=mean)
colnames(ROTS_agg) <- c("gene_name","logFC","rank","fc")

PePr <- merge(PePr,toptable, by.x=option,by.y="ID")
PePr <- PePr[order(PePr$FDR, decreasing = F),]
PePr_agg <- aggregate(x=PePr[,c("fc")],by = list(PePr$gene_name,PePr$logFC,PePr$rank), FUN=mean)
colnames(PePr_agg) <- c("gene_name","logFC","rank","fc")

DB <- merge(DB,toptable, by.x=option,by.y="ID")
DB <- DB[order(DB$FDR, decreasing = F),]
DB_agg <- aggregate(x=DB[,c("fc")],by = list(DB$gene_name,DB$logFC,DB$rank), FUN=mean)
colnames(DB_agg) <- c("gene_name","logFC","rank","fc")

DR <- merge(DR,toptable, by.x=option,by.y="ID")
DR <- DR[order(DR$FDR, decreasing = F),]
DR_agg <- aggregate(x=DR[,c("fc")],by = list(DR$gene_name,DR$logFC,DR$rank), FUN=mean)
colnames(DR_agg) <- c("gene_name","logFC","rank","fc")

THOR <- merge(THOR,toptable, by.x=option,by.y="ID")
THOR <- THOR[order(THOR$FDR, decreasing = F),]
THOR_agg <- aggregate(x=THOR[,c("fc")],by = list(THOR$gene_name,THOR$logFC,THOR$rank), FUN=mean)
colnames(THOR_agg) <- c("gene_name","logFC","rank","fc")

#remove the peaks that are annotated to the same genes (takes the first occurence)
deduplicate <- function(test){
  test <- test[!duplicated(test$gene_name),]
  return(test)
}


### load to save time, the previous step is quite long
load("RData/IFN_ATAC/ROTS_narrow_pp_anno.R")
IFN_ROTS_anno <- ROTS_anno
load("RData/IFN_ATAC/PePr_anno.R")
IFN_PePr_anno <- PePr_anno
load("RData/IFN_ATAC/THOR_anno.R")
IFN_THOR_anno <- THOR_anno
load("RData/IFN_ATAC/DB_narrow_pp_anno.R")
IFN_DB_anno <- DB_anno
load("RData/IFN_ATAC/MAnorm_narrow_pp_anno.R")
IFN_MAnorm_anno <- MAnorm_anno
load("RData/IFN_ATAC/DR_anno.R")
IFN_DR_anno <- DR_anno
load("RData/IFN_ATAC/toptableGR02GR01.RData")
IFN_toptable <- toptable
### load to save time, the previous step is quite long
load("RData/YF_ATAC/ROTS_narrow_pp_anno.R")
YF_ROTS_anno <- ROTS_anno
load("RData/YF_ATAC/PePr_anno.R")
YF_PePr_anno <- PePr_anno
load("RData/YF_ATAC/THOR_anno.R")
YF_THOR_anno <- THOR_anno
load("RData/YF_ATAC/DB_narrow_pp_anno.R")
YF_DB_anno <- DB_anno
load("RData/YF_ATAC/MAnorm_narrow_pp_anno.R")
YF_MAnorm_anno <- MAnorm_anno
load("RData/YF_ATAC/DR_anno.R")
YF_DR_anno <- DR_anno
load("RData/YF_ATAC/toptableGR02GR03.RData")
YF_toptable <- toptable

load("RData/RA_H3K36me3/ROTS_pp_anno.R")
K36_ROTS_anno <- ROTS_anno
load("RData/RA_H3K36me3/PePr_anno.R")
K36_PePr_anno <- PePr_anno
load("RData/RA_H3K36me3/THOR_anno.R")
K36_THOR_anno <- THOR_anno
load("RData/RA_H3K36me3/DB_pp_anno.R")
K36_DB_anno <- DB_anno
load("RData/RA_H3K36me3/MAnorm_pp_anno.R")
K36_MAnorm_anno <- MAnorm_anno
load("RData/RA_H3K36me3/DR_anno.R")
K36_DR_anno <- DR_anno
load("RData/RA_H3K36me3/toptableGR02GR01.RData")
K36_toptable <- toptable

load("RData/RA_H3K4me3/ROTS_narrow_pp_anno.R")
K4_ROTS_anno <- ROTS_anno
load("RData/RA_H3K4me3/PePr_anno.R")
K4_PePr_anno <- PePr_anno
load("RData/RA_H3K4me3/THOR_anno.R")
K4_THOR_anno <- THOR_anno
load("RData/RA_H3K4me3/DB_narrow_pp_anno.R")
K4_DB_anno <- DB_anno
load("RData/RA_H3K4me3/MAnorm_narrow_pp_anno.R")
K4_MAnorm_anno <- MAnorm_anno
load("RData/RA_H3K4me3/DR_anno.R")
K4_DR_anno <- DR_anno
load("RData/RA_H3K4me3/toptableGR02GR01.RData")
K4_toptable <- toptable


#function preprocessing the data for ATAC-seq and calling the plotting function
call_super_scatter <- function(ROTS, DB,MAnorm, DR, THOR, PePr, toptable, xlim, ylim, dataset){
  ROTS <- as.data.frame(ROTS)
  DB <- as.data.frame(DB)
  MAnorm <- as.data.frame(MAnorm)
  DR <- as.data.frame(DR)
  THOR <- as.data.frame(THOR)
  THOR <- THOR[order(THOR$FDR, decreasing = F),]
  THOR$rank <- 1:dim(THOR)[1]
  PePr <- as.data.frame(PePr)
  PePr <- PePr[order(PePr$FDR, decreasing = F),]
  PePr$rank <- 1:dim(PePr)[1]
  
  ROTS <- merge(ROTS,toptable, by.x="gene_name",by.y="ID")
  ROTS <- ROTS[order(ROTS$rank, decreasing = F),]
  ROTS_agg <- aggregate(x=ROTS[,c("fc")],by = list(ROTS$gene_name,ROTS$logFC,ROTS$rank), FUN=mean)
  colnames(ROTS_agg) <- c("gene_name","logFC","rank","fc")
  
  PePr <- merge(PePr,toptable, by.x="gene_name",by.y="ID")
  PePr <- PePr[order(PePr$FDR, decreasing = F),]
  PePr_agg <- aggregate(x=PePr[,c("fc")],by = list(PePr$gene_name,PePr$logFC,PePr$rank), FUN=mean)
  colnames(PePr_agg) <- c("gene_name","logFC","rank","fc")
  
  DB <- merge(DB,toptable, by.x="gene_name",by.y="ID")
  DB <- DB[order(DB$FDR, decreasing = F),]
  DB_agg <- aggregate(x=DB[,c("fc")],by = list(DB$gene_name,DB$logFC,DB$rank), FUN=mean)
  colnames(DB_agg) <- c("gene_name","logFC","rank","fc")
  
  MAnorm <- merge(MAnorm,toptable, by.x="gene_name",by.y="ID")
  MAnorm <- MAnorm[order(MAnorm$FDR, decreasing = F),]
  MAnorm_agg <- aggregate(x=MAnorm[,c("fc")],by = list(MAnorm$gene_name,MAnorm$logFC,MAnorm$rank), FUN=mean)
  colnames(MAnorm_agg) <- c("gene_name","logFC","rank","fc")
  
  DR <- merge(DR,toptable, by.x="gene_name",by.y="ID")
  DR <- DR[order(DR$FDR, decreasing = F),]
  DR_agg <- aggregate(x=DR[,c("fc")],by = list(DR$gene_name,DR$logFC,DR$rank), FUN=mean)
  colnames(DR_agg) <- c("gene_name","logFC","rank","fc")
  
  THOR <- merge(THOR,toptable, by.x="gene_name",by.y="ID")
  THOR <- THOR[order(THOR$FDR, decreasing = F),]
  THOR_agg <- aggregate(x=THOR[,c("fc")],by = list(THOR$gene_name,THOR$logFC,THOR$rank), FUN=mean)
  colnames(THOR_agg) <- c("gene_name","logFC","rank","fc")
  

  for( i in c(100,2000)){
    png(file=paste0("scatterplot_",i,"_toplist_",dataset,".png"),width =1000,height = 1000)
    layout(matrix(c(1,2,3,4,5,6), 2, 3), widths=c(1,1))
    par(mar = c(5, 5, 5, 5))
    
    super_scatter_toplist(TABLE = ROTS_agg,
                          xlim = xlim,
                          ylim = ylim,
                          size = i,
                          main = paste0("ROTS ",i),
                          method = "pearson"
    )
    super_scatter_toplist(TABLE = DB_agg,
                          xlim = xlim,
                          ylim = ylim,
                          size = i,
                          main = paste0("DiffBind ",i),
                          method = "pearson"
    )
    super_scatter_toplist(TABLE = MAnorm_agg,
                          xlim = xlim,
                          ylim = ylim,
                          size = i,
                          main = paste0("MAnorm2 ",i),
                          method = "pearson"
    )
    super_scatter_toplist(TABLE =  DR_agg,
                          xlim = xlim,
                          ylim = ylim,
                          size = i,
                          main = paste0("diffReps ",i),
                          method = "pearson"
    )
    super_scatter_toplist(TABLE =  THOR_agg,
                          xlim = xlim,
                          ylim = ylim,
                          size = i,
                          main = paste0("THOR ",i),
                          method = "pearson"
    )
    super_scatter_toplist(TABLE =  PePr_agg,
                          xlim = xlim,
                          ylim = ylim,
                          size = i,
                          main = paste0("PePr ",i),
                          method = "pearson"
    )
    dev.off()
  }
}

#function preprocessing the data for ChIP-seq and calling the plotting function
call_super_scatterRA <- function(ROTS, DB,MAnorm, DR, THOR, PePr, toptable, xlim, ylim, dataset){
  ROTS <- as.data.frame(ROTS)
  DB <- as.data.frame(DB)
  MAnorm <- as.data.frame(MAnorm)
  DR <- as.data.frame(DR)
  THOR <- as.data.frame(THOR)
  THOR <- THOR[order(THOR$FDR, decreasing = F),]
  THOR$rank <- 1:dim(THOR)[1]
  PePr <- as.data.frame(PePr)
  PePr <- PePr[order(PePr$FDR, decreasing = F),]
  PePr$rank <- 1:dim(PePr)[1]
  
  ROTS <- merge(ROTS,toptable, by.x="Ensembl",by.y="ID")
  ROTS <- ROTS[order(ROTS$rank, decreasing = F),]
  ROTS_agg <- aggregate(x=ROTS[,c("fc")],by = list(ROTS$gene_name,ROTS$logFC,ROTS$rank), FUN=mean)
  colnames(ROTS_agg) <- c("gene_name","logFC","rank","fc")
  
  PePr <- merge(PePr,toptable, by.x="Ensembl",by.y="ID")
  PePr <- PePr[order(PePr$FDR, decreasing = F),]
  PePr_agg <- aggregate(x=PePr[,c("fc")],by = list(PePr$gene_name,PePr$logFC,PePr$rank), FUN=mean)
  colnames(PePr_agg) <- c("gene_name","logFC","rank","fc")
  
  DB <- merge(DB,toptable, by.x="Ensembl",by.y="ID")
  DB <- DB[order(DB$FDR, decreasing = F),]
  DB_agg <- aggregate(x=DB[,c("fc")],by = list(DB$gene_name,DB$logFC,DB$rank), FUN=mean)
  colnames(DB_agg) <- c("gene_name","logFC","rank","fc")
  
  MAnorm <- merge(MAnorm,toptable, by.x="Ensembl",by.y="ID")
  MAnorm <- MAnorm[order(MAnorm$FDR, decreasing = F),]
  MAnorm_agg <- aggregate(x=MAnorm[,c("fc")],by = list(MAnorm$gene_name,MAnorm$logFC,MAnorm$rank), FUN=mean)
  colnames(MAnorm_agg) <- c("gene_name","logFC","rank","fc")
  
  DR <- merge(DR,toptable, by.x="Ensembl",by.y="ID")
  DR <- DR[order(DR$FDR, decreasing = F),]
  DR_agg <- aggregate(x=DR[,c("fc")],by = list(DR$gene_name,DR$logFC,DR$rank), FUN=mean)
  colnames(DR_agg) <- c("gene_name","logFC","rank","fc")
  
  THOR <- merge(THOR,toptable, by.x="Ensembl",by.y="ID")
  THOR <- THOR[order(THOR$FDR, decreasing = F),]
  THOR_agg <- aggregate(x=THOR[,c("fc")],by = list(THOR$gene_name,THOR$logFC,THOR$rank), FUN=mean)
  colnames(THOR_agg) <- c("gene_name","logFC","rank","fc")
  

  for( i in c(100,2000)){
    png(file=paste0("scatterplot_",i,"_toplist_",dataset,".png"),width =1000,height = 1000)
    layout(matrix(c(1,2,3,4,5,6), 2, 3), widths=c(1,1))
    par(mar = c(5, 5, 5, 5))
    
    super_scatter_toplist(TABLE = ROTS_agg,
                          xlim = xlim,
                          ylim = ylim,
                          size = i,
                          main = paste0("ROTS ",i),
                          method = "pearson"
    )
    super_scatter_toplist(TABLE = DB_agg,
                          xlim = xlim,
                          ylim = ylim,
                          size = i,
                          main = paste0("DiffBind ",i),
                          method = "pearson"
    )
    super_scatter_toplist(TABLE = MAnorm_agg,
                          xlim = xlim,
                          ylim = ylim,
                          size = i,
                          main = paste0("MAnorm2 ",i),
                          method = "pearson"
    )
    super_scatter_toplist(TABLE =  DR_agg,
                          xlim = xlim,
                          ylim = ylim,
                          size = i,
                          main = paste0("diffReps ",i),
                          method = "pearson"
    )
    super_scatter_toplist(TABLE =  THOR_agg,
                          xlim = xlim,
                          ylim = ylim,
                          size = i,
                          main = paste0("THOR ",i),
                          method = "pearson"
    )
    super_scatter_toplist(TABLE =  PePr_agg,
                          xlim = xlim,
                          ylim = ylim,
                          size = i,
                          main = paste0("PePr ",i),
                          method = "pearson"
    )
    dev.off()
  }
}




######main
call_super_scatter(ROTS = YF_ROTS_anno,
                   DB = YF_DB_anno,
                   MAnorm = YF_MAnorm_anno,
                   DR = YF_DR_anno,
                   THOR = YF_THOR_anno,
                   PePr = YF_PePr_anno,
                   toptable = YF_toptable,
                   xlim = c(-6,6), 
                   ylim = c(),
                   dataset = "YF",
                   inv = TRUE)

call_super_scatter(ROTS = IFN_ROTS_anno,
                   DB = IFN_DB_anno,
                   MAnorm = IFN_MAnorm_anno,
                   DR = IFN_DR_anno,
                   THOR = IFN_THOR_anno,
                   PePr = IFN_PePr_anno,
                   toptable = IFN_toptable,
                   xlim = c(-6,6), 
                   ylim = c(),
                   dataset = "IFN",
                   inv = FALSE)

call_super_scatterRA(ROTS = K4_ROTS_anno,
                   DB = K4_DB_anno,
                   MAnorm = K4_MAnorm_anno,
                   DR = K4_DR_anno,
                   THOR = K4_THOR_anno,
                   PePr = K4_PePr_anno,
                   toptable = K4_toptable,
                   xlim = c(-5,5), 
                   ylim = c(-5,5),
                   dataset = "K4",
                   inv = TRUE)

call_super_scatterRA(ROTS = K36_ROTS_anno,
                   DB = K36_DB_anno,
                   MAnorm = K36_MAnorm_anno,
                   DR = K36_DR_anno,
                   THOR = K36_THOR_anno,
                   PePr = K36_PePr_anno,
                   toptable = K36_toptable,
                   xlim = c(-3,3), 
                   ylim = c(-3,3),
                   dataset = "K36",
                   inv = TRUE)