deduplicate <- function(test){
  test <- test[!duplicated(test$gene_name),]
  return(test)
}

super_scatter_toplist <- function(TABLE, xlim=NULL, ylim=NULL, size=1, main="No title", method="pearson"){
  
  plot(deduplicate(TABLE)$fc[1:size],
       deduplicate(TABLE)$logFC[1:size],
       xlim=xlim,
       ylim=ylim,
       main=main,
       ylab="Differential RNA-seq logFC",
       xlab="Differential ChIP-seq logFC", cex.lab=2,cex.main=3,cex=2,pch=20)
  
  
  legend("bottomright",
         legend = c(
           paste0("global corr : ",
                  round(
                    cor(
                      deduplicate(TABLE)$fc[1:size],
                      deduplicate(TABLE)$logFC[1:size],
                      method = method
                    ),
                    digits = 3
                  ))
         ), 
         
         box.col = "transparent", 
         bg = "transparent", 
         text.font  = 2)
  legend("topright",
         legend = c(
           paste0("UR : ",
                  round(length(which(deduplicate(TABLE)[1:size,"logFC"] > 0 & deduplicate(TABLE)[1:size,"fc"] > 0))/size,
                        digits = 2)
           ),
           paste0("UL : ",
                  round(length(which(deduplicate(TABLE)[1:size,"logFC"] > 0 & deduplicate(TABLE)[1:size,"fc"] < 0))/size,
                        digits = 2)
           ),
           paste0("DR : ",
                  round(length(which(deduplicate(TABLE)[1:size,"logFC"] < 0 & deduplicate(TABLE)[1:size,"fc"] > 0))/size,
                        digits = 2)
           ),
           paste0("DL : ",
                  round(length(which(deduplicate(TABLE)[1:size,"logFC"] < 0 & deduplicate(TABLE)[1:size,"fc"] < 0))/size,
                        digits = 2)
           )
         ), 
         box.col = "transparent", 
         bg = "transparent", 
         text.font  = 2)
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


###################################################################################
#correlation x first of the top list with agregated fold change
mycorrelation_toplist <- function(table,size,method){
  mycor<-round(
    cor(
      deduplicate(table)$fc[1:size],
      deduplicate(table)$logFC[1:size],
      method = method
    ),
    digits = 3
  )
  return(mycor)
}

#remove the peaks that are annotated to the same genes (takes the first occurence)
deduplicate <- function(test){
  test <- test[!duplicated(test$gene_name),]
  return(test)
}

cor_ROTS <- c()
cor_DB <- c()
cor_DR <- c()
cor_THOR <- c()
cor_PePr <- c()
for(i in seq(100,2000,100)){
  cor_ROTS <- c(cor_ROTS,mycorrelation_toplist(ROTS_agg,i,"pearson"))
  cor_DB <- c(cor_DB,mycorrelation_toplist(DB_agg,i,"pearson"))
  cor_DR <- c(cor_DR,mycorrelation_toplist(DR_agg,i,"pearson"))
  cor_THOR <- c(cor_THOR,mycorrelation_toplist(THOR_agg,i,"pearson"))
  cor_PePr <- c(cor_PePr,mycorrelation_toplist(PePr_agg,i,"pearson"))
  
}

pdf(paste0(dataset,"correlation_toplist_summary_aggregated_2000.pdf"),width =12,height =12)
par(mar = c(4, 5, 4, 4))
#layout(matrix(c(1,2), 1, 2), widths=c(3,1))

plot(cor_PePr ,lwd=4, col = "red",type="l", ylim=c(0,1), ylab ="Pearson correlation",xlab = " ", main=" ", xaxt="n",cex.axis =3, cex.lab=3)
lines(cor_DB,  lwd=4, col="blue")
lines(cor_DR,  lwd=4, col="green")
lines(cor_THOR,  lwd=4,col="black")
lines(cor_ROTS,lwd=4,col="orange")
xtick<-seq(0, 20, by=5)
xlabs<-seq(0,2000,500)
axis(side=1, at=xtick,labels = FALSE)
text(x=xtick,  par("usr")[3], 
     labels = xlabs, pos = 1, xpd = TRUE, offset=2, cex=3)




#plot.new()
#legend("left", legend = c("ROTS", "DiffBind","PePr","diffReps","THOR"),
#      col = c("orange","blue","red","green","black"),
#      lty = 1, lwd = 2, box.col = "white", bg = "white", text.font  = 2)
#dev.off()
#}





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



call_super_scatter <- function(ROTS, DB,MAnorm, DR, THOR, PePr, toptable, xlim, ylim, dataset, inv){
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
  
  if (inv==TRUE){
    toptable$logFC <- toptable$logFC*(-1)
  }
  
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
  
  #YF ATAC xlim = c(-12,18)
  #YF ATAC ylim = c(-7,11)
  #IFN ATAC xlim = c(-12,18)
  #IFN ATAC ylim = c(-7,11)
  #RA K4 xlim = c(-5,5)
  #RA K4 ylim = c(-5,5)
  #RA K36 xlim = c(-3,3)
  #RA K36 ylim = c(-3,3)
  #c(100,200,300,400,500,1000,1500,2000)
  for( i in c(100,1000)){
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
call_super_scatterRA <- function(ROTS, DB,MAnorm, DR, THOR, PePr, toptable, xlim, ylim, dataset, inv){
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
  
  if (inv==TRUE){
    toptable$logFC <- toptable$logFC*(-1)
  }
  
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
  
  #YF ATAC xlim = c(-12,18)
  #YF ATAC ylim = c(-7,11)
  #IFN ATAC xlim = c(-12,18)
  #IFN ATAC ylim = c(-7,11)
  #RA K4 xlim = c(-5,5)
  #RA K4 ylim = c(-5,5)
  #RA K36 xlim = c(-3,3)
  #RA K36 ylim = c(-3,3)
  for( i in c(100,200,300,400,500,1000,1500,2000)){
    png(file=paste0("scatterplot_",i,"_toplist_",dataset,".png"),width =1000,height = 1000)
    layout(matrix(c(1,2,3,4,5,6), 2, 3), widths=c(1,1))
    par(mar = c(4, 4, 4, 4))
    
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