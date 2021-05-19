
plot_correlation <- function(ROTS_anno,DB_anno,MA_anno,DR_anno,PePr_anno,THOR_anno,toptable,dataset,option){
ROTS <- as.data.frame(ROTS_anno)
DB <- as.data.frame(DB_anno)
MAnorm <- as.data.frame(MA_anno)
DR <- as.data.frame(DR_anno)
DR <- DR[order(DR$FDR, decreasing = F),]
DR$rank <- 1:dim(DR)[1]
THOR <- as.data.frame(THOR_anno)
THOR <- THOR[order(THOR$FDR, decreasing = F),]
THOR$rank <- 1:dim(THOR)[1]
PePr <- as.data.frame(PePr_anno)
PePr <- PePr[order(PePr$FDR, decreasing = F),]
PePr$rank <- 1:dim(PePr)[1]

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

MAnorm <- merge(MAnorm,toptable, by.x=option,by.y="ID")
MAnorm <- MAnorm[order(MAnorm$rank, decreasing = F),]
MAnorm_agg <- aggregate(x=MAnorm[,c("fc")],by = list(MAnorm$gene_name,MAnorm$logFC,MAnorm$rank), FUN=mean)
colnames(MAnorm_agg) <- c("gene_name","logFC","rank","fc")

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
cor_MA <- c()
cor_DR <- c()
cor_THOR <- c()
cor_PePr <- c()

#fill the objects with correlation values
for(i in seq(100,2000,100)){
  cor_ROTS <- c(cor_ROTS,mycorrelation_toplist(ROTS_agg,i,"pearson"))
  cor_DB <- c(cor_DB,mycorrelation_toplist(DB_agg,i,"pearson"))
  cor_MA <- c(cor_MA,mycorrelation_toplist(MAnorm_agg,i,"pearson"))
  cor_DR <- c(cor_DR,mycorrelation_toplist(DR_agg,i,"pearson"))
  cor_THOR <- c(cor_THOR,mycorrelation_toplist(THOR_agg,i,"pearson"))
  cor_PePr <- c(cor_PePr,mycorrelation_toplist(PePr_agg,i,"pearson"))
  
}

#plot the correlations curves
pdf(paste0(dataset,"correlation_toplist_narrow_pooledpeaks_aggregated_2000.pdf"),width =12,height =12)
par(mar = c(4, 5, 4, 4))

plot(cor_PePr ,lwd=4, col = "red",type="l", ylim=c(0,1), ylab ="Pearson correlation",mgp = c(3,0.5,0),xlab = " ", main=" ", xaxt="n",cex.axis =3, cex.lab=3)
lines(cor_DB,  lwd=4, col="blue")
lines(cor_MA,  lwd=4, col="violet")
lines(cor_DR,  lwd=4, col="green")
lines(cor_THOR,  lwd=4,col="orange")
lines(cor_ROTS,lwd=4,col="black")
xtick<-seq(0, 20, by=5)
xlabs<-seq(0,2000,500)
axis(side=1, at=xtick,labels = FALSE)
text(x=xtick,  par("usr")[3], 
     labels = xlabs, pos = 1, xpd = TRUE, offset=2 , cex=3)

dev.off()
}





### load to save time, the previous step is quite long
load("RData/IFN_ATAC/ROTS_narrow_pp_anno.R")
IFN_ROTS_anno <- ROTS_anno
load("RData/IFN_ATAC/MAnorm_narrow_pp_anno.R")
IFN_MAnorm_anno <- MAnorm_anno
load("RData/IFN_ATAC/PePr_anno.R")
IFN_PePr_anno <- PePr_anno
load("RData/IFN_ATAC/THOR_anno.R")
IFN_THOR_anno <- THOR_anno
load("RData/IFN_ATAC/DB_narrow_pp_anno.R")
IFN_DB_anno <- DB_anno
load("RData/IFN_ATAC/DR_anno.R")
IFN_DR_anno <- DR_anno
load("RData/IFN_ATAC/toptableGR02GR01.RData")
IFN_toptable <- toptable
### load to save time, the previous step is quite long
load("RData/YF_ATAC/ROTS_narrow_pp_anno.R")
YF_ROTS_anno <- ROTS_anno
load("RData/YF_ATAC/MAnorm_narrow_pp_anno.R")
YF_MAnorm_anno <- MAnorm_anno
load("RData/YF_ATAC/PePr_anno.R")
YF_PePr_anno <- PePr_anno
load("RData/YF_ATAC/THOR_anno.R")
YF_THOR_anno <- THOR_anno
load("RData/YF_ATAC/DB_narrow_pp_anno.R")
YF_DB_anno <- DB_anno
load("RData/YF_ATAC/DR_anno.R")
YF_DR_anno <- DR_anno
load("RData/YF_ATAC/toptableGR02GR03.RData")
YF_toptable <- toptable

load("RData/RA_H3K36me3/ROTS_pp_anno.R")
K36_ROTS_anno <- ROTS_anno
load("RData/RA_H3K36me3/MAnorm_pp_anno.R")
K36_MAnorm_anno <- MAnorm_anno
load("RData/RA_H3K36me3/PePr_anno.R")
K36_PePr_anno <- PePr_anno
load("RData/RA_H3K36me3/THOR_anno.R")
K36_THOR_anno <- THOR_anno
load("RData/RA_H3K36me3/DB_pp_anno.R")
K36_DB_anno <- DB_anno
load("RData/RA_H3K36me3/DR_anno.R")
K36_DR_anno <- DR_anno
load("RData/RA_H3K36me3/toptableGR02GR01.RData")
K36_toptable <- toptable

load("RData/RA_H3K4me3/ROTS_narrow_pp_anno.R")
K4_ROTS_anno <- ROTS_anno
load("RData/RA_H3K4me3/MAnorm_narrow_pp_anno.R")
K4_MAnorm_anno <- MAnorm_anno
load("RData/RA_H3K4me3/PePr_anno.R")
K4_PePr_anno <- PePr_anno
load("RData/RA_H3K4me3/THOR_anno.R")
K4_THOR_anno <- THOR_anno
load("RData/RA_H3K4me3/DB_narrow_pp_anno.R")
K4_DB_anno <- DB_anno
load("RData/RA_H3K4me3/DR_anno.R")
K4_DR_anno <- DR_anno
load("RData/RA_H3K4me3/toptableGR02GR01.RData")
K4_toptable <- toptable

###main
plot_correlation(YF_ROTS_anno,YF_DB_anno,YF_MAnorm_anno,YF_DR_anno,YF_PePr_anno,YF_THOR_anno,YF_toptable,"YF","gene_name")
plot_correlation(IFN_ROTS_anno,IFN_DB_anno,IFN_MAnorm_anno,IFN_DR_anno,IFN_PePr_anno,IFN_THOR_anno,IFN_toptable,"IFN","gene_name")
plot_correlation(K4_ROTS_anno,K4_DB_anno,K4_MAnorm_anno,K4_DR_anno,K4_PePr_anno,K4_THOR_anno,K4_toptable,"K4","Ensembl")
plot_correlation(K36_ROTS_anno,K36_DB_anno,K36_MAnorm_anno,K36_DR_anno,K36_PePr_anno,K36_THOR_anno,K36_toptable,"K36","Ensembl")



pdf("legend.pdf",width = 12, height = 12)
plot.new()
legend("left", legend = c("ROTS", "DiffBind","MAnorm2","diffReps","PePr","THOR"),
     col = c("black","blue","violet","green","red","orange"),
     lty = 1, lwd = 2, box.col = "white", bg = "white", text.font  = 2)
dev.off()