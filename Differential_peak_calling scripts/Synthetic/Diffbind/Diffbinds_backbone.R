library(DiffBind)
library(DESeq2)
library(edgeR)

sheet <- dba(sampleSheet="diffbinds_H3K36me3.csv")
print("data loading done")

olap.rate <- dba.overlap(sheet,mode=DBA_OLAP_RATE)
#jpeg("overlapRateH3K36me3.jpeg")
#plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets')
#dev.off()

peakset <- dba.peakset(sheet, consensus=DBA_CONDITION)
#jpeg("peakOverlapH3K36me3.jpeg")
#dba.plotVenn(peakset,peakset$masks$Consensus)
#dev.off()

H3K36me3.count <- dba.count(sheet)
save(H3K36me3.count,file="H3K36me3.count.RData")
print("count done")

H3K36me3.contrast <- dba.contrast(H3K36me3.count, categories=DBA_CONDITION)
save(H3K36me3.contrast,file="H3K36me3.contrast.RData")
print("contrast done")

H3K36me3.analysed <- dba.analyze(H3K36me3.contrast,bParallel=FALSE,method=c(DBA_EDGER, DBA_DESEQ2))
save(H3K36me3.analysed,file="H3K36me3.analysed.RData")
print("analyze done")

#report <- dba.report(H3K36me3.analysed, th=1, DataType=DBA_DATA_FRAME,method=DBA_EDGER)
#score <- -10*(log10(report$FDR))
#write.table(cbind(report[,1:3],rownames(report),score),
#              "test_DBsitesH3K36me3_edger.bed", quote=FALSE, sep="\t",
#             row.names=FALSE, col.names=FALSE)

report <- dba.report(H3K36me3.analysed, th=1, DataType=DBA_DATA_FRAME,method=DBA_DESEQ2)
FDR <- report$FDR
pval <- report$`p-value`
#row.names(report) <- c("chr", "start","end","size","pvalue","FDR")
write.table(cbind(report[,1:3],rownames(report),pval,FDR),
              "DBsitesH3K36me3_deseq2.bed", quote=FALSE, sep="\t",
             row.names=FALSE, col.names=TRUE)

#report <- dba.report(H3K36me3.analysed, th=.05, DataType=DBA_DATA_FRAME,method=DBA_EDGER,bNormalized = F)
#score <- -10*(log10(report$FDR))
#write.table(cbind(report[,1:3],rownames(report),score),
#            "DBsitesH3K36me3_edger_no_normalization.bed", quote=FALSE, sep="\t",
#            row.names=FALSE, col.names=FALSE)

#report <- dba.report(H3K36me3.analysed, th=.05, DataType=DBA_DATA_FRAME,method=DBA_DESEQ2,bNormalized = F)
#score <- -10*(log10(report$FDR))
#write.table(cbind(report[,1:3],rownames(report),score),
#            "DBsitesH3K36me3_deseq2_no_normalization.bed", quote=FALSE, sep="\t",
#            row.names=FALSE, col.names=FALSE)

print("report done")


# par(mfrow=c(2,2))
#     dba.plotMA(H3K36me3.analysed, contrast=1, method=DBA_EDGER, bNormalized=FALSE)
#     dba.plotMA(H3K36me3.analysed, contrast=1, method=DBA_EDGER)
#     dba.plotMA(H3K36me3.analysed, contrast=1, method=DBA_DESEQ2,bNormalized=FALSE)
#     dba.plotMA(H3K36me3.analysed, contrast=1, method=DBA_DESEQ2)