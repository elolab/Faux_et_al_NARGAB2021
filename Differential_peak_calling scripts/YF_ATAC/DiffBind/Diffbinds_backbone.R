library(DiffBind)
library(DESeq2)
library(edgeR)

sheet <- dba(sampleSheet="diffbinds_ATAC.csv")
print("data loading done")

olap.rate <- dba.overlap(sheet,mode=DBA_OLAP_RATE)
#jpeg("overlapRateATAC.jpeg")
#plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets')
#dev.off()

peakset <- dba.peakset(sheet, consensus=DBA_CONDITION)
#jpeg("peakOverlapATAC.jpeg")
#dba.plotVenn(peakset,peakset$masks$Consensus)
#dev.off()

ATAC.count <- dba.count(sheet)
save(ATAC.count,file="ATAC.count.RData")
print("count done")

ATAC.contrast <- dba.contrast(ATAC.count, categories=DBA_CONDITION)
save(ATAC.contrast,file="ATAC.contrast.RData")
print("contrast done")

ATAC.analysed <- dba.analyze(ATAC.contrast,bParallel=FALSE,method=c(DBA_DESEQ2,DBA_EDGER))
save(ATAC.analysed,file="ATAC.analysed_narrow.RData")
print("analyze done")

#report <- dba.report(ATAC.analysed, th=.05, DataType=DBA_DATA_FRAME,method=DBA_EDGER)
#score <- -10*(log10(report$FDR))
#write.table(cbind(report[,1:3],rownames(report),score),
#              "DBsites_YF_edger_narrow.bed", quote=FALSE, sep="\t",
#             row.names=FALSE, col.names=FALSE)

report <- dba.report(ATAC.analysed, th=1, DataType=DBA_DATA_FRAME,method=DBA_DESEQ2)
score <- report$FDR
FC <- report$Fold
write.table(cbind(report[,1:3],rownames(report),score,FC),
            "DBsites_YF_deseq2_narrow.bed", quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)
print("report done")