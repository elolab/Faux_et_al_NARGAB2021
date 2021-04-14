library(GenomicRanges)
library(ChIPpeakAnno)
##############
#Data loading
##############
files <- list.files(path="comparison_final_results/")
load("RData/H3K36me3_bindingEvents_data.RData")
myfunction <- function(grange){
  ####
  #### Counts the number of overlap for each downsampling level and return a vector of the values
  ####
  ol1 <- findOverlapsOfPeaks(grange,binding.events$`0.1sampledEvents`,connectedPeaks = "merge")
  ol2 <- findOverlapsOfPeaks(grange,binding.events$`0.2sampledEvents`,connectedPeaks = "merge")
  ol3 <- findOverlapsOfPeaks(grange,binding.events$`0.3sampledEvents`,connectedPeaks = "merge")
  ol4 <- findOverlapsOfPeaks(grange,binding.events$`0.4sampledEvents`,connectedPeaks = "merge")
  ol5 <- findOverlapsOfPeaks(grange,binding.events$`0.5sampledEvents`,connectedPeaks = "merge")
  ol6 <- findOverlapsOfPeaks(grange,binding.events$`0.6sampledEvents`,connectedPeaks = "merge")
  ol7 <- findOverlapsOfPeaks(grange,binding.events$`0.7sampledEvents`,connectedPeaks = "merge")
  ol8 <- findOverlapsOfPeaks(grange,binding.events$`0.8sampledEvents`,connectedPeaks = "merge")
  ol9 <- findOverlapsOfPeaks(grange,binding.events$`0.9sampledEvents`,connectedPeaks = "merge")
  ol10 <- findOverlapsOfPeaks(grange,binding.events$`1.0sampledEvents`,connectedPeaks = "merge")
  ol11 <- findOverlapsOfPeaks(grange,binding.events$commonBindingEvents,connectedPeaks = "merge")
  vec <- c(length(ol10$mergedPeaks),
           length(ol9$mergedPeaks),
           length(ol8$mergedPeaks),
           length(ol7$mergedPeaks),
           length(ol6$mergedPeaks),
           length(ol5$mergedPeaks),
           length(ol4$mergedPeaks),
           length(ol3$mergedPeaks),
           length(ol2$mergedPeaks),
           length(ol1$mergedPeaks),
           length(ol11$mergedPeaks))
  return(vec)
}

temp <- read.table(paste0("source/DBsitesH3K36me3_deseq2.bed"), header=F)
temp <- temp[,1:3]
colnames(temp) <- c("space","start","end") 
temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("DBdeseq2",temp)

temp <- read.table(paste0("../../Results_comparison/Synthetic/comparison_final_results/DBsitesH3K36me3_edger.bed"), header=F)
temp <- temp[,1:3]
colnames(temp) <- c("space","start","end") 
temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("DBedgeR",temp)

temp <- read.table(paste0("source/THOR-synthetic-diffpeaks.bed"), header=F)
temp <- temp[,1:3]
colnames(temp) <- c("space","start","end") 
temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("THOR",temp)

temp <- read.table(paste0("source/diffReps_simulated_nb.txt"), header=T)
temp <- temp[,1:3]
colnames(temp) <- c("space","start","end") 
temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("diffReps",temp)

temp <- read.table(paste0("source/myresults_RWV0_run.bed"), header=T)
temp <- temp[,2:4]
colnames(temp) <- c("space","start","end") 
temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("ROTS",temp)

temp <- read.table(paste0("source/Simulated__PePr_chip2_peaks.bed"), header=F)
temp <- temp[,1:3]
colnames(temp) <- c("space","start","end") 
temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("PePr",temp)


temp <- read.table(paste0("source/MAnorm2_nonorm.bed"), header=F)
temp <- temp[,1:3]
colnames(temp) <- c("space","start","end") 
temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("MAnorm2",temp)



ideal <- c(1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,0)
counts <- cbind(ideal,files)

DiffBind <- DBdeseq2
counts<- cbind(ideal,
               ROTS,
               DiffBind,
               MAnorm2,
               diffReps,
               PePr,
               THOR
               )

####### Short rerun
#load("RData/counts.RData")
#load("RData/counts_inversed")
row.names(counts) <- c("100", "90", "80","70", "60", "50","40","30", "20","10","FP")
#row.names(counts) <- c("10", "20", "30","40", "50", "60","70","80", "90","100","FP")

#pdf("Figure_1_DBEdgeR.pdf", width = 18,height = 18)
par(cex.axis=2)
par(mar = c(5,5,5,5))
barplot(counts, ylim=c(0,20000),main="",width = c(rep(1,6)),
        ylab = "Number of significant differential peaks", col=c("#800000","#e6194B","#f58231","#ffe119","#bfef45","#3cb44b","#42d4f4","#4363d8","#911eb4","#f032e6","grey"),cex.lab=2,cex.axis =2)
legend(5,15000,"topright",legend=rownames(counts),ncol = 3, fill = c("#800000","#e6194B","#f58231","#ffe119","#bfef45","#3cb44b","#42d4f4","#4363d8","#911eb4","#f032e6","grey"),cex=2)
#dev.off()

system.time()
