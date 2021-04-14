library(pROC)

library(GenomicRanges)
library(ChIPpeakAnno)
##############
#Data loading
##############
setwd("/Users/thfaux/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/khgfgyt-1/wrk/asta/epouta/thomas_projects/B18070_Differential_ChIPseq_peak_calling/Global_Rerun/Results_comparison/Synthetic")
files <- list.files(path="comparison_final_results/")
load("H3K36me3_bindingEvents_data.RData")
myfunction <- function(grange){
  ####
  #### Counts the number of overlap for each downsampling level and return a vector of the values
  ####
  temp <- grange
  grange <- toGRanges(temp, format="BED", header=FALSE) 
  
#  ol1 <- findOverlaps(grange,binding.events$`0.1sampledEvents`)
#  ol2 <- findOverlaps(grange,binding.events$`0.2sampledEvents`)
#  ol3 <- findOverlaps(grange,binding.events$`0.3sampledEvents`)
#  ol4 <- findOverlaps(grange,binding.events$`0.4sampledEvents`)
#  ol5 <- findOverlaps(grange,binding.events$`0.5sampledEvents`)
#  ol6 <- findOverlaps(grange,binding.events$`0.6sampledEvents`)
#  ol7 <- findOverlaps(grange,binding.events$`0.7sampledEvents`)
#  ol8 <- findOverlaps(grange,binding.events$`0.8sampledEvents`)
#  ol9 <- findOverlaps(grange,binding.events$`0.9sampledEvents`)
#  ol10 <- findOverlaps(grange,binding.events$`1.0sampledEvents`)
#  ol11 <- findOverlaps(grange,binding.events$commonBindingEvents)
#  vec <- c(sort(unique(subjectHits(ol1))),
#           sort(unique(subjectHits(ol2)))+1000,
#           sort(unique(subjectHits(ol3)))+2000,
#           sort(unique(subjectHits(ol4)))+3000,
#           sort(unique(subjectHits(ol5)))+4000,
#           sort(unique(subjectHits(ol6)))+5000,
#           sort(unique(subjectHits(ol7)))+6000,
#           sort(unique(subjectHits(ol8)))+7000,
#           sort(unique(subjectHits(ol9)))+8000,
#           sort(unique(subjectHits(ol10)))+9000,
#           sort(unique(subjectHits(ol11)))+10000)

  ol1 <- findOverlaps(grange,binding.events$`0.1sampledEvents`)
  ol2 <- findOverlaps(grange,binding.events$`0.2sampledEvents`)
  ol3 <- findOverlaps(grange,binding.events$`0.3sampledEvents`)
  ol4 <- findOverlaps(grange,binding.events$`0.4sampledEvents`)
  ol5 <- findOverlaps(grange,binding.events$`0.5sampledEvents`)
  ol6 <- findOverlaps(grange,binding.events$`0.6sampledEvents`)
  ol7 <- findOverlaps(grange,binding.events$`0.7sampledEvents`)
  ol8 <- findOverlaps(grange,binding.events$`0.8sampledEvents`)
  ol9 <- findOverlaps(grange,binding.events$`0.9sampledEvents`)
  ol10 <- findOverlaps(grange,binding.events$`1.0sampledEvents`)
  ol11 <- findOverlaps(grange,binding.events$commonBindingEvents)
  
  dfol1 <- data.frame(queryHits(ol1),subjectHits(ol1),duplicated(subjectHits(ol1)))
  dfol2 <- data.frame(queryHits(ol2),subjectHits(ol2),duplicated(subjectHits(ol2)))
  dfol3 <- data.frame(queryHits(ol3),subjectHits(ol3),duplicated(subjectHits(ol3)))
  dfol4 <- data.frame(queryHits(ol4),subjectHits(ol4),duplicated(subjectHits(ol4)))
  dfol5 <- data.frame(queryHits(ol5),subjectHits(ol5),duplicated(subjectHits(ol5)))
  dfol6 <- data.frame(queryHits(ol6),subjectHits(ol6),duplicated(subjectHits(ol6)))
  dfol7 <- data.frame(queryHits(ol7),subjectHits(ol7),duplicated(subjectHits(ol7)))
  dfol8 <- data.frame(queryHits(ol8),subjectHits(ol8),duplicated(subjectHits(ol8)))
  dfol9 <- data.frame(queryHits(ol9),subjectHits(ol9),duplicated(subjectHits(ol9)))
  dfol10 <- data.frame(queryHits(ol10),subjectHits(ol10),duplicated(subjectHits(ol10)))
  dfol11 <- data.frame(queryHits(ol11),subjectHits(ol11),duplicated(subjectHits(ol11)))
  
  dfol1 <- dfol1[!dfol1$duplicated.subjectHits.ol1..,]
  dfol1$subjectHits.ol1. <- dfol1$subjectHits.ol1.+9000
  colnames(dfol1) <- c("prediction", "rank", "duplicated")
  
  dfol2 <- dfol2[!dfol2$duplicated.subjectHits.ol2..,]
  dfol2$subjectHits.ol2. <- dfol2$subjectHits.ol2.+8000
  colnames(dfol2) <- c("prediction", "rank", "duplicated")
  
  dfol3 <- dfol3[!dfol3$duplicated.subjectHits.ol3..,]
  dfol3$subjectHits.ol3. <- dfol3$subjectHits.ol3.+7000
  colnames(dfol3) <- c("prediction", "rank", "duplicated")
  
  dfol4 <- dfol4[!dfol4$duplicated.subjectHits.ol4..,]
  dfol4$subjectHits.ol4. <- dfol4$subjectHits.ol4.+6000
  colnames(dfol4) <- c("prediction", "rank", "duplicated")
  
  dfol5 <- dfol5[!dfol5$duplicated.subjectHits.ol5..,]
  dfol5$subjectHits.ol5. <- dfol5$subjectHits.ol5.+5000
  colnames(dfol5) <- c("prediction", "rank", "duplicated")
  
  dfol6 <- dfol6[!dfol6$duplicated.subjectHits.ol6..,]
  dfol6$subjectHits.ol6. <- dfol6$subjectHits.ol6.+4000
  colnames(dfol6) <- c("prediction", "rank", "duplicated")
  
  dfol7 <- dfol7[!dfol7$duplicated.subjectHits.ol7..,]
  dfol7$subjectHits.ol7. <- dfol7$subjectHits.ol7.+3000
  colnames(dfol7) <- c("prediction", "rank", "duplicated")
  
  dfol8 <- dfol8[!dfol8$duplicated.subjectHits.ol8..,]
  dfol8$subjectHits.ol8. <- dfol8$subjectHits.ol8.+2000
  colnames(dfol8) <- c("prediction", "rank", "duplicated")
  
  dfol9 <- dfol9[!dfol9$duplicated.subjectHits.ol9..,]
  dfol9$subjectHits.ol9. <- dfol9$subjectHits.ol9.+1000
  colnames(dfol9) <- c("prediction", "rank", "duplicated")
  
  dfol10 <- dfol10[!dfol10$duplicated.subjectHits.ol10..,]
  colnames(dfol10) <- c("prediction", "rank", "duplicated")
  
  dfol11 <- dfol11[!dfol11$duplicated.subjectHits.ol11..,]
  dfol11$subjectHits.ol11. <- dfol11$subjectHits.ol11.+10000
  colnames(dfol11) <- c("prediction", "rank", "duplicated")
  
  dfol <- rbind(dfol10,dfol9,dfol8,dfol7,dfol6,dfol5,dfol4,dfol3,dfol2,dfol1,dfol11)
  dfol$duplicated <- c(1:dim(dfol)[1])
  vec <- rep(0,20000)
  
  for(i in 1:dim(dfol)[1]){
    vec[dfol$rank[i]] <- dfol$prediction[i] 
  }
  
  
  
  
  #vec <- c(dfol10$queryHits.ol10.,
  #         dfol9$queryHits.ol9.,
  #         dfol8$queryHits.ol8.,
  #         dfol7$queryHits.ol7.,
  #         dfol6$queryHits.ol6.,
  #         dfol5$queryHits.ol5.,
  #        dfol4$queryHits.ol4.,
  #         dfol3$queryHits.ol3.,
  ##         dfol2$queryHits.ol2.,
  #         dfol1$queryHits.ol1.,
  #         dfol11$queryHits.ol1.
  #)  
  
    return(vec)
}

temp <- read.table(paste0("comparison_final_results/DBsitesH3K36me3_deseq2.bed"), header=F)
temp <- temp[,1:3]
temp <- cbind(temp,c(1:7965))
colnames(temp) <- c("space","start","end", "rank") 
#temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("DBdeseq2",temp)

temp <- read.table(paste0("comparison_final_results/DBsitesH3K36me3_edger.bed"), header=F)
temp <- temp[,1:3]
temp <- cbind(temp,c(1:10332))
colnames(temp) <- c("space","start","end", "rank")  
#temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("DBedger",temp)

temp <- read.table(paste0("comparison_final_results/THOR/THOR-synthetic-diffpeaks.bed"), header=F)
temp <- temp[,1:3]
temp <- cbind(temp,c(1:17597))
colnames(temp) <- c("space","start","end", "rank") 
#temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("THOR",temp)

temp <- read.table(paste0("comparison_final_results/diffReps_simulated_nb.txt"), header=T)
temp <- temp[,1:3]
temp <- cbind(temp,c(1:9743))
colnames(temp) <- c("space","start","end", "rank") 
#temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("diffReps",temp)

temp <- read.table(paste0("comparison_final_results/myresults_RWV0_run.bed"), header=T)
temp <- temp[,2:4]
temp <- cbind(temp,c(1:11112))
colnames(temp) <- c("space","start","end", "rank") 
#temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("ROTS",temp)

temp <- read.table(paste0("comparison_final_results/Simulated__PePr_chip2_peaks.bed"), header=F)
temp <- temp[,1:3]
temp <- cbind(temp,c(1:6660))
colnames(temp) <- c("space","start","end", "rank") 
#temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("PePr",temp)


temp <- read.table(paste0("comparison_final_results/MAnorm2_nonorm.bed"), header=F)
colnames(temp) <- c("space","start","end","name","score","strand","pval","qval","fc")
temp <- temp[order(temp$pval),]
temp <- temp[,1:3]
temp <- cbind(temp,c(1:32543))
colnames(temp) <- c("space","start","end","rank") 
#temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("MAnorm2",temp)

truth <- c(rep(1,10000),rep(0,10000))
#predictions_PePr <- c(rep(0,20000))
#predictions_PePr[PePr] <- 1
#predictions_THOR <- c(rep(0,20000))
#predictions_THOR[THOR] <- 1
#predictions_diffReps <- c(rep(0,20000))
#predictions_diffReps[diffReps] <- 1
#predictions_DBdeseq2 <- c(rep(0,20000))
#predictions_DBdeseq2[DBdeseq2] <- 1
#predictions_DBedger <- c(rep(0,20000))
#predictions_DBedger[DBedger] <- 1
#predictions_ROTS <- c(rep(0,20000))
#predictions_ROTS[ROTS] <- 1


results <- pROC::roc(truth~ROTS,percent=TRUE)
pROC::plot.roc(results, print.auc=TRUE, col="black",lwd=2,lty=6,print.auc.y = .9, print.auc.x = .85, direction = "<", smooth=FALSE)
pROC::plot.roc(truth,  THOR,percent=TRUE,print.auc=TRUE, col="orange",lwd=2,lty=6,print.auc.y = .8, print.auc.x = .9,add = TRUE, direction = "<", smooth=FALSE)
pROC::plot.roc(truth,  PePr,percent=TRUE,print.auc=TRUE, col="red",lwd=2,lty=1,print.auc.y = .675, print.auc.x = .8,add = TRUE, direction = "<", smooth=FALSE)
pROC::plot.roc(truth,  diffReps,percent=TRUE,print.auc=TRUE,col="green",lwd=2,lty=6,print.auc.y = .55, print.auc.x = .95,add = TRUE, direction = "<", smooth=FALSE)
pROC::plot.roc(truth,  DBdeseq2,percent=TRUE,print.auc=TRUE,col="blue",lwd=2,lty=2,print.auc.y = .7, print.auc.x = .8,add = TRUE, direction = "<", smooth=FALSE)
pROC::plot.roc(truth,  MAnorm2,percent=TRUE,print.auc=TRUE,col="violet",lwd=2,lty=2,print.auc.y = .775, print.auc.x = .95,add = TRUE, direction = "<", smooth=FALSE)

legend("bottomright", legend=c("ROTS", "DiffBind_DEseq2","MAnorm2","diffReps","PePr","THOR"), col=c("black", "blue","violet","green","red","orange"), lwd=2)

resultsROTS <- pROC::roc(truth~ROTS,percent=TRUE)
resultsTHOR <- pROC::roc(truth~THOR,percent=TRUE)
resultsPePr <- pROC::roc(truth~PePr,percent=TRUE)
resultsDBdeseq2 <- pROC::roc(truth~DBdeseq2,percent=TRUE)
resultsMAnorm2 <- pROC::roc(truth~MAnorm2,percent=TRUE)
resultsdiffReps <- pROC::roc(truth~diffReps,percent=TRUE)
plot(precision ~ recall,
     coords(resultsROTS, "all", ret = c("recall", "precision"), transpose = FALSE),
     type="l",lwd=3,lty=1, ylim = c(50, 100))
lines(precision ~ recall,
     coords(resultsTHOR, "all", ret = c("recall", "precision"), transpose = FALSE),
     col="orange",lwd=3,lty=1, ylim = c(50, 100))
lines(precision ~ recall,
      coords(resultsDBdeseq2, "all", ret = c("recall", "precision"), transpose = FALSE),
      col="blue",lwd=3,lty=1, ylim = c(50, 100))
lines(precision ~ recall,
      coords(resultsdiffReps, "all", ret = c("recall", "precision"), transpose = FALSE),
      col="green",lwd=3,lty=1, ylim = c(50, 100))
lines(precision ~ recall,
      coords(resultsMAnorm2, "all", ret = c("recall", "precision"), transpose = FALSE),
      col="violet",lwd=3,lty=1, ylim = c(50, 100))
lines(precision ~ recall,
      coords(resultsPePr, "all", ret = c("recall", "precision"), transpose = FALSE),
      col="red",lwd=3,lty=3, ylim = c(50, 100))

