library(pROC)

library(GenomicRanges)
library(ChIPpeakAnno)
##############
#Data loading
##############
setwd("/Users/thfaux/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/khgfgyt-1/wrk/asta/epouta/thomas_projects/B18070_Differential_ChIPseq_peak_calling/Global_Rerun/Results_comparison/Synthetic")
load("H3K36me3_bindingEvents_data.RData")
setwd("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/khgfgyt-1/wrk/asta/epouta/thomas_projects/B18070_Differential_ChIPseq_peak_calling/Global_Rerun/Final_plots/ROC_curves")

myfunction <- function(grange){
  ####
  #### Counts the number of overlap for each downsampling level and return a vector of the values
  ####
  temp <- grange
  grange <- toGRanges(temp, format="BED", header=FALSE) 
  
  #combines the binding events in one GRanges object
  bindings <- do.call("c",List(binding.events$`0.1sampledEvents`,
       binding.events$`0.2sampledEvents`,
       binding.events$`0.3sampledEvents`,
       binding.events$`0.4sampledEvents`,
       binding.events$`0.5sampledEvents`,
       binding.events$`0.6sampledEvents`,
       binding.events$`0.7sampledEvents`,
       binding.events$`0.8sampledEvents`,
       binding.events$`0.9sampledEvents`,
       binding.events$`1.0sampledEvents`,
       binding.events$commonBindingEvents))

  #finds the overlap between the given regions and the regions contained in the truth
  ol1 <- findOverlaps(grange,bindings)
  
  #organize the results of the queries into dataframes
  dfol1 <- data.frame(queryHits(ol1),subjectHits(ol1),duplicated(subjectHits(ol1)))
  
  # remove the duplicated peaks
  dfol1 <- dfol1[!dfol1$duplicated.subjectHits.ol1..,]
  
  # replces the binary values of the duplicated column by numbers
  dfol1$duplicated <- c(1:dim(dfol1)[1])
  colnames(dfol1) <- c("prediction", "rank", "duplicated")
  
  #create the return object
  vec <- data.frame(rank=rep(1,20000),pvalue=rep(1,20000))
  
  #Put the rank and pvalue for each detected peak
  for(i in 1:dim(dfol1)[1]){
    vec[dfol1$rank[i],"rank"] <- dfol1$prediction[i] 
    vec[dfol1$rank[i],"pvalue"] <- temp[dfol1$prediction[i],"pvalue"]
    }
  
  #Put a high rank to the peaks that are not found

  for(i in 1:dim(vec)[1]){
    if(vec$rank[i] ==1 && vec$pvalue[i] ==1){
      vec[i,"rank"] <- max(vec$rank) + 1
    }
      
  }
  
  #recalculate the rank to get 1 to 20000
  vec$trueorder <- 1:20000
  vec <- vec[order(vec$rank),]
  vec$rank <- 1:20000
  vec <- vec[order(vec$trueorder),]
  vec <- vec[,c(1,2)]
  
  return(vec)
}



#load clean and call the main function to prepare the data for the pROC calls
temp <- read.table(paste0("Full_results/DiffBind_synthetic_full.bed"), header=T)
temp <- temp[,1:4]
temp <- cbind(temp,c(1:28647))
colnames(temp) <- c("space","start","end","pvalue", "rank") 
#temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("DiffBind",temp)

temp <- read.table(paste0("Full_results/THOR_synthetic_full.bed"), header=T)
temp <- temp[,1:4]
temp <- cbind(temp,c(1:169399))
colnames(temp) <- c("space","start","end","pvalue", "rank") 
#temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("THOR",temp)

temp <- read.table(paste0("Full_results/diffReps_synthetic_full.bed"), header=T)
temp <- temp[,1:4]
temp <- cbind(temp,c(1:72450))
colnames(temp) <- c("space","start","end","pvalue", "rank") 
#temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("diffReps",temp)

temp <- read.table(paste0("Full_results/ROTS_synthetic_full.bed"), header=T)
temp <- temp[,1:4]
temp <- cbind(temp,c(1:31563))
colnames(temp) <- c("space","start","end","pvalue", "rank") 
#temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("ROTS",temp)

temp <- read.table(paste0("Full_results/PePr_synthetic_full.bed"), header=T)
temp <- temp[,1:4]

temp <- cbind(temp,c(1:656027))
colnames(temp) <- c("space","start","end","pvalue", "rank") 
temp <- temp[order(temp$pvalue),]
#temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("PePr",temp)


temp <- read.table(paste0("Full_results/MAnorm2_synthetic_full.bed"), header=T)
#temp <- temp[which(temp$V8 <= 0.05),]
#colnames(temp) <- c("space","start","end","name","score","strand","pval","qval","fc")
#temp <- temp[order(temp$pval),]
temp <- temp[,1:4]
temp <- cbind(temp,c(1:32543))
colnames(temp) <- c("space","start","end","pvalue","rank") 
#temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("MAnorm2",temp)

truth <- c(rep(1,10000),rep(0,10000))


##ROC curve call
results <- pROC::roc(truth~ROTS$rank,percent=TRUE)
pROC::plot.roc(results, print.auc=TRUE, col="black",lwd=2,lty=6,print.auc.y = .9, print.auc.x = .85, direction = "<", smooth=FALSE)
pROC::plot.roc(truth,  THOR$rank,percent=TRUE,print.auc=TRUE, col="orange",lwd=2,lty=6,print.auc.y = .8, print.auc.x = .9,add = TRUE, direction = "<", smooth=FALSE)
pROC::plot.roc(truth,  PePr$rank,percent=TRUE,print.auc=TRUE, col="red",lwd=2,lty=1,print.auc.y = .675, print.auc.x = .8,add = TRUE, direction = "<", smooth=FALSE)
pROC::plot.roc(truth,  diffReps$rank,percent=TRUE,print.auc=TRUE,col="green",lwd=2,lty=6,print.auc.y = .55, print.auc.x = .95,add = TRUE, direction = "<", smooth=FALSE)
pROC::plot.roc(truth,  DiffBind$rank,percent=TRUE,print.auc=TRUE,col="blue",lwd=2,lty=2,print.auc.y = .7, print.auc.x = .8,add = TRUE, direction = "<", smooth=FALSE)
pROC::plot.roc(truth,  MAnorm2$rank,percent=TRUE,print.auc=TRUE,col="violet",lwd=2,lty=2,print.auc.y = .775, print.auc.x = .95,add = TRUE, direction = "<", smooth=FALSE)

legend("bottomright", legend=c("ROTS", "DiffBind_DEseq2","MAnorm2","diffReps","PePr","THOR"), col=c("black", "blue","violet","green","red","orange"), lwd=2)

results <- pROC::roc(truth~ROTS$pvalue,percent=TRUE)
pROC::plot.roc(results, print.auc=TRUE, col="black",lwd=2,lty=6,print.auc.y = .9, print.auc.x = .85, direction = "<", smooth=FALSE)
pROC::plot.roc(truth,  THOR$pvalue,percent=TRUE,print.auc=TRUE, col="orange",lwd=2,lty=6,print.auc.y = .8, print.auc.x = .9,add = TRUE, direction = "<", smooth=FALSE)
pROC::plot.roc(truth,  PePr$pvalue,percent=TRUE,print.auc=TRUE, col="red",lwd=2,lty=1,print.auc.y = .675, print.auc.x = .8,add = TRUE, direction = "<", smooth=FALSE)
pROC::plot.roc(truth,  diffReps$pvalue,percent=TRUE,print.auc=TRUE,col="green",lwd=2,lty=6,print.auc.y = .55, print.auc.x = .95,add = TRUE, direction = "<", smooth=FALSE)
pROC::plot.roc(truth,  DiffBind$pvalue,percent=TRUE,print.auc=TRUE,col="blue",lwd=2,lty=2,print.auc.y = .7, print.auc.x = .8,add = TRUE, direction = "<", smooth=FALSE)
pROC::plot.roc(truth,  MAnorm2$pvalue,percent=TRUE,print.auc=TRUE,col="violet",lwd=2,lty=2,print.auc.y = .775, print.auc.x = .95,add = TRUE, direction = "<", smooth=FALSE)


#precision vs recall call
resultsROTS <- pROC::roc(truth~ROTS$rank,percent=TRUE)
resultsTHOR <- pROC::roc(truth~THOR$rank,percent=TRUE)
resultsPePr <- pROC::roc(truth~PePr$rank,percent=TRUE)
resultsDiffBind <- pROC::roc(truth~DiffBind$rank,percent=TRUE)
resultsMAnorm2 <- pROC::roc(truth~MAnorm2$rank,percent=TRUE)
resultsdiffReps <- pROC::roc(truth~diffReps$rank,percent=TRUE)
par(mar = c(5,5,5,5))
plot(precision ~ recall,
     coords(resultsROTS, "all", ret = c("recall", "precision"), transpose = FALSE),
     type="l",lwd=3,lty=1, ylim = c(50, 100),cex.lab=2,cex.axis =2)
lines(precision ~ recall,
     coords(resultsTHOR, "all", ret = c("recall", "precision"), transpose = FALSE),
     col="orange",lwd=3,lty=1, ylim = c(50, 100))
lines(precision ~ recall,
      coords(resultsDiffBind, "all", ret = c("recall", "precision"), transpose = FALSE),
      col="blue",lwd=3,lty=1, ylim = c(50, 100))
lines(precision ~ recall,
      coords(resultsdiffReps, "all", ret = c("recall", "precision"), transpose = FALSE),
      col="green",lwd=3,lty=1, ylim = c(50, 100))
lines(precision ~ recall,
      coords(resultsMAnorm2, "all", ret = c("recall", "precision"), transpose = FALSE),
      col="violet",lwd=3,lty=1, ylim = c(50, 100))
lines(precision ~ recall,
      coords(resultsPePr, "all", ret = c("recall", "precision"), transpose = FALSE),
      col="red",lwd=3,lty=1, ylim = c(50, 100))


#write the pROC input files in one file
bindings <- do.call("c",List(binding.events$`0.1sampledEvents`,
                             binding.events$`0.2sampledEvents`,
                             binding.events$`0.3sampledEvents`,
                             binding.events$`0.4sampledEvents`,
                             binding.events$`0.5sampledEvents`,
                             binding.events$`0.6sampledEvents`,
                             binding.events$`0.7sampledEvents`,
                             binding.events$`0.8sampledEvents`,
                             binding.events$`0.9sampledEvents`,
                             binding.events$`1.0sampledEvents`,
                             binding.events$commonBindingEvents))

df <- data.frame(chr=seqnames(bindings),
                 starts=start(bindings)-1,
                 ends=end(bindings),
                 peak_region_truth_status=truth,
                 ROTS.rank=ROTS$rank,
                 ROTS.pvalue=ROTS$pvalue,
                 DiffBind.rank=DiffBind$rank,
                 DiffBind.pvalue=DiffBind$pvalue,
                 MAnorm2.rank=MAnorm2$rank,
                 MAnorm2.pvalue=MAnorm2$pvalue,
                 diffReps.rank=diffReps$rank,
                 diffReps.pvalue=diffReps$pvalue,
                 PePr.rank=PePr$rank,
                 PePr.pvalue=PePr$pvalue,
                 THOR.rank=THOR$rank,
                 THOR.pvalue=THOR$pvalue)

                 
                 
write.table(x = df, 
            file = "truth_table_full.bed", 
            quote = F , 
            col.names = T , 
            row.names = F , 
            sep = "\t" )
summary(df)
