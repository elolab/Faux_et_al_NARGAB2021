library(pROC)

library(GenomicRanges)
library(ChIPpeakAnno)
##############
#Data loading
##############
#setwd("/Users/thfaux/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/khgfgyt-1/wrk/asta/epouta/thomas_projects/B18070_Differential_ChIPseq_peak_calling/Global_Rerun/Results_comparison/Synthetic")
#This object contains the position of the simulated peaks
load("H3K36me3_bindingEvents_data.RData")
#setwd("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/khgfgyt-1/wrk/asta/epouta/thomas_projects/B18070_Differential_ChIPseq_peak_calling/Global_Rerun/Final_plots/ROC_curves")

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
  ol1 <- findOverlaps(query = grange, subject = bindings)
  
  #organize the results of the queries into dataframes
  dfol1 <- data.frame(queryHits(ol1),subjectHits(ol1),duplicated(subjectHits(ol1)),duplicated(queryHits(ol1)))
  
  # remove the duplicated peaks
  dfol1 <- dfol1[!dfol1$duplicated.subjectHits.ol1..,]
  
  #dfol1 <- dfol1[!dfol1$duplicated.queryHits.ol1..,]
  colnames(dfol1) <- c("prediction", "rank", "duplicatedquery","duplicatedsubject")
  
  #create the return object
  vec <- data.frame(rank=rep(1,20000),pvalue=rep(1,20000),FDR=rep(1,20000))
  
  #Put the rank, pvalue and FDR for each detected peak
  for(i in 1:dim(dfol1)[1]){
    vec[dfol1$rank[i],"rank"] <- dfol1$prediction[i] 
    vec[dfol1$rank[i],"pvalue"] <- temp[dfol1$prediction[i],"pvalue"]
    vec[dfol1$rank[i],"FDR"] <- temp[dfol1$prediction[i],"FDR"]
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
  vec <- vec[,c(1,2,3)]
  
  return(vec)
}



#load clean and call the main function to prepare the data for the pROC calls
temp <- read.table(paste0("Full_results/DiffBind_synthetic_full.bed"), header=T)
temp <- temp[,1:5]
temp <- cbind(temp,c(1:28647))
colnames(temp) <- c("space","start","end","pvalue","FDR", "rank") 
temp <- temp[order(temp$pvalue),]
#temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("DiffBind",temp)

temp <- read.table(paste0("Full_results/THOR_synthetic_full.bed"), header=T)
temp <- temp[,1:5]
temp <- cbind(temp,c(1:169399))
colnames(temp) <- c("space","start","end","pvalue","FDR", "rank") 
temp <- temp[order(temp$pvalue),]
#temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("THOR",temp)

temp <- read.table(paste0("Full_results/diffReps_synthetic_full.bed"), header=T)
temp <- temp[,1:5]
temp <- temp[order(temp$pvalue),]
temp <- cbind(temp,c(1:72450))
colnames(temp) <- c("space","start","end","pvalue","FDR", "rank") 
temp <- temp[order(temp$pvalue),]
#temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("diffReps",temp)

temp <- read.table(paste0("Full_results/ROTS_synthetic_full.bed"), header=T)
temp <- temp[,1:5]
temp <- cbind(temp,c(1:31563))
colnames(temp) <- c("space","start","end","pvalue","FDR", "rank") 
temp <- temp[order(temp$pvalue),]
#temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("ROTS",temp)

temp <- read.table(paste0("Full_results/PePr_synthetic_full.bed"), header=T)
temp <- temp[,1:5]
temp <- cbind(temp,c(1:656027))
colnames(temp) <- c("space","start","end","pvalue","FDR", "rank") 
temp <- temp[order(temp$pvalue),]
#temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("PePr",temp)


temp <- read.table(paste0("Full_results/MAnorm2_synthetic_full.bed"), header=T)
#temp <- temp[which(temp$V8 <= 0.05),]
#colnames(temp) <- c("space","start","end","name","score","strand","pval","qval","fc")
temp <- temp[order(temp$pval),]
temp <- temp[,1:5]
temp <- cbind(temp,c(1:32543))
colnames(temp) <- c("space","start","end","pvalue","FDR","rank") 
#temp <- toGRanges(temp, format="BED", header=FALSE) 
temp <- myfunction(temp)
assign("MAnorm2",temp)

truth <- c(rep(1,10000),rep(0,10000))




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

intensities <- c(rep(x = 10,1000),
                 rep(20,1000),
                 rep(30,1000),
                 rep(40,1000),
                 rep(50,1000),
                 rep(60,1000),
                 rep(70,1000),
                 rep(80,1000),
                 rep(90,1000),
                 rep(100,1000),
                 rep("FP",10000))
data <- data.frame(chr=seqnames(bindings),
                 starts=start(bindings)-1,
                 ends=end(bindings),
                 intensity=intensities,
                 peak_region_truth_status=truth,
                 ROTS.rank=ROTS$rank,
                 ROTS.pvalue=ROTS$pvalue,
                 ROTS.FDR=ROTS$FDR,
                 DiffBind.rank=DiffBind$rank,
                 DiffBind.pvalue=DiffBind$pvalue,
                 DiffBind.FDR=DiffBind$FDR,
                 MAnorm2.rank=MAnorm2$rank,
                 MAnorm2.pvalue=MAnorm2$pvalue,
                 MAnorm2.FDR=MAnorm2$FDR,
                 diffReps.rank=diffReps$rank,
                 diffReps.pvalue=diffReps$pvalue,
                 diffReps.FDR=diffReps$FDR,
                 PePr.rank=PePr$rank,
                 PePr.pvalue=PePr$pvalue,
                 PePr.FDR=PePr$FDR,
                 THOR.rank=THOR$rank,
                 THOR.pvalue=THOR$pvalue,
                 THOR.FDR=THOR$FDR)


write.table(x = data, 
            file = "Full_results/truth_table_full.bed", 
            quote = F , 
            col.names = T , 
            row.names = F , 
            sep = "\t" )

#############################################################
###             Supplementary figure 1                    ###
#############################################################

png(file="20212911_ROC_pvalues.png",height = 600,width = 600)
set.seed(1234)

roc(response=data$peak_region_truth_status, predictor=rank(data$DiffBind.pvalue, ties.method="random"), plot=T, transpose=F, direction=">", col="blue",cex.axis =2,cex.lab=2,xlab="",ylab="")

roc(response=data$peak_region_truth_status, predictor=rank(data$MAnorm2.pvalue, ties.method="random"), plot=T, transpose=F, direction=">", add=T, col="violet")

roc(response=data$peak_region_truth_status, predictor=rank(data$diffReps.pvalue, ties.method="random"), plot=T, transpose=F, direction=">", add=T, col="green")

roc(response=data$peak_region_truth_status, predictor=rank(data$PePr.pvalue, ties.method="random"), plot=T, transpose=F, direction=">", add=T, col="red")

roc(response=data$peak_region_truth_status, predictor=rank(data$THOR.pvalue, ties.method="random"), plot=T, transpose=F, direction=">", add=T, col="orange")

roc(response=data$peak_region_truth_status, predictor=rank(data$ROTS.pvalue, ties.method="random"), plot=T, transpose=F, direction=">", add=T, col="black")

legend("bottomright", c("ROTS","DiffBind","MAnorm2","diffReps","PePr","THOR"), fill=c("black","blue","violet","green","red","orange"), bty="n",cex=1.5)

dev.off()



#############################################################
###                  figure 1 A                           ###
#############################################################


data <- read.table("Full_results/truth_table_full.bed",header = TRUE,sep="\t")

# Count the peaks for each intensity of difference
ROTS <- c(dim(data[which(data[9001:10000,"ROTS.FDR"] <= 0.05),])[1],
          dim(data[which(data[8001:9000,"ROTS.FDR"] <= 0.05),])[1],
          dim(data[which(data[7001:8000,"ROTS.FDR"] <= 0.05),])[1],
          dim(data[which(data[6001:7000,"ROTS.FDR"] <= 0.05),])[1],
          dim(data[which(data[5001:6000,"ROTS.FDR"] <= 0.05),])[1],
          dim(data[which(data[4001:5000,"ROTS.FDR"] <= 0.05),])[1],
          dim(data[which(data[3001:4000,"ROTS.FDR"] <= 0.05),])[1],
          dim(data[which(data[2001:3000,"ROTS.FDR"] <= 0.05),])[1],
          dim(data[which(data[1001:2000,"ROTS.FDR"] <= 0.05),])[1],
          dim(data[which(data[1:1000,"ROTS.FDR"] <= 0.05),])[1],
          dim(data[which(data[10001:20000,"ROTS.FDR"] <= 0.05),])[1])



DiffBind <- c(dim(data[which(data[9001:10000,"DiffBind.FDR"] <= 0.05),])[1],
          dim(data[which(data[8001:9000,"DiffBind.FDR"] <= 0.05),])[1],
          dim(data[which(data[7001:8000,"DiffBind.FDR"] <= 0.05),])[1],
          dim(data[which(data[6001:7000,"DiffBind.FDR"] <= 0.05),])[1],
          dim(data[which(data[5001:6000,"DiffBind.FDR"] <= 0.05),])[1],
          dim(data[which(data[4001:5000,"DiffBind.FDR"] <= 0.05),])[1],
          dim(data[which(data[3001:4000,"DiffBind.FDR"] <= 0.05),])[1],
          dim(data[which(data[2001:3000,"DiffBind.FDR"] <= 0.05),])[1],
          dim(data[which(data[1001:2000,"DiffBind.FDR"] <= 0.05),])[1],
          dim(data[which(data[1:1000,"DiffBind.FDR"] <= 0.05),])[1],
          dim(data[which(data[10001:20000,"DiffBind.FDR"] <= 0.05),])[1])

MAnorm2 <- c(dim(data[which(data[9001:10000,"MAnorm2.FDR"] <= 0.05),])[1],
          dim(data[which(data[8001:9000,"MAnorm2.FDR"] <= 0.05),])[1],
          dim(data[which(data[7001:8000,"MAnorm2.FDR"] <= 0.05),])[1],
          dim(data[which(data[6001:7000,"MAnorm2.FDR"] <= 0.05),])[1],
          dim(data[which(data[5001:6000,"MAnorm2.FDR"] <= 0.05),])[1],
          dim(data[which(data[4001:5000,"MAnorm2.FDR"] <= 0.05),])[1],
          dim(data[which(data[3001:4000,"MAnorm2.FDR"] <= 0.05),])[1],
          dim(data[which(data[2001:3000,"MAnorm2.FDR"] <= 0.05),])[1],
          dim(data[which(data[1001:2000,"MAnorm2.FDR"] <= 0.05),])[1],
          dim(data[which(data[1:1000,"MAnorm2.FDR"] <= 0.05),])[1],
          dim(data[which(data[10001:20000,"MAnorm2.FDR"] <= 0.05),])[1])

diffReps <- c(dim(data[which(data[9001:10000,"diffReps.FDR"] <= 0.05),])[1],
          dim(data[which(data[8001:9000,"diffReps.FDR"] <= 0.05),])[1],
          dim(data[which(data[7001:8000,"diffReps.FDR"] <= 0.05),])[1],
          dim(data[which(data[6001:7000,"diffReps.FDR"] <= 0.05),])[1],
          dim(data[which(data[5001:6000,"diffReps.FDR"] <= 0.05),])[1],
          dim(data[which(data[4001:5000,"diffReps.FDR"] <= 0.05),])[1],
          dim(data[which(data[3001:4000,"diffReps.FDR"] <= 0.05),])[1],
          dim(data[which(data[2001:3000,"diffReps.FDR"] <= 0.05),])[1],
          dim(data[which(data[1001:2000,"diffReps.FDR"] <= 0.05),])[1],
          dim(data[which(data[1:1000,"diffReps.FDR"] <= 0.05),])[1],
          dim(data[which(data[10001:20000,"diffReps.FDR"] <= 0.05),])[1])

PePr <- c(dim(data[which(data[9001:10000,"PePr.FDR"] <= 0.05),])[1],
          dim(data[which(data[8001:9000,"PePr.FDR"] <= 0.05),])[1],
          dim(data[which(data[7001:8000,"PePr.FDR"] <= 0.05),])[1],
          dim(data[which(data[6001:7000,"PePr.FDR"] <= 0.05),])[1],
          dim(data[which(data[5001:6000,"PePr.FDR"] <= 0.05),])[1],
          dim(data[which(data[4001:5000,"PePr.FDR"] <= 0.05),])[1],
          dim(data[which(data[3001:4000,"PePr.FDR"] <= 0.05),])[1],
          dim(data[which(data[2001:3000,"PePr.FDR"] <= 0.05),])[1],
          dim(data[which(data[1001:2000,"PePr.FDR"] <= 0.05),])[1],
          dim(data[which(data[1:1000,"PePr.FDR"] <= 0.05),])[1],
          dim(data[which(data[10001:20000,"PePr.FDR"] <= 0.05),])[1])

THOR <- c(dim(data[which(data[9001:10000,"THOR.FDR"] <= 0.05),])[1],
          dim(data[which(data[8001:9000,"THOR.FDR"] <= 0.05),])[1],
          dim(data[which(data[7001:8000,"THOR.FDR"] <= 0.05),])[1],
          dim(data[which(data[6001:7000,"THOR.FDR"] <= 0.05),])[1],
          dim(data[which(data[5001:6000,"THOR.FDR"] <= 0.05),])[1],
          dim(data[which(data[4001:5000,"THOR.FDR"] <= 0.05),])[1],
          dim(data[which(data[3001:4000,"THOR.FDR"] <= 0.05),])[1],
          dim(data[which(data[2001:3000,"THOR.FDR"] <= 0.05),])[1],
          dim(data[which(data[1001:2000,"THOR.FDR"] <= 0.05),])[1],
          dim(data[which(data[1:1000,"THOR.FDR"] <= 0.05),])[1],
          dim(data[which(data[10001:20000,"THOR.FDR"] <= 0.05),])[1])


# Create the final object containing the counts
ideal <- c(1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,0)
mycounts<- cbind(ideal,
               ROTS,
               DiffBind,
               MAnorm2,
               diffReps,
               PePr,
               THOR
)

row.names(counts) <- c("100", "90", "80","70", "60", "50","40","30", "20","10","FP")
#row.names(counts) <- c("10", "20", "30","40", "50", "60","70","80", "90","100","FP")

#png("Figure_1_A.png", width = 1800,height = 1800)
par(cex.axis=2)
par(mar = c(5,5,5,5))
barplot(mycounts, ylim=c(0,11000),main="",width = c(rep(1,6)),
        ylab = "Number of significant differential peaks", col=c("#800000","#e6194B","#f58231","#ffe119","#bfef45","#3cb44b","#42d4f4","#4363d8","#911eb4","#f032e6","grey"),cex.lab=2,cex.axis =2)
legend(3.3,11000,"topright",legend=rownames(counts),ncol = 6, fill = c("#800000","#e6194B","#f58231","#ffe119","#bfef45","#3cb44b","#42d4f4","#4363d8","#911eb4","#f032e6","grey"),cex=2)
#dev.off()