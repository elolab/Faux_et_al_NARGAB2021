library(ROTS)
library(GenomicRanges)
library(rtracklayer)
library(Rsubread)
library(Rsamtools)
library(ade4)
library(made4)

###################################
##
## load sources
##
###################################
source("code/differential_call.R")
source("code/readcount.R")
source("code/loaddata.R")
source("code/makepeakset.R")
source("code/Normalization.R")
source("code/filtering.R")
source("code/plots.R")
###################################
##
## Main
##
###################################
#Reading in the complete data
#Setting the filesPath and files variable for complete data

# To import narrowPeak files
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")

# To import broadPeak files
1#extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric")

Object_holder <- createDataObject("ROTS_input.csv",extraCols = extraCols_narrowPeak)

#create consensus peaks
Object_holder <- makePeakSet(DB_object = Object_holder,lower = 1)

#get the consensus to gtf format
Object_holder <- formatConsensus(Object_holder)

### create readers
Object_holder <- readcount(Object_holder)
#save(Object_holder, file = "RData/object.RData")

#rowWiseVarience <- apply(Object_holder$CountsNorm,MARGIN = 1,FUN = var)
#Object_holder <- filter(Object_holder,filtering = "RWV",rowWiseVarience)

normMethod <- "DESeq"
Object_holder <- normalization(Object_holder,normMethod=normMethod)

#MAplot(Object_holder,normalised=TRUE,"plots/MAplot_RLE_matrix_lower1_pfalse.png")
#PCAplot(Object_holder,normalised=TRUE,"plots/pca_RLE_data_lower1_pfalse.png")

  
print("ROTS call")
Object_holder2 <- differentialCall(Object_holder, B=100, K=floor(nrow(Object_holder$Filtered_NormCounts)/2), seed = 14,paired = FALSE, normalized = TRUE)
save(Object_holder2, file = "RData/object_lower1_pfalse.RData")
print("save results object")
#consensusConditions <- read.table("consensus.bed", header=TRUE, sep="\t")
results <-outputGeneration(Object_holder2,fdr=1)
write.table(results,"results/myresults_lower1_pfalse.bed",quote = FALSE,row.names = FALSE)
  
save(Object_holder2,file = "object_lower1_pfalse.RData")
#rots_summary <- summary(Object_holder2$ROTSout,fdr=0.05)
#pca(Object_holder2$ROTSout$data[rownames(rots_summary),],Object_holder2$Info$V3,"plots/pca_final_lower1_pfalse.png")
  
