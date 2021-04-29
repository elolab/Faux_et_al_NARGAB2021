############################################################################################################################################
##
## makepeakset.R -- make the consensus and overlap rate analysis
## 28 november 2018
## Thomas Faux
## Medical Bioinformatics Centre
############################################################################################################################################

## readcount -- read the counts for each peaks of each replicate and make a matrix of binding affinity

############################################################################################################################################
##
## readcount -- read the counts for each peaks of each replicate and make a matrix of binding affinity
##
## DB_object -- DB_object : a differential binding object created earlier
## nbthreads -- integer   : Number of threads allocated to each iteration of read counting
##
## returns a DB object containing the info in the csv files, the peaks for each sample, the consensus, the coverage and the matrix of binding affinity
############################################################################################################################################

readcount <- function(DB_object,nthreads=4){
  countMatrix <- Rsubread::featureCounts(files = DB_object$Info$V2,annot.ext = DB_object$Consensus,nthreads = nthreads,isPairedEnd=TRUE)
  return_object <- list(DB_object$Peaks,DB_object$Consensus,DB_object$Coverage,DB_object$Info,countMatrix$counts)
  names(return_object) <- c("Peaks","Consensus","Coverage","Info","Counts")
  return(return_object)
}

