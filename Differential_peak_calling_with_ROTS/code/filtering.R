############################################################################################################################################
##
## filtering.R -- apply filtering on the matrix of counts stored in the DB object
## 28 november 2018
## Thomas Faux
## Medical Bioinformatics Centre
############################################################################################################################################

## filter -- Wrapper function that let you choose the type of filtering needed

## RWVfiltering-- apply Row Wise Variance filtering on the matrix of counts

## removeZeroLines -- Remove the zero lines contained in the matrix of counts

############################################################################################################################################
##
## RWVfiltering-- apply Row Wise Variance filtering on the matrix of counts and remove all the lines under a certain amount of variance
##
## matrix      -- matrix    : matrix of binding affinity
## threshold   -- integer   : Determine the threshold of filtering (amount of variance)
##
## returns a DB object containing the info in the csv files, the peaks for each sample, the consensus and the coverage
############################################################################################################################################
RWVfiltering <- function(matrix,threshold){
  size <- dim(matrix)[1]
  rowWiseVarience <- apply(matrix,MARGIN = 1,FUN = var)
  matrix <- matrix[which(rowWiseVarience >= threshold),]

  print(paste0("a total of ", size - dim(matrix)[1]," have been removed by Row wise variance filtering"))

  return(matrix)
}

############################################################################################################################################
##
## removeZeroLines -- Remove the zero lines contained in the matrix of counts
##
## matrix      -- matrix    : matrix of binding affinity
## threshold   -- boolean   : Determine the threshold of filtering (number of zeros among the samples)
##
## returns a DB object containing the info in the csv files, the peaks for each sample, the consensus and the coverage
############################################################################################################################################

removeZeroLines <- function(matrix, threshold){
  size <- dim(matrix)[1]
  numberZero <- apply(matrix,MARGIN = 1,FUN = function(x){length(which(x==0))})
  matrix <- matrix[-which(numberZero > threshold),]
  
  print(paste0("a total of ", size - dim(matrix)[1]," have been removed by removing the zero lines"))
  
  return(matrix)
}

############################################################################################################################################
##
## filter -- Wrapper function that let you choose the type of filtering needed
##
## DB_object   -- DB_object : created previously
## filtering   -- filtering : are available "RWV" (row wise variance) and "RZL" (remove zero line)
## threshold   -- boolean   : Determine the threshold of filtering (number of zeros among the samples or amount of variance)
##
## returns a DB object containing the info in the csv files, the peaks for each sample, the consensus and the coverage
############################################################################################################################################

filter <-function(DB_object,filtering,threshold){
  
  if(filtering %!in% c("RWV","RZL")){
    stop(" The type of filtering selected is not implemented, please select RWV or RZL")
  }
  switch(filtering,
         RWV = {
           if("CountsNorm" %in% names(DB_object)){
             matrix <- RWVfiltering(DB_object$CountsNorm,threshold)
           }
           else{
             matrix <- RWVfiltering(DB_object$Counts,threshold)
           }
         },
         RZL = {
           if("CountsNorm" %in% names(DB_object)){
             matrix <- removeZeroLines(DB_object$CountsNorm,threshold)
           }
           else{
             matrix <-removeZeroLines(DB_object$Counts,threshold)
           }
         }
  )
  
  return_object <- list(DB_object$Peaks,DB_object$Consensus,DB_object$Coverage,DB_object$Info,DB_object$Counts,matrix)
  names(return_object) <- c("Peaks","Consensus","Coverage","Info","Counts","CountsNorm")
  return(return_object)
}
