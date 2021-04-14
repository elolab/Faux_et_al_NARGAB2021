fd############################################################################################################################################
##
## makepeakset.R -- make the consensus and overlap rate analysis
## 28 november 2018
## Thomas Faux
## Medical Bioinformatics Centre
############################################################################################################################################

## differentialCall -- performs the differential analysis with ROTS

## outputGeneration -- Generates a human readable output 

############################################################################################################################################
##
## differentialCall -- performs the differential analysis
##
## DB_object   -- DB_object   : a differential binding object created earlier
## B           -- integer     : Number of iteration for bootstrap
## K           -- integer     : Size of the top list to compare
## seed        -- integer     : seed for reproducibility
## paired      -- boolean     : performing paired analysis or not
## 
##
## returns a DB object containing the info in the csv files, the peaks for each sample, the consensus and the coverage, binding affinity matrix and normalized matrix
############################################################################################################################################
differentialCall <- function(DB_object,B,K,seed,paired,normalized=TRUE){
  groups <- as.numeric(DB_object$Info$V3)
  if(normalized){
    rots_out <- ROTS(DB_object$CountsNorm, groups, B=B, K=K, seed = seed,paired = paired, log = FALSE,progress = TRUE)
  }
  else{
    rots_out <- ROTS(DB_object$Counts, groups, B=B, K=K, seed = seed,paired = paired, log = FALSE,progress = TRUE)
  }
  
  return_object <- list(DB_object$Peaks,DB_object$Consensus,DB_object$Coverage,DB_object$Info,DB_object$Counts,DB_object$CountsNorm,rots_out)
  names(return_object) <- c("Peaks","Consensus","Coverage","Info","Counts","CountsNorm","ROTSout")
  return(return_object)
}

############################################################################################################################################
##
## outputGeneration -- Generates a human readable output 
##
## DB_object   -- DB_object   : a differential binding object created earlier
## fdr         -- double      : a threshold value for the output
## 
##
## returns a dataframe of the significant regions
############################################################################################################################################
outputGeneration <- function(DB_object, fdr){
  rots_summary<-summary.ROTS(DB_object$ROTSout,fdr=fdr)
  consensusConditions <- DB_object$Consensus
  #results <- cbind(data.frame(consensusConditions[rots_summary[,1],]),rots_summary[,2:4])
  #test the new result making, the old one was messing with the row order
  results <- cbind(consensusConditions[rownames(rots_summary),],rownames(rots_summary),rots_summary,DB_object$ROTSout$logfc[rownames(rots_summary)])
  colnames(results) <- c("GeneID", "Chr", "Start", "End", "width", "Strand","rownames", "Row", "ROTS-statistic", "pvalue", "FDR", "logfc")
  return(results)
}


