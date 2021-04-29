############################################################################################################################################
##
## loaddata.R -- create a data object with the input files given
## 28 november 2018
## Thomas Faux
## Medical Bioinformatics Centre
############################################################################################################################################

## createDataObject -- construct the differential binding analysis initial object

## %!in% -- Utility function inverse of %in%

############################################################################################################################################
##
## createDataObject -- construct the differential binding analysis initial object
##
## path      -- character          : path to the csv file containing files informations
## extraCols -- list of characters : needed to import bed files with more than 4 columns
##           -- exemple for MACS2 narrow peaks output c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")
##
## returns a list object containing the info in the csv files and the peaks for each sample
############################################################################################################################################


createDataObject <- function(path,extraCols=NULL){
  files <- read.table(file = "ROTS_input.csv",sep = ";",colClasses = c("character","character","character"))
  # reads the input csv file and create the necessary objects for running further steps
  {
    completeDATA <- lapply(files$V1,import,format = "BED",extraCols = extraCols)
    names(completeDATA) <- basename(path = files$V1)
  }
  return_object <- list(GRangesList(completeDATA),files)
  names(return_object) <- c("Peaks","Info")
  return(return_object)
}


############################################################################################################################################
##
## %!in% -- Utility function inverse of %in%
##
## x      -- character
## y      -- character
##           
############################################################################################################################################
{
  '%!in%' <- function(x,y)!('%in%'(x,y)) 
}



