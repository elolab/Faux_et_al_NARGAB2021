############################################################################################################################################
##
## plots.R -- make the consensus and overlap rate analysis
## 28 november 2018
## Thomas Faux
## Medical Bioinformatics Centre
############################################################################################################################################

## MAplot -- plots the MA plot of the data contained in the binding affinity matrix 

## PCAplot -- plots the PCA of the data contained in the binding affinity matrix 

## OverlapRateplot -- plot the overlap rate 

############################################################################################################################################
##
## MAplot -- plots an MAplot of the values contained in a matrix
##
## DB_object    -- DB_object   : a differential binding object created earlier
## normalised   -- boolean     : Do you want to plot the normalised values or not
## name         -- character   : name of the plot
##
############################################################################################################################################

MAplot <- function(DB_object, normalised=TRUE,name){
  
  groups <- DB_object$Info$V3
  
  if(normalised == TRUE){
    matrix <- DB_object$CountsNorm+1
  }
  else{
    matrix <- DB_object$Counts+1
  }
  
  png(height = 900,width = 900,filename =  name)
  ma.plot( rowMeans(log2(matrix)), log2(apply(matrix[,which(groups==1)],1,mean))-log2(apply(matrix[,which(groups==2)],1,mean)), cex=1 )
  dev.off()

}

############################################################################################################################################
##
## PCAplot -- plots the PCA of the data contained in the binding affinity matrix 
##
## DB_object    -- DB_object   : a differential binding object created earlier
## normalised   -- boolean     : Do you want to plot the normalised values or not
## name         -- character   : name of the plot
##
############################################################################################################################################

PCAplot <- function(DB_object, normalised=TRUE,name){
  
  if(normalised == TRUE){
    matrix <- DB_object$CountsNorm
  }
  else{
    matrix <- DB_object$Counts
  }
  groups <- DB_object$Info$V3
  pca <- dudi.pca(matrix,scan=F,nf=3)
  
  png(height = 900,width = 900,filename =  name)
  rotate3d(pca$co, classvec=groups)
  dev.off()
  
}

############################################################################################################################################
##
## plotOverlapRate -- plot the overlap rate 
##
## DB_object   -- DB_object : a differential binding object created earlier
## pos         -- integer   : a position specifier for the text. If specified this overrides any adj value given.
##                            Values of 1, 2, 3 and 4, respectively indicate positions below, to the left of, 
##                            above and to the right of the specified (x,y) coordinates.
## cex         -- double    : numeric character expansion factor; multiplied by par("cex") yields the final character size. NULL and NA are equivalent to 1.0.
############################################################################################################################################

OverlapRateplot <- function(DB_object,pos=3,cex=0.7){
  
  plot(DB_object$Coverage,type="o", ylab = "# peaks", xlab = "Overlap in at least this many samples")
  text(DB_object$Coverage,cex = cex ,pos=pos,labels = DB_object$Coverage)
  
}


