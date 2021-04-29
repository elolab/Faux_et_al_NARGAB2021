############################################################################################################################################
##
## makepeakset.R -- make the consensus and overlap rate analysis
## 28 november 2018
## Thomas Faux
## Medical Bioinformatics Centre
############################################################################################################################################

## makePeakSet -- create the consensus by overlaping the peaks of all the samples

## formatConsensus -- convert the consensus to GTF format

############################################################################################################################################
##
## makePeakSet -- create the consensus by overlaping the peaks of all the samples
##
## DB_object   -- DB_object : a differential binding object created earlier
## peak.center -- boolean   : Determine if the peaks will be recentered or not
## peak.ext    -- integer   : nb of bases to add upstream and downstream after recentering (determine the size of the peaks)
## lower       -- integer   : keep the peaks that overlap in at least this many samples
## merge       -- integer   : minimum gap width for merging of the peaks, if 0 then no merging done
## 
##
## returns a DB object containing the info in the csv files, the peaks for each sample, the consensus and the coverage
############################################################################################################################################

makePeakSet<-function(DB_object,peak.center=FALSE,peak.ext=0,lower=1, merge=0){
  peaks <- DB_object$Peaks
  if(peak.center){
    for(i in 1:length(peaks)){
      start(peaks[[i]])<-floor((start(peaks[[i]])+end(peaks[[i]]))/2)-peak.ext
      end(peaks[[i]])<-floor((start(peaks[[i]])+end(peaks[[i]]))/2)+peak.ext
      start(peaks[[i]])[start(peaks[[i]])<=0] <- 1
    }
  }
  peak_coverage <- coverage(peaks)
  coverage <- unlist(sapply(1:length(peaks),simplify = F, function(f) {length(GRanges(slice(peak_coverage, lower=f, rangesOnly=T)))}))
  covered_ranges <- slice(peak_coverage, lower=lower, rangesOnly=T)
  covered_granges <- GRanges(covered_ranges)
  if(merge > 0){
    reduce(covered_granges, min.gapwidth=merge,drop.empty.ranges = TRUE)
  }
  
  return_object <- list(peaks,covered_granges,coverage,DB_object$Info)
  names(return_object) <- c("Peaks","Consensus","Coverage","Info")

  return(return_object)
}

############################################################################################################################################
##
## formatConsensus -- convert the consensus to GTF format
##
## DB_object   -- DB_object : a differential binding object created earlier
##
## returns the DB object with the consensus in GTF format
############################################################################################################################################

formatConsensus <- function(DB_object){
  df <- as.data.frame(DB_object$Consensus)
  prefix <- "peak"
  suffix <- seq(1:length(DB_object$Consensus))
  peaknames <- paste(prefix,suffix,sep = "_")
  width<- end(DB_object$Consensus)-start(DB_object$Consensus)
  consensus <- cbind(peaknames,df[,c(1,2,3)],width,df[,5])
  rownames(consensus) <- consensus[,1]
  colnames(consensus) <-c("GeneID","Chr","Start","End","width","Strand")
  
  return_object <- list(DB_object$Peaks,consensus,DB_object$Coverage,DB_object$Info)
  names(return_object) <- c("Peaks","Consensus","Coverage","Info")
  
  return(return_object)
  
}


#old function that center the peaks on peak summit instead of peak center
#should be implemented latter on

# makePeakSetTemp<-function(peaks,peak.center=FALSE,peak.ext=0){
#   
#   message("Making peak list......\n")
#   n=length(peaks)
#   peak.list=GRangesList()
#   pmat=import(peaks[1])
#   if(peak.center){
#     start(pmat)=(start(pmat)+end(pmat))/2-peak.ext
#     end(pmat)=(start(pmat)+end(pmat))/2+peak.ext
#     start(pmat)[start(pmat)<=0]=1
#   }	
#   peak.list[[1]]=pmat
#   for(i in 2:n){
#     mat=import(peaks[i])
#     peak.list[[i]]=mat
#     if(peak.center){
#       start(pmat)=(start(pmat)+end(pmat))/2-peak.ext
#       end(pmat)=(start(pmat)+end(pmat))/2+peak.ext
#       start(pmat)[start(pmat)<=0]=1
#     }
#     pmat=union(pmat,mat)
#   }
#   tmp=findOverlaps(peak.list[[1]],pmat)
#   oidx=unique(subjectHits(tmp))	
#   for(i in 2:n){
#     tmp=findOverlaps(peak.list[[i]],pmat)
#     oidx=intersect(oidx, unique(subjectHits(tmp)))	
#   }
#   peakSet=list(peak.list=peak.list,pmat=pmat,oidx=oidx)
#   peakSet
# }

# makePeakSetCenterSummit<-function(peaks,peak.center=FALSE,peak.ext=0){
#   
#   message("Making peak list......\n")
#   n=length(peaks)
#   peak.list=GRangesList()
#   pmat=import(peaks[1])
#   summits <- read.table(paste0(strsplit(peaks[1],split=".bed"),".narrowPeak"))
#   summits <- summits[,9]
#   pmat$summits <- summits
#   
#   if(peak.center){
#     start(pmat)=floor((start(pmat)+pmat$summits)-peak.ext)
#     end(pmat)=floor((start(pmat)+pmat$summits)+peak.ext)
#     start(pmat)[start(pmat)<=0]=1
#   } 
#   peak.list[[1]]=pmat
#   for(i in 2:n){
#     mat=import(peaks[i])
#     summits <- read.table(paste0(strsplit(peaks[i],split=".bed"),".narrowPeak"))
#     summits <- summits[,9]
#     mat$summits <- summits
#     peak.list[[i]]=mat
#     
#     if(peak.center){
#       start(mat)=floor((start(mat)+mat$summits)-peak.ext)
#       end(mat)=floor((start(mat)+mat$summits)+peak.ext)
#       start(mat)[start(mat)<=0]=1
#     }
#     pmat=union(pmat,mat)
#   }
#   tmp=findOverlaps(peak.list[[1]],pmat)
#   oidx=unique(subjectHits(tmp)) 
#   for(i in 2:n){
#     tmp=findOverlaps(peak.list[[i]],pmat)
#     oidx=intersect(oidx, unique(subjectHits(tmp)))  
#   }
#   peakSet=list(peak.list=peak.list,pmat=pmat,oidx=oidx)
#   peakSet
# }

