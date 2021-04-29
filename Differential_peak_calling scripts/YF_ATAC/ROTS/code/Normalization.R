############################################################################################################################################
##
## normalization.R -- choose and apply normalizations
## 28 november 2018
## Thomas Faux
## Medical Bioinformatics Centre
############################################################################################################################################

## normalization -- wrapper function for normalisation, CPM, transformation, and quantile

## NormalisationStep -- plot the overlap rate 

## CountPerMilions -- Transform the count matrix in CPM values

## transformation -- perform various data transformations

## QuantileNormalization -- performs quantile normalisation of a matrix

############################################################################################################################################
##
## normalization -- wrapper function for normalisation, CPM, transformation, and quantile
##
## DB_object   -- DB_object   : a differential binding object created earlier
## normMethod  -- character   : action to do with normalisation ("TMM","UQ","DESeq","Median","RLE","RUVg","RUVs","RUVr","none")
## CPM         -- character   : action to do with CPM transformation ("CPMOutlierRemoved","CPM","none")
## transform   -- character   : action to do with transformations ("logTransform","voomTransform","vstTransform","clrTransform","none")
## quantile    -- character   : action to do with quantile normalisation ("quantile","none")
## 
##
## returns a DB object containing the info in the csv files, the peaks for each sample, the consensus and the coverage, binding affinity matrix and normalized matrix
############################################################################################################################################
normalization <- function(DB_object, normMethod="none" , CPM="none", transform="none", quantile="none"){
  
  if(normMethod %!in% c("TMM","UQ","DESeq","Median","RLE","RUVg","RUVs","RUVr","none")){
    stop(" The type of normalisation method selected is not implemented, please select TMM,UQ,DESeq,Median,RLE,RUVg,RUVs,RUVr,none")
  }
  if(transform %!in% c("logTransform","voomTransform","vstTransform","clrTransform","none")){
    stop(" The type of transform method selected is not implemented, please select logTransform,voomTransform,vstTransform,clrTransform,none")
  }
  if(CPM %!in% c("CPMOutlierRemoved","CPM","none")){
    stop(" The type of normalisation method selected is not implemented, please select CPMOutlierRemoved,CPM,none")
  }
  if(quantile %!in% c("quantile","none")){
    stop(" The type of normalisation method selected is not implemented, please select TMM,UQ,DESeq,Median,RLE,RUVg,RUVs,RUVr,none")
  }
  

  matrix<- NormalizationStep(DB_object = DB_object, 
                             firstGroup = table(DB_object$Info$V3)[[1]], 
                             secondGroup = table(DB_object$Info$V3)[[2]], 
                             normMethod = normMethod)
  
  matrix <- CountPerMilion(matrix,
                           CPM, 
                           lib.size = colSums(DB_object$Counts))
  
  matrix <- transformation(matrix,
                    transform,
                    firstGroup = table(DB_object$Info$V3)[[1]],
                    secondGroup = table(DB_object$Info$V3)[[2]])
  
  matrix <- QuantileNormalization(matrix,
                                  quantile)
  
  return_object <- list(DB_object$Peaks,DB_object$Consensus,DB_object$Coverage,DB_object$Info,DB_object$Counts,matrix)
  names(return_object) <- c("Peaks","Consensus","Coverage","Info","Counts","CountsNorm")
  
  return(return_object)
  
}
  

############################################################################################################################################
##
## NormalizationStep -- Perform the selected normalisation
##
## DB_object    -- DB_object   : a differential binding object created earlier
## firstGroup   -- vector      : vector of numbers of the size of the group 1
## secondGroup  -- character   : vector of numbers of the size of the group 2
## normMethod   -- character   : Normalisation method choosen ("TMM","UQ","DESeq","Median","RLE","RUVg","RUVs","RUVr","none")
##
############################################################################################################################################

#Stepwise functions
{
  #Step 1
  {
    NormalizationStep <- function(DB_object, firstGroup, secondGroup, normMethod) {
      
      if("CountsNorm" %in% names(DB_object)){
        rawCountsDF <- DB_object$CountsNorm
      }
      else{
        rawCountsDF <- DB_object$Counts
      }
      
      #rawCountsDF <- DB_object$Counts
      genes <- DB_object$Consensus$GeneID
      Treat <- c(rep("GRP1", firstGroup), rep("GRP2", secondGroup))
      switch(normMethod,
             TMM = {
               #Normalization factors for TMM
               {
                 y <- edgeR::DGEList(counts=rawCountsDF, group=Treat)
                 y <- edgeR::calcNormFactors(object = y, method = "TMM")
               }
               #Applying the normalization factors
               {
                 TMM.Normalized <- t(t(rawCountsDF)/y$samples$norm.factors)
                 rownames(TMM.Normalized) <- genes
               }
               return(TMM.Normalized)
             },
             UQ = {
               y <- DGEList(counts=rawCountsDF, group=Treat)
               y <- calcNormFactors(object = y, method = "upperquartile")
               Normalized <- t(t(rawCountsDF)/y$samples$norm.factors)
               rownames(Normalized) <- genes
               return(Normalized)
             },
             DESeq = {
               condition <- Treat
               type <- Treat
               Treat <-cbind(condition, type)
               rownames(Treat) <- c(1:(firstGroup+secondGroup))
               colnames(Treat) <- c("condition","type")
               #Treat[,"condition"] <- factor(Treat[,"condition"])
               #Treat[,"type"] <- factor(Treat[,"type"])
               
               
               cds <- DESeqDataSetFromMatrix(countData = rawCountsDF,
                              colData = Treat,design = ~ condition)
               #cds = newCountDataSet(rawCountsDF, conditions=Treat)
               cds = estimateSizeFactors(cds)
               DESeq.Normalized = as.data.frame(round(counts(cds, normalized=TRUE)))
               rownames(DESeq.Normalized) <- genes
               return(DESeq.Normalized)
             },
             Median = {
               Median.Normalized <- as.data.frame(round(limma::normalizeMedianAbsValues(rawCountsDF)))
               rownames(Median.Normalized) <- genes
               return(Median.Normalized)
             },
             RLE = {
               {
                 y <- DGEList(rawCountsDF, group=Treat)
                 y <- calcNormFactors(object = y, method = "RLE")
                 #return(edgeR::cpm(y))
               }
               #Applying the normalization factors
               {
                 Normalized <- t(t(rawCountsDF)/y$samples$norm.factors)
                 rownames(Normalized) <- genes
                 return(Normalized)
               }
             },
             RUVg = {
                 x <- Treat
                 set <- newSeqExpressionSet(rawCountsDF,
                                            phenoData = data.frame(x, row.names=colnames(rawCountsDF)))
                 
                 design <- model.matrix(~x, data=pData(set))
                 y <- DGEList(counts=counts(set), group=x)
                 y <- calcNormFactors(y, method="upperquartile")
                 y <- estimateGLMCommonDisp(y, design)
                 y <- estimateGLMTagwiseDisp(y, design)
                 fit <- glmFit(y, design)
                 lrt <- glmLRT(fit, coef=2)
                 genes <- as.character(genes)
                 #RUVg
                 {
                   top <- topTags(lrt, n=nrow(set))$table
                   empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:(floor(nrow(rawCountsDF) * (15/100)))]))]
                   set <- RUVg(set, empirical, k=1)
                   Normalized <- EDASeq::normCounts(set)
                   rownames(Normalized) <- genes
                   return(Normalized)
                 }
             },
             RUVs={
               x <- Treat
               set <- newSeqExpressionSet(rawCountsDF,
                                          phenoData = data.frame(x, row.names=colnames(rawCountsDF)))
               
               design <- model.matrix(~x, data=pData(set))
               y <- DGEList(counts=counts(set), group=x)
               y <- calcNormFactors(y, method="upperquartile")
               y <- estimateGLMCommonDisp(y, design)
               y <- estimateGLMTagwiseDisp(y, design)
               fit <- glmFit(y, design)
               lrt <- glmLRT(fit, coef=2)
               genes <- as.character(genes)
               
               differences <- makeGroups(x)
               differences
               
               set <- RUVs(set, genes, k=1, differences)
               Normalized <- normCounts(set)
               rownames(Normalized) <- genes
               return(Normalized)
               
             },
             RUVr ={
               x <- Treat
               set <- newSeqExpressionSet(rawCountsDF,
                                          phenoData = data.frame(x, row.names=colnames(rawCountsDF)))
               
               design <- model.matrix(~x, data=pData(set))
               y <- DGEList(counts=counts(set), group=x)
               y <- calcNormFactors(y, method="upperquartile")
               y <- estimateGLMCommonDisp(y, design)
               y <- estimateGLMTagwiseDisp(y, design)
               fit <- glmFit(y, design)
               lrt <- glmLRT(fit, coef=2)
               genes <- as.character(genes)
               
               differences <- makeGroups(x)
               differences
               
               res <- residuals(fit, type="deviance")
               set <- RUVr(x = set, cIdx = genes, k=1, res)
               
               Normalized <- normCounts(set)
               rownames(Normalized) <- genes
               return(Normalized)
               },
             none = {
               rownames(rawCountsDF) <- genes
               return(rawCountsDF)
             }
      )
    }
    
    #For testing
    # Normalized <- NormalizationStep(rawCountsDF = simulatedData, 
    #                                 firstGroup = 5, 
    #                                 secondGroup = 5, 
    #                                 genes = rownames(simulatedData), 
    #                                 normMethod = "RLE")
  }

  ############################################################################################################################################
  ##
  ## CountPerMilion -- Perform the selected CPM transformation action
  ##
  ## Normalized      -- matrix      : matrix to perform the action on
  ## outlierAction   -- character   : action choosen ("CPMOutlierRemoved","CPM","none")
  ## lib.size        -- vector      : vector of integer of the size of the number of samples
  ##
  ############################################################################################################################################  

  {
    CountPerMilion <- function(Normalized, outlierAction, lib.size) {
      lib.size <- 1e-6*lib.size
      switch(outlierAction,
             CPMOutlierRemoved = {
               ###Remove outliers here
               flat <- as.numeric(unlist(Normalized))
               NinetiethPercentile <- as.numeric(quantile(flat, c(.9)))
               keep <- as.logical(unlist(apply(X = Normalized, MARGIN = 1, FUN = function(r) any(r < NinetiethPercentile)), use.names = F))
               return(t(t(Normalized[keep,])/lib.size))
             },
             CPM = {
               return(t(t(Normalized)/lib.size))
             }, 
             none = {
               return(Normalized)
             }
      )
    }
    #For testing 
    # test <- CountPerMilion(Normalized = Normalized, outlierAction = "yes")
  }

  
  
  ############################################################################################################################################
  ##
  ## transformation  -- Perform the selected data transformation action
  ##
  ## input        -- matrix      : matrix to perform the action on
  ## type         -- character   : action choosen ("logTransform","voomTransform","vstTransform","clrTransform","none")
  ## firstGroup   -- vector      : vector of numbers of the size of the group 1
  ## secondGroup  -- character   : vector of numbers of the size of the group 2
  ##
  ############################################################################################################################################  
  {
    transformation <- function(input, type, firstGroup, secondGroup) {
      Treat <- factor(c(rep("GRP1", firstGroup), rep("GRP2", secondGroup)))
      switch(type,
             logTransform = {
               return(log2(floor(input) + 1))
             },
             voomTransform = {
               return(limma::voom(counts = floor(input), design = as.matrix(as.numeric(Treat)))$E)
             },
             vstTransform = {
               return(DESeq2::vst(object = as.matrix((round(input, digits = 0)))))
             },
             clrTransform = {
               return(ALDEx2::aldex.clr(reads = floor(as.data.frame(input)), conds = as.character(Treat))@reads)
             },
             none = {
               return(input)
             }
      )
    }
    #For testing
    # test <- transformation(input = test, type = "voomTransform", firstGroup = length(samps$GRP1), secondGroup = length(samps$GRP5))
  }
  
  
  ############################################################################################################################################
  ##
  ## QuantileNormalization  -- Perform the selected quantile normalisation action
  ##
  ## input        -- matrix      : matrix to perform the action on
  ## action       -- character   : action choosen ("quantile","none")
  ##
  ############################################################################################################################################  
  {
    QuantileNormalization <- function(input, action) {
      switch(action,
             quantile = {
               return(limma::normalizeQuantiles(input))
             },
             none = {
               return(input)
             }
      )
    }
    # test <- QuantlieNormalization(input = test, action = "quantile")
  }
  
}


#Loading libraries
{
  #The libraries that needs to be used are loaded here
  #SummarizedExperiment
  {
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("SummarizedExperiment")
    library("SummarizedExperiment")
    
  }
  
  #survival
  {
    #install.packages("survival")
    library("survival")
  }
  
  #easyRNASeq
  {
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("easyRNASeq")
    #library("easyRNASeq")
  }
  #limma
  # Data analysis, linear models and differential expression for microarray data.
  {
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("limma")
    #library(limma)
  }
  #edgeR
  # edgeRUsersGuide()
  {
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("edgeR")
    #library("edgeR")
  }
  #preprocessCore
  {
    #source('http://bioconductor.org/biocLite.R')
    #biocLite('preprocessCore')
    #library(preprocessCore)
  }
  #DESeq
  # Estimate variance-mean dependence in count data 
  # from high-throughput sequencing assays and test 
  # for differential expression based on a model using 
  # the negative binomial distribution
  # browseVignettes("DESeq")
  {
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("DESeq")
    #library(DESeq)
  }
  
  {
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("DESeq2")
    library(DESeq2)
  }
  # browseVignettes("seqc")
  {
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("seqc")
    #library(seqc)
  }
  
  {
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("affy")
    #library(affy)
  }
  #Another dataset
  {
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("NBPSeq")
    #library(NBPSeq)
  }
  {
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("ROTS")
    #library(ROTS)
  }
  {
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("RUVSeq")
    #library("RUVSeq")
  }
  {
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("ALDEx2")
    #library("ALDEx2")
  }
  #library(parallel)
}

