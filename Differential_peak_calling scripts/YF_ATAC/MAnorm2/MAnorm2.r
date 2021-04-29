library(MAnorm2)


## Apply MAnorm2 to raw read counts.

# Read data.
profile_bins <- read.table("profile_bins.xls",
                           header = TRUE, stringsAsFactors = FALSE)
head(profile_bins)
n1 <- 16
n <- 16

# Perform within-group normalization.
norm <- MAnorm2::normalize(profile_bins, count = 4:19, occupancy = 20:35)


# Construct a bioCond for each group of samples.
conds <- list(group1 = bioCond(norm[4:11], norm[20:27], name = "naive"),
              group2 = bioCond(norm[12:19], norm[28:35], name = "vaccinated"))

# Perform between-group normalization.
conds <- normBioCond(conds)

# Fit a mean-variance curve.
conds <- fitMeanVarCurve(conds, method = "parametric", init.coef = c(0.1, 10))

# Perform statistical tests for identifying differential signals.
res <- diffTest(conds[[1]], conds[[2]])
head(res)

#calculate fold change
FC <- log2((rowMeans(profile_bins[,4:11])+0.1)/(rowMeans(profile_bins[,12:19])+0.1))
results <- cbind(profile_bins[,1:3],res,FC)

names <- seq(dim(results)[1])
names <- paste0("peak_",names)
strand <- rep("*",dim(results)[1])
formatedresults <- cbind(results[,1:3],names,results[,3]-results[,2],strand,results$pval,results$padj,results$FC)

tab <- formatedresults[order(results$padj),]

write.table(x = tab, 
            file = "MAnorm2_pooledpeaks.bed", 
            quote = F , 
            col.names = F , 
            row.names = F , 
            sep = "\t" )




