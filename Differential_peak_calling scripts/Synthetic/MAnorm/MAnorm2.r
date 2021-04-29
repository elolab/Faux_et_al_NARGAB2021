## Using the differential analysis of H3K4me3 ChIP-seq data between GM12891 and GM12892
## LCLs (lymphoblastoid cell lines) to illustrate the workflow of MAnorm2.
library(MAnorm2)


## Apply MAnorm2 to raw read counts.

# Read data.
profile_bins <- read.table("Simulated_profile_bins.xls",
                           header = TRUE, stringsAsFactors = FALSE)

# trying to get the data normally distributed
profile_bins <- normalizeBySizeFactors(profile_bins,count = c("s11.read_cnt",
                                                                   "s21.read_cnt",
                                                                   "s31.read_cnt",
                                                                   "s41.read_cnt",
                                                                   "s51.read_cnt",
                                                                   "s12.read_cnt",
                                                                   "s22.read_cnt",
                                                                   "s32.read_cnt",
                                                                   "s42.read_cnt",
                                                                   "s52.read_cnt"))

# Perform within-group normalization.
#norm <- normalize(profile_bins, count = 4:19, occupancy = 20:35)


# Construct a bioCond for each group of samples.
group1 <- bioCond(profile_bins[4:8], profile_bins[14:17], name = "control")
group2 <- bioCond(profile_bins[9:13], profile_bins[18:22], name = "ChIP")
conds <- list(group1 = bioCond(profile_bins[4:8], profile_bins[14:17], name = "control"),
              group2 = bioCond(profile_bins[9:13], profile_bins[18:22], name = "ChIP"))

# Perform between-group normalization.
#conds <- normBioCond(conds)

# Fit a mean-variance curve.
conds <- fitMeanVarCurve(conds, method = "parametric", init.coef = c(0.1, 10))

# Perform statistical tests for identifying differential signals.
res <- diffTest( conds[[2]],conds[[1]])
MAplot.bioCond(group1,group2)
MAplot.diffBioCond(res)

head(res)

FC <- log2((rowMeans(profile_bins[,4:8])+0.1)/(rowMeans(profile_bins[,9:13])+0.1))
results <- cbind(profile_bins[,1:3],res,FC)

names <- seq(dim(results)[1])
names <- paste0("peak_",names)
strand <- rep("*",dim(results)[1])
formatedresults <- cbind(results[,1:3],names,results[,3]-results[,2],strand,results$pval,results$padj,results$FC)

tab <- formatedresults[order(results$padj),]

write.table(x = tab, 
            file = "MAnorm2_medianOfRatios.bed", 
            quote = F , 
            col.names = F , 
            row.names = F , 
            sep = "\t" )
## Apply MAnorm2 to normalized log2 read counts.

# Read data.
#normalized <- read.table("data/GM12891_cmp_GM12892.H3K4me3_hieMA.xls",
#                         header = TRUE, stringsAsFactors = FALSE)
#head(normalized)
#norm2 <- profile_bins 
#norm2[4:(3 + n)] <- normalized

# Follow the original workflow for differential analysis.
#conds2 <- list(GM12891 = bioCond(norm2[4:(3 + n1)], norm2[4:(3 + n1) + n], name = "GM12891"),
#               GM12892 = bioCond(norm2[(4 + n1):(3 + n)], norm2[(4 + n1):(3 + n) + n], name = "GM12892"))
#conds2 <- fitMeanVarCurve(conds2, method = "parametric", init.coef = c(0.1, 10))
#res2 <- diffTest(conds2[[1]], conds2[[2]])
#head(res2)










####test


MAplot <- function(my_matrix, name){
  
  groups <- c(1,1,1,1,1,2,2,2,2,2)
  
  matrix <- my_matrix+1

  png(height = 900,width = 900,filename =  name)
  plot( rowMeans(log2(my_matrix)), log2(apply(my_matrix[,which(groups==1)],1,mean))-log2(apply(my_matrix[,which(groups==2)],1,mean)), cex=1 )
  dev.off()
  
}
