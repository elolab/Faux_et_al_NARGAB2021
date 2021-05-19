library(GenomicRanges)
library(RepViz)
library(ChIPpeakAnno)

#This script is using the bioconductor package RepViz (Faux et al. 2019) to plot the read densities 

###################################################
##
##           Figure 4
##
###################################################



region<- GRanges("chr1:2082000-2084000")
jpeg("Figure_4A_chr1:2082000-2083000.jpeg", width = 1400, height = 1000)
RepViz(region = region,
           genome = "hg19",
           BAM = "BAM_input.csv",
           BED = "BED_input.csv",
           avgTrack = T,
           geneTrack = T,
           verbose=T,
           col = c("orange","red","green","violet","blue","black","black"),cex = 2)

dev.off()

region<- GRanges("chr1:11789000-11793000")
jpeg("Figure_4B_chr1:11789000-11793000.jpeg", width = 1400, height = 1000)
RepViz(region = region,
       genome = "hg19",
       BAM = "BAM_input.csv",
       BED = "BED_input.csv",
       avgTrack = T,
       geneTrack = T,
       verbose=T,
       col = c("orange","red","green","violet","blue","black","black"),cex = 2)

dev.off()



region<- GRanges("chr11:61594000-61599000")
jpeg("Fig_4C_chr11:61594000-61599000.jpeg", width = 1400, height = 1000)
RepViz(region = region,
           genome = "hg19",
           BAM = "BAM_input.csv",
           BED = "BED_input.csv",
           avgTrack = T,
           geneTrack = T,
           verbose=T,
       col = c("orange","red","green","violet","blue","black","black"),cex = 2)
dev.off()


###################################################
##
##           Supplementary figure 4
##
###################################################

# genomic region as a GRanges object
region<- GRanges("chr14:25102000-25105000")

#RepViz call
jpeg("Supplementary_3A_chr14:25102000-25105000.jpeg", width = 1400, height = 1000)
RepViz(region = region,
           genome = "hg19",
           BAM = "BAM_input.csv",
           BED = "BED_input.csv",
           avgTrack = T,
           geneTrack = T,
           verbose=T,
       col = c("orange","red","green","violet","blue","black","black"),cex = 2)
dev.off()


region<- GRanges("chr1:25390500-25393650")
jpeg("supplementary_3B_chr1:25390500-25393650.jpeg", width = 1400, height = 1000)
RepViz(region = region,
       genome = "hg19",
       BAM = "BAM_input.csv",
       BED = "BED_input.csv",
       avgTrack = T,
       geneTrack = T,
       verbose=T,
       col = c("orange","red","green","violet","blue","black","black"),cex = 2)
dev.off()