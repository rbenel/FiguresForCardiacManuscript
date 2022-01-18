library(DESeq2)
library(vsn)
library(dplyr)
library(ggplot2)
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library(limma)
library(hexbin)

###############
#LOAD THE DATA
###############
#load the data from the previous script, we want to utilize the deseq object here 

load("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/ddsTxi_obj.RData"))

#Since common statistical methods for exploratory analysis of multidimensional data work best for data that 
#generally has the same range of variance at different ranges of the mean values. This kind of data is said to be homoskedastic. 
 
#However RNA-seq counts the variance grows with the mean.To correct for this, the DESeq package uses two types of transformations
#The VST or rlog-transformed data then become approximately homoskedastic (more flat trend in the 
#meanSdPlot), and can be used directly for computing distances between samples, making PCA plots, 
#or as input to downstream methods which perform best with homoskedastic data.

###################
#TRANSORM THE DATA
###################

#as this data contains < than 30 samples we use rlog as reccomended in tutorial as the transformation 
rld <- rlog(ddsTxi, blind = FALSE)
head(assay(rld), 3)

meanSdPlot(assay(rld))

mat <- assay(vsd)
mat2 <- limma::removeBatchEffect(mat, vsd$condition)
#save(mat2, file = paste0(local_path, "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/vsd.norm.counts.RData"))

###################
#HEATMAP OF SAMPLES
###################
#in order to see how similar the samples are to one another we use euclidian distance 
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists ) #turn into a matrix 
rownames(sampleDistMatrix) <- paste( rld$day, rld$replicate, sep = " - " ) #give useful rownames 
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

#Heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists, #use the sample distance because default is that the matrix contains the data values which it doesnt...
         clustering_distance_cols = sampleDists,
         col = colors)

# #if we want to see how dissimilar the samples are
# # This measure of dissimilarity between counts also takes the inherent variance structure of 
# # counts into consideration when calculating the distances between samples. The PoissonDistance 
# # function takes the original count matrix (not normalized) with samples as rows instead of columns, 
# # so we need to transpose the counts in dds.
# 
# poisd <- PoissonDistance(t(counts(ddsTxi)))
# samplePoisDistMatrix <- as.matrix( poisd$dd )
# rownames(samplePoisDistMatrix) <- paste( ddsTxi$day, ddsTxi$replicate, sep=" - " )
# colnames(samplePoisDistMatrix) <- NULL
# 
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# 
# #heatmap
# pheatmap(samplePoisDistMatrix,
#          clustering_distance_rows = poisd$dd,
#          clustering_distance_cols = poisd$dd,
#          col = colors)

##########
#PLOT PCA
##########
pcaData <- plotPCA(vsd, intgroup=c("day", "replicate"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
pca_rlog <- ggplot(pcaData, aes(PC1, PC2, color=day, shape=replicate)) +
            geom_point(size=3) +
            xlab(paste0("PC1: ",percentVar[1],"% variance")) +
            ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
            coord_fixed() + 
            scale_color_manual(values=c("#f9bdbd", "#f74a4a", "#ffaa56", "#f9e67a", "#a6dd6a", "#05aa60", 
                                                  "#b0e5df",  "#4bb7e5","#648cea", "#caa9e5","#9a25d1" )) + #rainbow colors 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
            theme(panel.background = element_blank()) +
            theme(axis.line.x = element_line(color="black", size = 0.5),
                  axis.line.y = element_line(color="black", size = 0.5)) 





