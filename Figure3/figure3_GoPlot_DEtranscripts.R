#install.packages('GOplot')
library(GOplot)
library(viridis)
library(dplyr)
library(tidyr)

#this functions comes from this package 
#https://wencke.github.io/

#local_path <- getwd()

#load data
#load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/transcript_results/March20_OnlyDTUtranscriptsAnnot.dfLoweredFilterN16.RData"))

load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/transcript_results/March20_OnlyDTUtranscriptsAnnotN16WithGeneName.df.RData"))

finalAnnot.df <- annot.df

#March 24th 2020
#load res def padj values from LRT to add FC or padj to the final graph
#here if we use the FC this is for Day 60. 
load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/updatedPadjTSL_LRT_lincResDF.RData"))

#March 24th 2020
#updated lnc, linc, and funcLncRNA annotations
load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/Practice 18.12.16/data/Tx.linc_lncRNA_FunctlncRNA_release92TSL.RData"))

#this is different from gencode_DTU heatmap, here we **only** have transcripts that cross borders!

diffClusGenes <- data.frame()
#now we want to know which lincRNA genes have transcripts in different clusters... 
#so we only need to focus on those with unique elements in the cluster list column
for(i in unique(finalAnnot.df$geneID)){
  subAnotGene <- subset(finalAnnot.df, finalAnnot.df$geneID == i)
  cluster.list <- data.frame(subAnotGene$cluster)
  if(length(unique(cluster.list[,1])) == 1){  #checks if unique elements of the column
    newRow <- data.frame(gene_name = i, cluster = cluster.list[1,1]) #if this is true we don't need to focus on this take first element bec they are the same 
    next #if we don't want to look at the ones that stay, for the graph!!!
    #diffClusGenes <- rbind(diffClusGenes, newRow)
  } else {
    #this is before we started to look at everything 
    #newRowElse <- c(i) #paste(i, "F", sep = ",") #this would allow us to see T and F 
    #diffClusGenes <- append(diffClusGenes, newRowElse
    clean.cluster <- cluster.list %>% distinct() #save only unique entry of cluster name 
    clean.cluster <- as.data.frame(t(clean.cluster)) #transform bec we can collapse a row 
    clean.cluster$V3 <- apply(clean.cluster, 1, paste, collapse="_") #collapse the row to one word
    newRowElse <- data.frame(gene_name = i, cluster = clean.cluster$V3) #paste(i, "F", sep = ",") #this would allow us to see T and F 
    diffClusGenes <- rbind(diffClusGenes, newRowElse)
  }
}

#this is a better way to do it, what is not commented, so the collapse "+" went back to "_"

#string seperate all of the clusters so that we have each cluster in a seperate columns
#diffClusGenes$firstCluster <- sapply(strsplit(as.character(diffClusGenes$cluster), "+", fixed = TRUE), '[', 1)
#diffClusGenes$secondCluster <- sapply(strsplit(as.character(diffClusGenes$cluster), "+", fixed = TRUE), '[', 2)
#diffClusGenes$thirdCluster <- sapply(strsplit(as.character(diffClusGenes$cluster), "+", fixed = TRUE), '[', 3)
#diffClusGenes$fourthCluster <- sapply(strsplit(as.character(diffClusGenes$cluster), "+", fixed = TRUE), '[', 4)

#You can build the matrix on your own or you use the implemented function chord_dat which does the job for you
#correct for previous plot march 2020
#NumClusters <- c(3, 4, 5, 2, 1, 7) #this is the order of the colors on the heatmap?

#aug 2021 
NumClusters <- c(1,2,3)
for(i in NumClusters){
  
  ClusterCheck <- as.data.frame(ifelse(grepl(paste0("cluster_", i, ".*"), diffClusGenes$cluster), 1, 0))
  colnames(ClusterCheck) <- paste0("cluster_", i)
  diffClusGenes <- cbind(diffClusGenes, ClusterCheck)
  
}

# diffClusGenes$cluster1 <- ifelse(grepl("cluster_1.*", diffClusGenes$cluster), 1, 0)
# diffClusGenes$cluster4 <- ifelse(grepl("cluster_4.*", diffClusGenes$cluster), 1, 0)
# diffClusGenes$cluster3 <- ifelse(grepl("cluster_3.*", diffClusGenes$cluster), 1, 0)
# diffClusGenes$cluster6 <- ifelse(grepl("cluster_6.*", diffClusGenes$cluster), 1, 0)
# diffClusGenes$cluster2 <- ifelse(grepl("cluster_2.*", diffClusGenes$cluster), 1, 0)
# diffClusGenes$cluster5 <- ifelse(grepl("cluster_5.*", diffClusGenes$cluster), 1, 0)
# diffClusGenes$cluster7 <- ifelse(grepl("cluster_7.*", diffClusGenes$cluster), 1, 0)

######################################################
#FROM LRT DATA - padj or FC - but FC is from first day
######################################################
###Add the padj or FC to the plot
#this is the data from normFC_analysis
padj_linc <- padj_linc %>% dplyr::select(log2FoldChange) 

#padj_linc <- -log10(padj_linc) #this doesn't really help the visualization 
#change the column names so that the function thinks it is receiving FC 
colnames(padj_linc) <- "logFC"

#join the LRT data with the diffClustGenes
joinedGoChordPlot <- merge(diffClusGenes, padj_linc, by.x = "gene_name", by.y = "row.names", all.x = T)

#turn into rownames, also for the function. 
rownames(joinedGoChordPlot) <- joinedGoChordPlot$gene_name
GoChordPlotDF <- joinedGoChordPlot[, -1]

#if padj, then the NA needs to be 1, if FC needs to be 0
GoChordPlotDF$logFC[is.na(GoChordPlotDF$logFC)] <- 0

#match their external gene name for ensembl gene name
rownames(GoChordPlotDF) <- Tx.lincRNA[ match(rownames(GoChordPlotDF), Tx.lincRNA[['gene_id']] ) , 'gene_name']
###################
#remove the initial column called "cluster" that had the combination of the clusters, for the graph
GoChordPlotDF <- GoChordPlotDF[ , !colnames(GoChordPlotDF) %in% "cluster"] 

#chose ribbon colors that match our heatmap, add the green as cluster 7, the "undiffrentiatied"
ribbon_colors <- c(viridis::magma(20))
#OrderedRibbonColors <- ribbon_colors[c(4, 16, 8, 6, 12, 18, 21)]
#march 2020 colors
#OrderedRibbonColors <- ribbon_colors[c(6, 8, 18, 16, 12, 21)]
#aug 2021 colors
OrderedRibbonColors <- ribbon_colors[c(12, 16, 6)] #for three clusters
#OrderedRibbonColors <- c(rep("grey", 3)) #saved a grey version :)

pdf(file = file.path("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/transcript_results/TxPlots_lincRNA/otherPlots/N16_3ClusterMarch2021GoChordPlotWithFCNoCluster4.pdf"), width = 17, height = 17)

#plot
GOChord(GoChordPlotDF, space = 0.02, gene.space = 0.25, gene.size = 5.3,
        ribbon.col	= OrderedRibbonColors,  gene.order = 'logFC', 
        lfc.min = -4.5 , lfc.max = 5, 
        limit = c(2,3))

dev.off()
