library(GOplot)
library(viridis)
library(dplyr)
library(tidyr)

#this function comes from this package 
#https://wencke.github.io/


#load annotation DF from heatmap script 
load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/transcript_results/March20_OnlyDTUtranscriptsAnnotN16WithGeneName.df.RData"))

finalAnnot.df <- annot.df

#load res def padj values from LRT in order to add FC to the final graph
#the FC here would be Wald for Day0 day 60 becuase LRT doesnt give FC
load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/updatedPadjTSL_LRT_lincResDF.RData"))

#updated lnc, linc, and funcLncRNA annotations
load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/Practice 18.12.16/data/Tx.linc_lncRNA_FunctlncRNA_release92TSL.RData"))

#our biological question:
#find which transcripts move from one cluster to another 

diffClusGenes <- data.frame()
#so we only need to focus on those with unique elements in the cluster list column
for(i in unique(finalAnnot.df$geneID)){
  subAnotGene <- subset(finalAnnot.df, finalAnnot.df$geneID == i)
  cluster.list <- data.frame(subAnotGene$cluster)
  if(length(unique(cluster.list[,1])) == 1){  #checks if unique elements of the column
    newRow <- data.frame(gene_name = i, cluster = cluster.list[1,1]) #if this is true we don't need to focus on this take first element bec they are the same 
    next #if we don't want to look at the ones that stay, for the graph!!!
  } else {
    clean.cluster <- cluster.list %>% distinct() #save only unique entry of cluster name 
    clean.cluster <- as.data.frame(t(clean.cluster)) #transform so we can collapse a row 
    clean.cluster$V3 <- apply(clean.cluster, 1, paste, collapse="_") #collapse the row to one word
    newRowElse <- data.frame(gene_name = i, cluster = clean.cluster$V3) #paste(i, "F", sep = ",") #this would allow us to see T and F 
    diffClusGenes <- rbind(diffClusGenes, newRowElse)
  }
}


#You can build the matrix on your own or you use the implemented function chord_dat which does the job for you
NumClusters <- c(1,2,3)
for(i in NumClusters){
  
  ClusterCheck <- as.data.frame(ifelse(grepl(paste0("cluster_", i, ".*"), diffClusGenes$cluster), 1, 0))
  colnames(ClusterCheck) <- paste0("cluster_", i)
  diffClusGenes <- cbind(diffClusGenes, ClusterCheck)
  
}


######################################################
#Add LRT FC which is really from Wald test of Day0 Day60
######################################################
#Add the FC to the plot
#this is the data from normFC_analysis
padj_linc <- padj_linc %>% dplyr::select(log2FoldChange) 

#change the column names so that the function knows it is receiving FC 
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

#remove the initial column called "cluster" that had the combination of the clusters, for the graph
GoChordPlotDF <- GoChordPlotDF[ , !colnames(GoChordPlotDF) %in% "cluster"] 

###############
#Actual Plot
###############
#chose ribbon colors that match our heatmap 
ribbon_colors <- c(viridis::magma(20))

#ordered colors
OrderedRibbonColors <- ribbon_colors[c(12, 16, 6)] #for three clusters

#file to save output
pdf(file = file.path("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/transcript_results/TxPlots_lincRNA/otherPlots/N16_3ClusterMarch2021GoChordPlotWithFCNoCluster4.pdf"), width = 17, height = 17)

GOChord(GoChordPlotDF, space = 0.02, gene.space = 0.25, gene.size = 5.3,
        ribbon.col	= OrderedRibbonColors,  gene.order = 'logFC', 
        lfc.min = -4.5 , lfc.max = 5, 
        limit = c(2,3))

dev.off()
