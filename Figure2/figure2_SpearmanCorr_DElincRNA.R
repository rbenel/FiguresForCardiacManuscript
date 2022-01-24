library(dplyr)
library(tidyr)
library(pheatmap)
library(reshape)
library(data.table)
library(ggplot2)
library(tidyverse)
library(corrr)
library(igraph)
library(ggraph)
library(ggrepel)
library(viridis)

#all DE lincRNA
load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/updatedPadjTSLnormFiltLRT_linc.RData"))


#Updated March 2021
load(file = "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/Practice 18.12.16/data/NoCRISPRupdatedFunc_Tx.linc_lncRNA_FunctlncRNA_release92TSL.RData")

##################################
#Spearman Correlation all lincRNA
##################################
#spearman correlation 
correlation_spearman <-cor(t(norm_FC_lincRNA), t(norm_FC_lincRNA),  method = c("spearman"))

#This is the same thing
# cor(t(norm_FC_lincRNA), method = "spearman")

correlation_spearman = as.matrix(correlation_spearman)

color <- colorRampPalette(c("#579af2", "white", "#cc1a35"))

hist <- hist(correlation_spearman, breaks=20, xlab= "Correlation Values", main="Distribution of Spearman Correlation",
             col= color(20))

pheatmap(correlation_spearman ,
         # breaks = seq(from=-thr, to=thr, length=101),
         scale = "none",
         #kmeans_k = 4,
         border_color = NA,
         clustering_method = "complete",
         cluster_cols = T,
         cluster_rows = T,
         show_rownames = F,
         show_colnames = F,
         cellwidth = NA,
         cellheight = NA,
         revC = F,
         treeheight_row = 0,
         treeheight_col = 0,
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdBu")))(256),
         main = "Correlation of DE lincRNAs"
)

#save(correlation_spearman, file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/April20_lincRNA_lincRNA_CorrScores.RData"))

#Added on May 5 2021
############################
#Clustering and annotations
############################
#Make some visual changes to the heatmap 
numClusters <- 5
lincRNAheatmap <- pheatmap(correlation_spearman ,
                   # breaks = seq(from=-thr, to=thr, length=101),
                   scale = "none",
                   #kmeans_k = 4,
                   cutree_rows = numClusters,
                   border_color = NA,
                   clustering_method = "complete",
                   cluster_cols = T,
                   cluster_rows = T,
                   show_rownames = F,
                   show_colnames = F,
                   cellwidth = NA,
                   cellheight = NA,
                   revC = F,
                   treeheight_row = 20,
                   treeheight_col = 0,
                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdBu")))(256),
                   main = "Correlation of DE lincRNAs"
          )

####################
#Heatmap Annotation 
###################
#cut the heatmap tree
heatmapCutree <- cutree(lincRNAheatmap$tree_row, k = numClusters)

#on May 30th 2021 - turn into letters so less confusing in text
heatmapCutree <- sapply(heatmapCutree, function(i) toupper(letters[i]))


#choose colors for the heatmap 
heatmapColors <- RColorBrewer::brewer.pal(5, 'Set2')

#annotated df of rownames + name of cluster 
annotHeatmap <- data.frame(cluster = paste("cluster", heatmapCutree, sep = "_"))
rownames(annotHeatmap) <- rownames(correlation_spearman)

#create a df with functional and lincRNA
FunctionalLincRNA <- subset(Tx.lncRNAfunctional, Tx.lncRNAfunctional$tx_biotype == "lincRNA")

#see how many are in each cluster
table(annotHeatmap)

#how many are functional
length(intersect(rownames(annotHeatmap), Tx.lncRNAfunctional$gene_id))

#add functional to be part of the annnotation DF
annotHeatmap$functional <- as.factor(ifelse(rownames(annotHeatmap) %in% Tx.lncRNAfunctional$gene_id, "func", ""))

#annotHeatmap$functional <- relevel(annotHeatmap$functional, ref = 0)
#assing colors to functional 
functionalColors = c( "white", "black")

table(annotHeatmap)

#add colors to the genes 
names(heatmapColors) <- unique(annotHeatmap$cluster)
names(functionalColors) <- unique(annotHeatmap$functional)
my.colors <- list(cluster = heatmapColors,
                  functional = functionalColors)



lincRNAheatmapColored <- pheatmap(correlation_spearman ,
                           # breaks = seq(from=-thr, to=thr, length=101),
                           scale = "none",
                           #kmeans_k = 4,
                           cutree_rows = numClusters,
                           border_color = NA,
                           clustering_method = "complete",
                           cluster_cols = T,
                           cluster_rows = T,
                           show_rownames = F,
                           show_colnames = F,
                           cellwidth = NA,
                           cellheight = NA,
                           revC = F,
                           treeheight_row = 20,
                           treeheight_col = 0,
                           annotation_row = annotHeatmap,
                           annotation_colors = my.colors,
                           angle_col = 315, #tilts the word cluster 
                           color = colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdBu")))(256),
                           main = "Correlation of DE lincRNAs")

#save the output for Yael 
#add the gene names 
annotHeatmap$genename <- Tx.lincRNA$gene_name[match(rownames(annotHeatmap), Tx.lincRNA $gene_id)]

#run this before so that the save is being printed to the right thing. 

dev.off()
for (i in toupper(letters[1:5])){
  print(i)
  IndividualCluster <- subset(annotHeatmap, annotHeatmap$cluster == paste0("cluster_", i))
  #write.csv(IndividualCluster, file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/IndivudalClusterLists/", "cluster", i),
  #          row.names = TRUE, quote = F)
  
  smallHeatmap <- norm_FC_lincRNA[rownames(norm_FC_lincRNA) %in% rownames(IndividualCluster), ]
  
  print(i)
  tiff(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/IndivudalClusterLists/", "HeatmapOfCluster", i, ".tiff"))
  
  print(i)
  ClusterHeatmap <- pheatmap(smallHeatmap ,
           # breaks = seq(from=-thr, to=thr, length=101),
           scale = "row",
           #kmeans_k = 4,
           border_color = NA,
           clustering_method = "complete",
           cluster_cols = F,
           cluster_rows = T,
           show_rownames = F,
           show_colnames = T,
           cellwidth = NA,
           cellheight = NA,
           revC = F,
           treeheight_row = 20,
           treeheight_col = 0,
           color = colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdBu")))(256),
           main = paste("Expression of Cluster", i))
  dev.off()
}