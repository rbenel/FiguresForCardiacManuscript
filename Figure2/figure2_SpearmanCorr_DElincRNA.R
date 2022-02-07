library(pheatmap)
library(data.table)
library(corrr)

#all DE lincRNA
load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/updatedPadjTSLnormFiltLRT_linc.RData"))


#Updated March 2021
load(file = "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/Practice 18.12.16/data/NoCRISPRupdatedFunc_Tx.linc_lncRNA_FunctlncRNA_release92TSL.RData")

##################################
#Spearman Correlation all lincRNA
##################################

correlation_spearman <-cor(t(norm_FC_lincRNA), t(norm_FC_lincRNA),  method = c("spearman"))

#turn into a matrix 
correlation_spearman = as.matrix(correlation_spearman)

#color scheme (blue/red)
color <- colorRampPalette(c("#579af2", "white", "#cc1a35"))

#view the distribution of the scores 
hist <- hist(correlation_spearman, breaks=20, xlab= "Correlation Values", main="Distribution of Spearman Correlation",
             col= color(20))

pheatmap(correlation_spearman ,
         scale = "none",
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

############################
#Heatmap with Clusters
############################
#Cut the graph into X number of clusters and add annotations  
numClusters <- 5

#heatmap with clusters and annotations 
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
                   main = "Correlation of DE lincRNAs")
      

############
#Annotation 
############
#cut the heatmap tree
heatmapCutree <- cutree(lincRNAheatmap$tree_row, k = numClusters)

#turn into letters so less confusing in text
heatmapCutree <- sapply(heatmapCutree, function(i) toupper(letters[i]))

#choose colors for the annotations of the clusters 
heatmapColors <- RColorBrewer::brewer.pal(5, 'Set2')

#annotated df of rownames + name of cluster 
annotHeatmap <- data.frame(cluster = paste("cluster", heatmapCutree, sep = "_"))
rownames(annotHeatmap) <- rownames(correlation_spearman)

#create a df with functional and lincRNA to add as an annotation 
FunctionalLincRNA <- subset(Tx.lncRNAfunctional, Tx.lncRNAfunctional$tx_biotype == "lincRNA")

#see how many are in each cluster
table(annotHeatmap)

#how many are functional
length(intersect(rownames(annotHeatmap), Tx.lncRNAfunctional$gene_id))

#add functional to be part of the annotation DF
annotHeatmap$functional <- as.factor(ifelse(rownames(annotHeatmap) %in% Tx.lncRNAfunctional$gene_id, "func", ""))

#annotHeatmap$functional <- relevel(annotHeatmap$functional, ref = 0)
#assign colors to functional 
functionalColors = c( "white", "black")

#look at the table of the annotation
table(annotHeatmap)

#add colors to the genes 
names(heatmapColors) <- unique(annotHeatmap$cluster)
names(functionalColors) <- unique(annotHeatmap$functional)
my.colors <- list(cluster = heatmapColors,
                  functional = functionalColors)

########################
#Annotation (continued)
########################
lincRNAheatmapColored <- pheatmap(correlation_spearman ,
                           scale = "none",
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

#add the gene names 
annotHeatmap$genename <- Tx.lincRNA$gene_name[match(rownames(annotHeatmap), Tx.lincRNA $gene_id)]
