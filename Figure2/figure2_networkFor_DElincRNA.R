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

###########
#Load Data
##########

#results from DESeq2 for Day 0 v Day 60 -- Wald Test 
load(file = "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/resDay0Day60.RData")

#turn into a data frame and decide in which direction the fold change is in
res_df <- as.data.frame(resWald)
#add rownames?
res_df$rownames <- rownames(res_df)

#add FC direction for coloring of the correlation network
#Unregulated genes at Day 60 are positive FC numbers 

#one liner to get a vector of the status of all FCs
res_df$FCstatus <- cut(res_df$log2FoldChange, breaks = c(min(res_df$log2FoldChange), -1, 
                                            1, max(res_df$log2FoldChange)), labels = c("downregulated", "neither", "upregulated"))


###########
#NETWORK
###########
#see here for a tutorial 
#in the example they use correlate and stretch(), which removes NAs and the upper triangle
#https://drsimonj.svbtle.com/how-to-create-correlation-network-plots-with-corrr-and-ggraph

#load data if you want all of the lincRNAs together... 
load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/April20_lincRNA_lincRNA_CorrScores.RData"))

###############################################
#Filter out the pairs that are one to another
###############################################
#to remove the diagnol use diag = TRUE 
#https://stackoverflow.com/questions/26377199/convert-a-matrix-in-r-into-a-upper-triangular-lower-triangular-matrix-with-those
correlation_spearman[upper.tri(correlation_spearman, diag = TRUE)] <- NA

#lets reshape the data for easy sub-setting. 
reshape_upperDiag <- reshape::melt(correlation_spearman) 

#rename the columns 
colnames(reshape_upperDiag) <- c("first_lincRNA", "second_lincRNA", "corScore")

#remove "NAs"
CleanReshape_upperDiag <- na.omit(reshape_upperDiag)

#now we want only the groups that are over 0.9? 
sigCorr_diag <- subset(CleanReshape_upperDiag, corScore >= 0.9)

table(sigCorr_diag$first_lincRNA)


################
#Filter Network
################
#filter the graph and save multiple options
set.seed(12345)

#in the end we chose a a filter of 0.9! 
#here on in this is the network we will be using 
corScores_forFilter <- seq(from = 0.7, to = 0.99, by = 0.05)

for(i in corScores_forFilter){
    
    #create the filtered graph, using only positive correlations 
    graphCorr <- as_tibble(reshape_upperDiag) %>% dplyr::filter(corScore >= i) %>% 
      graph_from_data_frame(directed = F)
    
    print(graphCorr)
   
    #save the graph
    #save(graphCorr, file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/graphScores0.7_1/GraphlincRNA_lincRNA_score", i, ".RData"))
    
    #to learn how to label only some genes in final visual see here :) 
    #https://stackoverflow.com/questions/47175541/plot-labels-for-specific-geom-node-text
    
    #open tiff/pdf document so we can write to it at the end of the plot
    tiff(paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/graphScores0.7_1Visual/BW_GraphVisual_", i, "_plot.tiff"))

    graph <- ggraph::ggraph(graphCorr) +
      geom_edge_link(color = "grey80", edge_width = 1) +
      geom_node_point(color = "black", size = 1.5) +
      theme_graph() +
      labs(x = "", y = "") +
      labs(title = paste0("lincRNA lincRNA Correlation Network (", i, ")"))
      #to print graph
      print(graph)
      dev.off()
      
      #############################
      #Calculate the Major Isolate
      #############################
      #As we can see from the visualization of the graph, there are some components that aren't connected. 
      #We will restrict our analysis only on the biggest component
      
      #calculate the maximal connected components of a graph
      graphComp <- components(graphCorr)
      
      #next lets find out which component is the max
      #csize = numeric vector giving the sizes of the clusters.
      LargestComp <- which.max(graphComp$csize)
      
      #Major isolate of the each graph
      MajorIsolate <- graphComp$csize[LargestComp]
      print(MajorIsolate)
      
      #subset the graph to include only the largest component
      MainGraphCorr <- induced_subgraph(graphCorr, which(graphComp$membership == LargestComp))
    
      
      #to learn how to add attributes to the network look here. 
      #https://rpubs.com/niyer/389851
      #matching the FCstatus with the nodes by including it as an attribute 
      V(MainGraphCorr)$FCstatus <- as.character(res_df$FCstatus[match(V(MainGraphCorr)$name, res_df$rownames)])

      
      #lets save this major isolate as well instead of calculating again
      #save(MainGraphCorr, file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/graphScores0.7_1/MajorIsolateGraphs/MajorIsolateGraph_score", i, ".RData"))

      print(MainGraphCorr)
      #open tiff/pdf document so we can write to it at the end of the plot
      tiff(paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/graphScores0.7_1Visual/Graphs0.9/RdBlu_MajorIsolate_Colored_GraphVisual_", i, "_plot.tiff"))

      IsolateGraph <- ggraph::ggraph(MainGraphCorr) +
        geom_edge_link(color = "grey80", edge_width = 1) +
        #geom_node_point(color = "black", size = 1.5) + 
        geom_node_point(aes(color = FCstatus), show.legend = F) + #turn to factor
        #scale_color_manual(values = c("#451077FF", "#4b4c4c", "#F1605DFF")) +
        scale_color_manual(values = c("#4393C3", "#4b4c4c", "#D6604D")) +
        theme_graph() +
        labs(x = "", y = "") + 
        labs(title = paste0("Major Isolate lincRNA lincRNA Correlation Network (", i, ")"))
      #to print graph
      print(IsolateGraph)
      dev.off()
      
      
}
