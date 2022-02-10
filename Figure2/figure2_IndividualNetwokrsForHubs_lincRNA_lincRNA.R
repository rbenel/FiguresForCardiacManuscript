#May 17th 2020 
#make individual networks for top 1% "hubs" i.e. lincRNAs with highest degree. 

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
library(RColorBrewer)

#load top degree output
#load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/networkFeatures/OnlyTopDegree/", "0.9", "Top", "10", "perc", ".RData"))

#April 2020: updated lnc lincRNA and funcLincRNA annotations
load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/Practice 18.12.16/data/Tx.linc_lncRNA_FunctlncRNA_release92TSL.RData"))

#decided to represent these three lincRNAs
ManuallyCuratedHubs <- c("ENSG00000251562", "ENSG00000225783", "ENSG00000249669")
HubGeneNames <- c("MALAT1", "MIAT", "CARMN")
#HubColors <- c("#D6604D", "#4393C3", "#D6604D")

#if we want to color the edges and not the nodes we need this 
RedColor <- brewer.pal(n = 9, name = "Reds")
CustomRed <- RedColor[c(3, 5, 6, 7, 9)]

BlueColor <- brewer.pal(n = 9, "Blues")
CustomBlue <- BlueColor[c(3, 5, 6, 7, 9)]
  
EdgeColors <- list(CustomRed, CustomBlue, CustomRed)

#added labels to certain nodes in the network. 
#Change the network to be based off gene names, to make this more reproducible 
MALATlabelsGeneName <- c("MALAT1", "AC010894.2", "LINC02016", "AC024361.2", "LINC00963", 
                        "AL513303.1", "LINC01290", "LINC01219", "CARMN", "BCAR4", "CAHM")


MIATlabelsGeneName <- c("MIAT", "SNHG19", "LINC02582", "AC020928.2", "AC002091.1", "LINC01108",
                        "SNHG3", "LNCPRESS1", "MIR124-2HG", "SNHG5", "LINC00545")
  

CARMNlabelsGeneName <- c("CARMN", "LINC01290", "BX255923.1", "C1orf143", "LINC00964", "LINC02274", 
                         "MALAT1", "LINC01324", "MIR99AHG", "LINC01940", "LINC01128")

NodeGeneNames <- list(MALATlabelsGeneName, MIATlabelsGeneName, CARMNlabelsGeneName)

########
#GRAPH?
########
#February 12th 2020 -> Make a Graph of Heatmap all lincRNA & lincRNA correlations
#see here for instructions 
#in the example they use correlate and stretch(), which removes Nas and the upper triangle
#https://drsimonj.svbtle.com/how-to-create-correlation-network-plots-with-corrr-and-ggraph

#load data 
load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/April20_lincRNA_lincRNA_CorrScores.RData"))

###############################################
#Filter out the pairs that are one to the other
###############################################
#to remove the diagnol use diag = TRUE 
#https://stackoverflow.com/questions/26377199/convert-a-matrix-in-r-into-a-upper-triangular-lower-triangular-matrix-with-those
correlation_spearman[upper.tri(correlation_spearman, diag = TRUE)] <- NA

#lets reshape the data for easy subsetting. 
#can't figure out how to do this properly using dplyr 
reshape_upperDiag <- reshape::melt(correlation_spearman) 

colnames(reshape_upperDiag) <- c("first_lincRNA", "second_lincRNA", "corScore")

CleanReshape_upperDiag <- na.omit(reshape_upperDiag)

#now we want only the groups that are over 0.9? 
sigCorr_diag <- subset(CleanReshape_upperDiag, corScore >= 0.9)

table(sigCorr_diag$first_lincRNA)

################
#Filter Network
################
#filter the graph and save mutliple options
set.seed(12345)
#top 1% degrees 
#Hubs <- c(rownames(OnlyTopDegree))

#not reporducible merging fromfigure_graphs_delincRNA for resDF  - yael wants list of hubs with up and down
# HubsDEres <- subset(res_df, rownames %in% Hubs)
# Top5percHubsUpDown <- merge(OnlyTopDegree, HubsDEres, by = "row.names")
# write.csv(Top5percHubsUpDown, file = "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/top1_hubs/Top10percHubsUpDown.csv", 
#           row.names = F)

for(i in c(1:3)){
  
  # #On June 8th I added this... see explanation
  # #create output to manually chose labels and to label them this might be more reproducible?
  
  HubCorrelations <- as_tibble(reshape_upperDiag) %>%
                    dplyr::filter((first_lincRNA == ManuallyCuratedHubs[i] & corScore >= 0.9) |
                    (second_lincRNA == ManuallyCuratedHubs[i] & corScore >= 0.9))
  #change to external name so yael can chose the labels
  HubCorrelations$first_lincRNA <- Tx.lincRNA[match(HubCorrelations$first_lincRNA,
                                              Tx.lincRNA[['gene_id']] ), 'gene_name']
  #change to external name so yael can chose the labels
  HubCorrelations$second_lincRNA <- Tx.lincRNA[match(HubCorrelations$second_lincRNA,
                                                     Tx.lincRNA[['gene_id']]), 'gene_name']
  
  #yael wants in the table to be hubs, corscore, chr#, tx start and tx end. 
  #save(HubCorrelations, file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/top1_hubs/ManuallyCuratedHubs/HubLists/HubCorrelationList", HubGeneNames[i], ".RData"))
  # write.csv(HubCorrelations, file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/top1_hubs/ManuallyCuratedHubs/HubLists/HubCorrelationList", HubGeneNames[i], ".csv"))
  
  #create the filtered graph, only want positive correlations
  # graphCorr <- as_tibble(reshape_upperDiag) %>%
  #   dplyr::filter((first_lincRNA == ManuallyCuratedHubs[i] & corScore >= 0.9) | 
  #                   (second_lincRNA == ManuallyCuratedHubs[i] & corScore >= 0.9))  %>%
  #   graph_from_data_frame(directed = F) 
    
  graphCorr <- HubCorrelations %>% graph_from_data_frame(directed = F) 
  
  print(graphCorr)
  
  #save the graph
  #save(graphCorr, file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/graphScores0.7_1/Graph", i, ".RData"))
  
  #U = undirected 
  #N = names graph 
  #the number of vertices and then the number of edges in the graph
  #vertices are total number of nodes
  
  #attr: name (v/c), corScore (e/n) 
  #v/c is the vertex level character attribute and e/n is edge level numeric attribute 
  
  #the entities are nodes or vertices, conenctions are edges or links
  
  #igraph is based on adjacency matries --> square matirx where column and rownames are the nodes of the network 
  #Edge List: translates the two "from" and "to" columns to binary matrix Y/N connection for each of the nodes?
  #edge list -> one column of nodes that are the source of a connection and the other column of nodes that are 
  #targets of the connection. 
  
  #example with the letters and cities:
  #source = cities where the correspondents wrote letters
  #destination = cities where Daniel received the letters 
  
  #node list = all cities both from source and destination 
  
  #edge list can also contain other columns that describe attributes of the edge i.e. magnitude i.e. weighted 
  
  NumEdges <- gsize(graphCorr)
  print(NumEdges)
  
  #this is were the idea of how to label only TFs came from :) 
  #https://stackoverflow.com/questions/47175541/plot-labels-for-specific-geom-node-text
  
  #open tiff/pdf document so we can write to it at the end of the plot
  #remove the margin on the plot to maximize room
  
  #tiff(paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/top1_hubs/ManuallyCuratedHubs/ColoredEdgesLabels", HubGeneNames[i], "_plot.tiff"))
  #pdf(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/top1_hubs/ManuallyCuratedHubs/ColoredEdgesLabels", HubGeneNames[i], "_plot.pdf"))
  
  par(mar = c(0,0,0,0))
  #June 1st 2020 - ggrpah isnt working. Will use igraph instead?
  # graph <- ggraph::ggraph(graphCorr) +
  #   geom_edge_link(color = "grey80", edge_width = 1) +
  #   geom_node_point(color = "black", size = 1.5) +
  #   theme_graph() +
  #   #labs(x = "", y = "") +
  #   labs(title = paste0(i, " Hub of the 0.9 Correlation Network"))
  # #to print graph
  # print(graph)
  #dev.off()
  
  
  #add edge color and weight to the graph 
  #check this out for color and weight 
  #https://stackoverflow.com/questions/49171958/igraph-edge-width-and-color-positive-and-negative-correlation
  #https://stackoverflow.com/questions/49076441/set-igraph-vertex-color-based-on-vertex-attribute-using-rcolorbrewer
  #add weight as an attribute to the graph by utilizing the rounded corScore
  E(graphCorr)$weight <- round(E(graphCorr)$corScore, 2)
  #add width as well by using the weight attribute 
  E(graphCorr)$width <- as.factor(E(graphCorr)$weight)
  #add color as well, this time add the color as a numeric of the factor of weight 
  E(graphCorr)$color <- EdgeColors[[i]][as.numeric(as.factor(cut(E(graphCorr)$weight, 
                                                                 seq(from = 0.89, to = .99, by = 0.02))))]
  
  
  #add color to the graph 
  #we want the central hub to be blank, and the rest of the nodes colored according to direction (pluri/diff)
  V(graphCorr)$color <- ifelse(V(graphCorr)$name == HubGeneNames[i], "white",  "grey80") #HubColors[i])
  
  #add a label to the central node and the 10 manually selected nodes
  #use gene name of the hub
  #V(graphCorr)$label <- ifelse(V(graphCorr)$name == ManuallyCuratedHubs[i], HubGeneNames[i], "")
  
  V(graphCorr)$label <- ifelse(V(graphCorr)$name %in% NodeGeneNames[[i]], V(graphCorr)$name,  "")
  #change label size so we can read it 
  V(graphCorr)$label.cex <- 1.5
  
  V(graphCorr)$label.color <- "black"
  
  V(graphCorr)$label.font <- 2
  
  igraphGraph <- plot.igraph(graphCorr, vertex.size = 15, 
                      width = 8, edge.arrow.width = 2, asp = 1) #margin = -0.2
  
  print(igraphGraph)
  
  dev.off()
}
