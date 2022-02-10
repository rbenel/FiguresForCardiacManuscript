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

#load data
#updated lnc lincRNA and funcLincRNA annotations
load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/Practice 18.12.16/data/Tx.linc_lncRNA_FunctlncRNA_release92TSL.RData"))

#decided to represent these three lincRNAs as hubs 
ManuallyCuratedHubs <- c("ENSG00000251562", "ENSG00000225783", "ENSG00000249669")
HubGeneNames <- c("MALAT1", "MIAT", "CARMN")


################
#Arrange colors
################

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

#####################
#Correlation scores
######################
#An individual network of lincRNA & lincRNA correlations for hubs 
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

#see table of lincRNAs
table(sigCorr_diag$first_lincRNA)

#for reproducible work 
set.seed(12345)



for(i in c(1:3)){
  
  #Filter the correlation data to include only the correlations with the specific hub we chose
  HubCorrelations <- as_tibble(reshape_upperDiag) %>%
                    dplyr::filter((first_lincRNA == ManuallyCuratedHubs[i] & corScore >= 0.9) |
                    (second_lincRNA == ManuallyCuratedHubs[i] & corScore >= 0.9))
  
  #add the gene names to the lincRNAs which will be used afterwards for labels 
  HubCorrelations$first_lincRNA <- Tx.lincRNA[match(HubCorrelations$first_lincRNA,
                                              Tx.lincRNA[['gene_id']] ), 'gene_name']
  #add the gene names to the lincRNAs which will be used afterwards for labels 
  HubCorrelations$second_lincRNA <- Tx.lincRNA[match(HubCorrelations$second_lincRNA,
                                                     Tx.lincRNA[['gene_id']]), 'gene_name']
  
  #save(HubCorrelations, file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/top1_hubs/ManuallyCuratedHubs/HubLists/HubCorrelationList", HubGeneNames[i], ".RData"))
  ################
  #Make the Graph 
  ###############
  graphCorr <- HubCorrelations %>% graph_from_data_frame(directed = F) 
  
  #visualzie it 
  print(graphCorr)
  
  #check the size of the graph 
  NumEdges <- gsize(graphCorr)
  print(NumEdges)
  
  #open tiff/pdf document so we can write to it at the end of the plot
  #tiff(paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/top1_hubs/ManuallyCuratedHubs/ColoredEdgesLabels", HubGeneNames[i], "_plot.tiff"))

  #remove the margin on the plot to maximize room
  par(mar = c(0,0,0,0))


  #tutorial for color and weight 
  #https://stackoverflow.com/questions/49171958/igraph-edge-width-and-color-positive-and-negative-correlation
  #https://stackoverflow.com/questions/49076441/set-igraph-vertex-color-based-on-vertex-attribute-using-rcolorbrewer
  
  #add edge color and weight to the graph 
  #add weight as an attribute to the graph by utilizing the rounded corScore
  E(graphCorr)$weight <- round(E(graphCorr)$corScore, 2)
  #add width as well by using the weight attribute 
  E(graphCorr)$width <- as.factor(E(graphCorr)$weight)
  #add color as well, this time add the color as a numeric of the factor of weight 
  E(graphCorr)$color <- EdgeColors[[i]][as.numeric(as.factor(cut(E(graphCorr)$weight, 
                                                                 seq(from = 0.89, to = .99, by = 0.02))))]
  
  
  #now add color to the graph 
  #we want the central hub to be blank, and the rest of the nodes colored according to direction (pluri/diff)
  V(graphCorr)$color <- ifelse(V(graphCorr)$name == HubGeneNames[i], "white",  "grey80")
  
  #add a label to the central node and the 10 manually selected nodes
  #use gene name of the hub

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
