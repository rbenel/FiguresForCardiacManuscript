library(igraph)
library(ggraph)
library(data.table)

#intro to graph theory
#https://www.ebi.ac.uk/training/online/course/network-analysis-protein-interaction-data-introduction/graph-theory-some-basic-definitions


#####################
#Topological Analysis
######################
#load network function - organized here 
#there are notes and explinations, see source code. 
source("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/scripts/Functions_TopolgicalFeatures.R")

#list of features
NetworkFeatures <- data.frame()


for(i in corScore){
  
  #load the major isolate in previous step instead of calculating each time 
  load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/graphScores0.7_1/MajorIsolateGraphs/MajorIsolateGraph_score", i, ".RData"))

  #individual Node features 
  btwnDegreeDF <- NodeFeatures(MainGraphCorr)
  
  #lets get the output the DF and work with it in a different scripts
  save(btwnDegreeDF, file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/networkFeatures/CorScore", i, "features.RData"))
  
  #run the graph features
  IndividualNetwork <- GraphFeatures(MainGraphCorr)
  
  #let's make a table for the Network Features
  if(nrow(NetworkFeatures) == 0)
  
    NetworkFeatures <- IndividualNetwork
  else {
    NetworkFeatures <- rbind(NetworkFeatures, IndividualNetwork)
  }
  
}

#add rownames to the df
rownames(NetworkFeatures) <- corScore

#export to .csv
write.csv(NetworkFeatures, file = "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/FeaturesOfNetworksModules.csv")
#save as RData
save(NetworkFeatures, file = "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/FeaturesOfNetworksModules.RData")
