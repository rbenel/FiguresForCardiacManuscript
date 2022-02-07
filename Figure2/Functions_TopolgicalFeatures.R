#build functions that calculates various graph features 

#see some sources here 
#https://www.youtube.com/watch?v=K2WF4pT5pFY
#https://www.e-education.psu.edu/geog597i_02/node/832


#the first one is for features that are PER node
NodeFeatures <- function(x){
  
  
  #Updated March 2021
  load(file = "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/Practice 18.12.16/data/NoCRISPRupdatedFunc_Tx.linc_lncRNA_FunctlncRNA_release92TSL.RData")
  
  ############
  #Betweenness
  ############
  #It is basically defined as the number of shortest paths in the graph that pass through the node 
  #Nodes with a high betweenness centrality are interesting because they lie on communication paths 
  #and can control information flow
  #divided by the total number of shortest paths.
  
  #Bnorm=2*B/(n*n-3*n+2)
  #where Bnorm is the normalized, B the raw betweenness, and n is the number of vertices in the graph.
  #Betweenness: Put simply, with normalization, the result is divided by the number of ordered or 
  #unordered vertex pairs in directed and undirected graphs, respectively.
  
  #calculate the betweenness centrality measures how important a node is to the shortest paths through the network.
  #this is a better "global" measure of the node because it takes into account the network size. 
  btwn <- data.frame(betweenness(MainGraphCorr, v = V(MainGraphCorr), directed = FALSE, normalized = TRUE))
  nonNormBtwn <- data.frame(betweenness(MainGraphCorr, v = V(MainGraphCorr), directed = FALSE, normalized = FALSE))
  
  ############
  #Degree
  ############
  #calculate the degree or the number of its adjacent edges
  #number of connections each node has == degree 
  
  #Degree: The result is divided by the highest possible degree in a simple graph, 
  #i.e. one less than the number of vertices.
  
  degree <- data.frame(degree(MainGraphCorr, v = V(MainGraphCorr), normalized = TRUE))
  nonNormDegree <- data.frame(degree(MainGraphCorr, v= V(MainGraphCorr), normalized = FALSE))
  
  
  ###########
  #Closeness
  ###########
  #Closeness centrality measures how short the shortest paths are from node i to all nodes. 
  #It is usually expressed as the normalised inverse of the sum of the topological distances in the graph
  #using the adjacent matrix, calculate the shortest path btwn each two nodes, then row sum. 
  #Calculate: num total nodes -1/row sum 
  
  #Logical scalar, whether to calculate the normalized closeness. Normalization is performed by multiplying 
  #the raw closeness by n-1, where n is the number of vertices in the graph.
  
  #Here igraph uses the inverse of the average length of the shortest path to/from all other vertices in the graph 
  #Therefore Higher scores are given to the nodes that appear more central in terms of distance 
  #(they can reach other in a few hops)
  #A vertex with a high closeness centrality would mean it has close relationships with many vertices
  
  #this warning is when we have disconnected elements, before taking largest component 
  #what happens when you try to calculate closeness for a graph that is disconnected..
  #WARNING: closeness centrality is not well-defined for disconnected graphs
  
  #Closeness: The result is divided by the number of “other vertices” that a distance is calculated 
  #to, i.e. one less than the total number of vertices.
  
  close <- data.frame(closeness(MainGraphCorr, v = V(MainGraphCorr), normalized = TRUE))
  nonNormClose <- data.frame(closeness(MainGraphCorr, v = V(MainGraphCorr), normalized = FALSE))
  
  ####################
  #Cluster Coefficient
  ####################
  #How likely are two nodes that are connected are part of a larger clique? 
  #Look at a particular node, and then count its degrees. 
  #Then look at the number of neighbour nodes that are connected to each other! 
  #Equation to calculate: 2 * Numb of neighbours that are connected/Num of degrees * (Num of degrees -1)
  #Possibe interactions between the neighbours of a particular node --> always a fraction between 0 and 1
  
  #this calculates local transitivity -- isolates can be Na or zero
  #If NaN then local transitivity is not reported in the averaging, defines how to treat vertices with degree zero and one
  #if there are no vertices with degree two or higher
  
  #a local clustering coefficient of 0 "Finally, none of the possible connections among the neighbours of the 
  #blue node are realised, producing a local clustering coefficient value of 0.
  #depends if we want them calculauted in the global measure 
  clusterCoef <- data.frame(transitivity(MainGraphCorr, type = "local", v = V(MainGraphCorr), isolates = "NaN"))

  ############
  #Page Rank 
  ############
  #calculates the google page rank for specified vertices
  pageRank <- data.frame(page_rank(MainGraphCorr, vids = V(MainGraphCorr))$vector)
  
  
  #combine all of the outputs
  #get DF and add rownames for labels 
  btwnDegreeDF <- setnames(cbind(btwn, nonNormBtwn, degree, nonNormDegree, 
                                 close, nonNormClose, clusterCoef, pageRank, rownames(btwn)), c("btwn", "NonNormBtwn", "degree", "NonNormDegree", "close", "NonNormClose", "clusterCoef", "pageRank", "names"))
  
  #change to external names, because too difficult to read ensembl IDS
  btwnDegreeDF$names <- Tx.lincRNA[match(btwnDegreeDF$names, Tx.lincRNA[['gene_id']]), 'gene_name']
  
  #if we just want to label the functional lincRNAs
  #this allows us to label just markers on our plot 
  btwnDegreeDF$funcLincRNA <- as.factor(ifelse(btwnDegreeDF$names %in% Tx.lncRNAfunctional$gene_name, paste(btwnDegreeDF$names), NA))
  
  return(btwnDegreeDF)
}


#the second function is for the entire graph enviroment 
GraphFeatures <- function(x){
  #how many edges does the graph have?
  NumEdges <- gsize(MainGraphCorr)
  #print(NumEdges)
  
  #how many nodes/vertices does the graph have?
  NumNodes <- gorder(MainGraphCorr)
  #print(NumNodes)
  
  #find the avg degree for the graph? 
  #need to first calcualte the degree for each node to get avg
  avgDegree <- mean(degree(MainGraphCorr, v= V(MainGraphCorr), normalized = FALSE))
  #print(avgDegree)
  
  #calculates the Global Cluster Coef i.e. for the entire graph
  GlobalClusterCoef <- transitivity(MainGraphCorr, type = "global", v = V(MainGraphCorr), isolates = "NaN")
  #print(GlobalClusterCoef)
  
  #calculates the Avg Cluster Coef i.e. for the entire graph
  AvgClusterCoef <- transitivity(MainGraphCorr, type = "average", v = V(MainGraphCorr), isolates = "NaN")
  #print(AvgClusterCoef)
  
  #We also want to look at average shortest path length
  #an distance and avg length receive the same value
  avgPathLength <- average.path.length(MainGraphCorr)
  #print(avgPathLength)
  
  #calculates the average path length in a graph, by calculating the shortest paths between all pairs of vertices 
  #(both ways for directed graphs). This function does not consider edge weights currently and uses a breadth-first search
  #meanDist <- mean_distance(MainGraphCorr)
  #print(paste("meanDist", meanDist, sep = " "))
  
  #build DF with the Network Features 
  IndividualNetwork <- data.frame(NumEdges, NumNodes, avgDegree, AvgClusterCoef, avgPathLength)
  
  return(IndividualNetwork)
}

