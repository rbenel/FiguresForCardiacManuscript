#script to set up data for XLmhg 

#load functional lncRNA lists
load(file = "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/Practice 18.12.16/data/ThreeSubsetsFunctlncRNA.RData")


for(i in modules){
  
  load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/networkFeatures/CorScore", i, "features.RData"))
  
  #three attributs we are focusing on
  threeAttributes <- c("NonNormBtwn", "NonNormDegree", "NonNormClose", "pageRank")
  
  for(k in threeAttributes){
  
    databases <- c(FinalGeneralFunctlncRNA, FinalDiseaseCancerIDs)
    names(databases) <- c("GeneralFunctional", "Disease_Cancer")
    
    for(l in c(1:2)){
      
    #order the DF according to each attribute, then convert to list of 0 and 1s for XLmhg
    orderedBtwnDegreeDF <- btwnDegreeDF[order(-btwnDegreeDF[k]),] 
    #we want to rank the ordered list
    orderedBtwnDegreeDF$rankedCol <- rank(-orderedBtwnDegreeDF[[k]], na.last = TRUE,
                                          ties.method = "average")
    
    #find the percent rank or each lincRNA
    #find the length of the list
    lenCol <- nrow(orderedBtwnDegreeDF)
    #divide the ranking by the length of the column, multiply by 100 and round 
    orderedBtwnDegreeDF$percRank <- round((orderedBtwnDegreeDF$rankedCol/lenCol)*100, 2)
    
    print(head(orderedBtwnDegreeDF))
    print(i)
    print(k)
    
    #if the lincRNA found in the functional DB? 
    orderedBtwnDegreeDF$XLmhg <- ifelse(rownames(orderedBtwnDegreeDF) %in% databases[[l]], 1, 0)
    
    
    #output for XLmhg
    XLmhgInput <- orderedBtwnDegreeDF$XLmhg
    
    write.csv(XLmhgInput, file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/XLmhg/SingleModuleUpdatedJuly2021_2SubsetGeneral_Disease/","tiesMethodMin",  i , k, names(databases)[l], ".csv"),
                row.names = F)
     write.csv(orderedBtwnDegreeDF, file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/XLmhg/SingleModuleUpdatedJuly2021_2SubsetGeneral_Disease/OrderedListsOfAttributes/", "TiesMethodMeanAllAttribute", i , k, names(databases)[l], ".csv"), #names(databases)[l]
               row.names = T)
    }
  }
}
