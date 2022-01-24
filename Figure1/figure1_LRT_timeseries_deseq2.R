library(DESeq2)
library(pheatmap)
library("RColorBrewer")
library(dplyr)
library(tidyr)
library(tibble)


#see tutorial here for time-series RNAseq analysis 
#https://www.bioconductor.org/help/course-materials/2015/LearnBioconductorFeb2015/B02.1.1_RNASeqLab.html

#############################
#Time-series experiments LRT
#############################
#test for DE using the likelihood ratio test as described in the following section

#ensembl version
load(paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/ddsTxi_obj.RData"))

#LRT compares two model, and the difference between them will be what you're statistically testing.

ddsTxi_lrt <- DESeq(ddsTxi, test="LRT", reduced = ~ replicate)

#see results... 
resultsNames(ddsTxi_lrt)

resLRT <- results(ddsTxi_lrt)
summary(resLRT)

head(resLRT[order(resLRT$pvalue),],4)

#save all results 
#save(resLRT, file = paste0(local_path, "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/allLRT_cardiac.RData"))

#############
#Plot counts
#############
# We can plot the counts for the groups over time using ggplot2, for the **gene** with the smallest p value, 

sigGene <- plotCounts(ddsTxi_lrt, which.min(resLRT$padj), 
                           intgroup = c("replicate","day"), returnData = TRUE)

sigGeneGgplot <- ggplot(sigGene,
            aes(x = as.numeric(day), y = count, color = day, group = replicate)) + #, label = samples$sampleID)) + 
  geom_point() + geom_smooth(se = FALSE, method = "loess") + scale_y_log10() 

sigGeneGgplot +  scale_color_manual(values=c("#f9bdbd", "#f74a4a", "#648cea", "#ffaa56", "#f9e67a", "#caa9e5", "#a6dd6a", "#05aa60", 
                                 "#b0e5df","#9a25d1" , "#4bb7e5")) + labs(x = "days", y = "count (log10)")


##################
#Output Sig genes
##################
#In LRT  genes with small p values from this test are those which at one 
#or more time points after time 0 showed a strain-specific effect. Therefore we can take all of
#the genes that have a padj of < 0.05.. and this can be our filtered list of DE genes.. 

#REMINDER:
#Although the output includes FC, this can not be taken into account for the entire comparison
#as it represents the FC from day0 to day60 whereas padj reflects a difference at any time point

sig.LRTres <- subset(as.data.frame(resLRT), padj < 0.05)
#write.csv(sig.LRTres, paste0(local_path, "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/LRTpadj_cardiac.txt"))

