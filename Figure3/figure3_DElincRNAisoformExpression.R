library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

############
#Load Data
############
#counts in TPM
load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/transcript_results/cts50_oneSample.RData"))

#diff expressed transcripts 
load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/transcript_results/dex50_padjLoweredFilterN16.RData"))

#read in design information 
samples <- read.table(paste0( "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/metadata_file.csv"), sep = ",", header = TRUE)
samples

#updated lnc, linc, and funcLncRNA annotations
load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/Practice 18.12.16/data/Tx.linc_lncRNA_FunctlncRNA_release92TSL.RData"))

#TF annotations
TFs_annot = read.csv(paste0( "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/Practice 18.12.16/data/tf_human_annotations2018.txt"),
                     header = TRUE, sep = ",", stringsAsFactors = F)
################
#Load Functions
################
source(paste0( "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/scripts/flattenCorrMatirxFunction.R"))
strp <- function(x) substr(x,1,15) #function to get rid of version numbers 
source(paste0( "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/functions/summarySE_function.R"))

#######################################
#subsets of all the "types" of transcripts we are interested in 
padjTx_Gene <- subset(dex.padj, dex.padj$transcript <= 0.05) #make sure that also the transcripts are the ones that overcome the 0.05

#get only lincRNA 
padjlincRNA <- subset(padjTx_Gene, padjTx_Gene$geneID %in% Tx.lincRNA$gene_id) #look at lincRNA

########################################################
rownames(cts) <- strp(rownames(cts)) #remove version numbers of count data

#read in design information 
samples <- read.table(paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/metadata_file.csv"), sep = ",", header = TRUE)

samples_name <- as.list(paste0(samples$day, samples$replicate))

###################
#DATA MANIPULATION 
###################
#turn data into avgs of the three replicates
colnames(cts) <- samples_name

cts2 <- tibble::rownames_to_column(as.data.frame(cts), "id")

cts2_long <- cts2 %>% tidyr::gather(sample, count, 2:ncol(cts2))

cts2_long$sample = gsub("Rep_1|Rep_2|Rep_3", "", cts2_long$sample)

long.df <- cts2_long %>% group_by(id, sample) %>%  summarize(count = mean(count))

#this will not work if plyr is loaded
tidy.wide <- tidyr::spread(long.df, key = sample, value = count)

clean.counts <-  as.data.frame(tibble::column_to_rownames(tidy.wide, var = "id"))

clean.counts <- clean.counts[c("day0", "day1", "day2", "day3", "day4", "day5", "day6",
                               "day8", "day15", "day30", "day60")] #to ensure we have the correct order of columns

#convert to log2 for downstream analysis 
clean.counts <- log2(clean.counts + 1) 

################################
#Get counts for DE transcripts
###############################
#counts for updated lincRNA transcripts that are DE 
padjCounts <- subset(clean.counts, rownames(clean.counts) %in% padjlincRNA$txID) #ALL DE TRANSCRIPT 


###########
#Heatmap
##########
#heatmap for the relevant transcripts 
num_clusters <- 3
heatmap <- pheatmap::pheatmap(padjCounts,
                              #kmeans_k = 6, #The function also allows to aggregate the rows using kmeans clustering. This is advisable if number of rows is so big that R cannot handle their hierarchical clustering anymore, 
                              show_rownames  = F, 
                              cluster_cols = F,
                              clustering_distance_rows = "correlation", 
                              scale = "row",
                              color = colorRampPalette( rev(RColorBrewer::brewer.pal(10, "RdBu")))(256),
                              border_color = NA,
                              cutree_rows = num_clusters,
                              treeheight_row = 0,
                              angle_col = 45, #tilts the col labels
                              clustering_method = "complete")

#cut the heatmap tree
heatmap.cut <- cutree(heatmap$tree_row, k = num_clusters)
#choose colors for the heatmap 
my.colors <- viridis::magma(20)

#use the large vector to select 3 colors var apart from each other 
my.colors <- my.colors[c(12, 16, 6)] 


#annotated df of rownames + name of cluster 
annot.df <- data.frame(cluster = (paste("cluster", heatmap.cut, sep = "_")))
rownames(annot.df) <- rownames(padjCounts)

#add colors to the genes 
names(my.colors) <- unique(annot.df$cluster)
my.colors <- list(cluster = my.colors)


#redraw heatmap with annotations and colors 
heatmap.colors <- pheatmap::pheatmap(padjCounts,
                                     show_rownames  = F, 
                                     cluster_cols = F,
                                     clustering_distance_rows = "correlation", 
                                     scale = "row",
                                     color = colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdBu")))(256),
                                     border_color = NA,
                                     clustering_method = "complete",
                                     annotation_row = annot.df,
                                     annotation_colors = my.colors,
                                     annotation_legend = T,
                                     treeheight_row = 0,
                                     angle_col = 315, #tilts the col labels 
                                     cex = 1,
                                     cutree_rows = num_clusters)

####################
#DF of the clusters
####################
#add more info to the annot.df 
annot.df <- merge(annot.df, padjlincRNA, by.x = "row.names", by.y = "txID")

#if we turn the cluster column into a character 
#than we can sort by order and maybe we wont need all of the gsub?
#this seems to work
annot.df$cluster <- as.character(annot.df$cluster)
annot.df <- annot.df[order(annot.df$cluster), ]
names(annot.df)[names(annot.df) == "Row.names"] <- "txID"


#add the gene name to the results table 
annot.df$GeneName <- Tx.lncRNA[match(annot.df$geneID, Tx.lncRNA[["gene_id"]]), "gene_name"]

#save as .RData
#save(annot.df, file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/transcript_results/March20_OnlyDTUtranscriptsAnnotN16WithGeneName.df.RData"))
