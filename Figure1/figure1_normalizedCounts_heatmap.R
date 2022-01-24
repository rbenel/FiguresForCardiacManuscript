library(dplyr)
library(tidyr)
library(pheatmap)
library(reshape)
library(data.table)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)

############
#Load data
############
#updated lncRNA and lincRNA and funcLncRNA
load(file = paste0( "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/Practice 18.12.16/data/Tx.linc_lncRNA_FunctlncRNA_release92TSL.RData"))


#results of norm_counts 
#count data is log2 + pseudocounts (vsd) of TPM...
load(paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/vsd.norm.counts.RData"))

#want to use the DE results from the previous script as a filter 
sigLRT <- read.csv(paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/LRTpadj_cardiac.txt"),
         header = TRUE, sep = ",",row.names = 1)


#read in design information 
samples <- read.table(paste0( "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/metadata_file.csv"), sep = ",", header = TRUE)

samples_name <- as.list(paste0(samples$day, samples$replicate))

###################
#Data Manipulation
###################
colnames(mat2) <- samples_name

mat2 <- tibble::rownames_to_column(as.data.frame(mat2), "id")

mat2_long <- mat2 %>% tidyr::gather(sample, count, 2:ncol(mat2))

mat2_long$sample = gsub("Rep_1|Rep_2|Rep_3", "", mat2_long$sample)

long.df <- mat2_long %>% group_by(id, sample) %>%  summarize(count = mean(count))

tidy.wide <- spread(long.df, key = sample, value = count)

clean.counts <-  as.data.frame(tibble::column_to_rownames(tidy.wide, var = "id"))

clean.counts <- clean.counts[c("day0", "day1", "day2", "day3", "day4", "day5", "day6",
                                     "day8", "day15", "day30", "day60")]

#save(clean.counts, file = paste0(local_path, "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/avgCountsAllData.RData"))

#####################
#Filter for DE genes
#####################
#use LRT filter 
padj_linc <-subset(sigLRT, rownames(sigLRT) %in% Tx.lincRNA$gene_id)

#combine the filter with the normalized data
norm_FC_lincRNA <- subset(clean.counts, rownames(clean.counts) %in% rownames(padj_linc))


#######################
#Create input for MFuzz
#######################
#this could be used for the input for Mfuzz instead of running all of the data. 

#if we want both TF and lincRNA together in same varibale 
norm.FC.TF.linc <- rbind(norm_FC_TF, norm_FC_lincRNA)

#save for later MFuzz
#save(norm.FC.TF.linc, file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/updatedPadjTSLnormFiltLRT_TF.linc.RData"))

#can also save seperately as linc norm counts and TF norm counts 
#save(norm_FC_TF, file = paste0( "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/updatedPadjTSLnormFiltLRT_TF.RData"))
#save(norm_FC_lincRNA, file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/updatedPadjTSLnormFiltLRT_linc.RData"))


#############################
#Heatmaps of lincRNA DE genes 
#############################
#standard lincRNA heatmap 
pheatmap(norm_FC_lincRNA,
         scale = "row",
         border_color = NA,
         clustering_method = "complete",
         cluster_cols = F,
         cluster_rows = T,
         show_rownames = F,
         show_colnames =T,
         cellwidth = NA,
         cellheight = NA,
         revC = F,
         angle_col = 45,
         #cutree_rows = 3,
         color = colorRampPalette( rev(RColorBrewer::brewer.pal(10, "RdBu")))(256))

#Heatmap with names 
#use ComplexHeatmap package

#more consice colnames
colnames(norm_FC_lincRNA) <- c("D0", "D1", "D2", 
                               "D3", "D4", "D5", 
                               "D6", "D8", "D15", 
                               "D30", "D60")


mat_scaled = t(apply(norm_FC_lincRNA, 1, scale))
colnames(mat_scaled) <- colnames(norm_FC_lincRNA)

#names of lincRNAs for our heatmap 
lincRNALabels <- c("TERC", "LINC00261", "RMST", "MALAT1", "CARMN", 
                    "DANCR", "MIAT")

#gene IDS
lincRNAGeneNames <- c("ENSG00000270141", "ENSG00000259974",
                      "ENSG00000255794", "ENSG00000251562", "ENSG00000249669",
                      "ENSG00000226950", "ENSG00000225783")

#match the correct names to the rownames of the data 
ind <- match(lincRNAGeneNames, rownames(mat_scaled))

##############
#Draw Heatmap
##############

ComplexHeatmap::Heatmap(mat_scaled, 
                        show_row_names = FALSE, 
                        show_column_names = T,
                        show_row_dend = FALSE, 
                        show_column_dend = FALSE, 
                        col = colorRampPalette( rev(RColorBrewer::brewer.pal(10, "RdBu")))(256),
                        #col = colorRamp2(seq(-3, 3, by =0.75), rev(RColorBrewer::brewer.pal(9, "RdBu"))),  
                        cluster_columns = FALSE, 
                        column_title_rot = 90,
                        column_names_rot = 45, 
                        heatmap_legend_param = list(title = "", legend_direction = "vertical", 
                                                    color_bar = "continous", labels = c(-3, -2, 0, 2, 3)),
                        clustering_method_rows = "complete", row_dend_side = c("left"), cluster_rows = TRUE) + 
                        rowAnnotation(link = anno_mark(at = ind, labels = lincRNALabels),
                        width = unit(0.05, "cm") + max_text_width(lincRNALabels))

