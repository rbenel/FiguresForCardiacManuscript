#a global picture of lincRNA V. protein coding genes expressed during the experiment 

library(dplyr)
library(tidyr)
library(pheatmap)
library(reshape)
library(data.table)
library(DESeq2)
library(ggplot2)
library(viridis)
library(ggpubr)
library(rstatix)

#The count data here is log2 + pseudocounts (vsd) of TPM...

#######################
#LOAD GENE ANNOTATIONS 
#######################
#read in the annotations compiled in a different script
#based on annotations from the ensmbl package 

#lists for lncRNA, lincRNA, and functional lincRNAs can be loaded from here. 
load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/Practice 18.12.16/data/Tx.linc_lncRNA_FunctlncRNA_release92TSL.RData"))

#a list of protein coding genes can be loaded here 
load(paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/Practice 18.12.16/data/proteinCoding_release92TSL.RData"))


######################
#LOAD DATA AND DESIGN
######################
#load the results of normalized counts for the experiment saved in previous script 
load("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/vsd.norm.counts.RData")

#read in design information for this exepriment 
samples <- read.table(paste0( "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/metadata_file.csv"), sep = ",", header = TRUE)

samples_name <- as.list(paste0(samples$day, samples$replicate))

###################
#DATA MANIPULATION
###################
#here turn the data from wide to long to find the mean of each day, then regroup and turn back to wide 
colnames(mat2) <- samples_name

mat2 <- tibble::rownames_to_column(as.data.frame(mat2), "id")

mat2_long <- mat2 %>% tidyr::gather(sample, count, 2:ncol(mat2))

mat2_long$sample = gsub("Rep_1|Rep_2|Rep_3", "", mat2_long$sample)

long.df <- mat2_long %>% group_by(id, sample) %>%  summarize(count = mean(count))

tidy.wide <- spread(long.df, key = sample, value = count)

clean.counts <-  as.data.frame(tibble::column_to_rownames(tidy.wide, var = "id"))

#rearrange the columns to non alphabetical order 
clean.counts <- clean.counts[c("day0", "day1", "day2", "day3", "day4", "day5", "day6",
                               "day8", "day15", "day30", "day60")]

###########################
##COUNTS FOR EACH CATEGORY
###########################
#lincRNAs 
lincRNAgenes <- subset(clean.counts, rownames(clean.counts) %in% Tx.lincRNA$gene_id)
#PCG
proteinCoding <- subset(clean.counts, rownames(clean.counts) %in% Tx.proteinCoding$gene_id)


#add type to each df so we can seperate them into groups. 
lincRNAgenes <- data.frame(type = "lincRNA", lincRNAgenes)

proteinCoding <- data.frame(type = "protein_coding", proteinCoding)

TotalGenes.df <- do.call(rbind, list(lincRNAgenes, proteinCoding))


#to find out the counts of each group for the legend of the boxplot
table(TotalGenes.df$type)

TotalGenes.df <- tibble::rownames_to_column(TotalGenes.df, "id")

#long data.. with the day, count, and type 
TotalGenes.long <- TotalGenes.df %>% tidyr::gather(days, count, 3:ncol(TotalGenes.df))

#so the days appear in order.. need to rearrange the factor 
TotalGenes.long$days <-factor(TotalGenes.long$days, levels=c("day0", "day1", "day2", "day3", "day4", "day5", "day6",
                                                       "day8", "day15", "day30", "day60"))


#stat test before we add pvalues to plot 
#this will give us the option to add manual pvalues, bec the padj and the ns dont show 
#are not compatible. this way should work much better!
stat.test <- TotalGenes.long %>%
  group_by(type) %>%
  wilcox_test(count ~ days, ref.group = "day0") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() 
 
stat.test <- stat.test %>%
  mutate(y.position = c(1, 21, rep(1,5), 22, 23, 24, rep(1,4), 25.5, rep(1, 3), 26.5, 27.5))
#filter(p.signif != "ns")

#write.csv(stat.test, file = "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/figures_paper/Figure1/March22_2020_figure1/StatCompareDataWilcoxTestRefDay0.csv",
#          row.names = F)

#boxplot of all days
BoxPlotSideBySide <- ggplot(TotalGenes.long, aes(x=days, y = count, fill = type, color = days)) +
                    geom_boxplot(color="black",position = position_dodge2(padding = .1), notch = T) +
                    labs(y = "count (log2)", x = " ") + 
                    theme_classic() + 
                    scale_fill_manual(values = c('#00BFC4', '#AB2F5EFF')) + #F8766D')) +  
                    #stat_compare_means(aes(group = days, label = ..p.signif..), method = "wilcox.test", 
                    #ref.group = "day0", hide.ns = TRUE) +
                    stat_pvalue_manual(stat.test,
                                       label  = "p.adj.signif", hide.ns = TRUE) +
                    # legend
                    theme(legend.title = element_blank()) +
                    theme(legend.text = element_text(colour="black", size = 20, face = "bold")) +
                    theme( axis.title.x = element_text(family="sans",size = 18, face="bold", hjust=0.5, vjust=-0.5),
                           axis.title.y = element_text(family="sans",size = 16, angle=90, face="bold", hjust=0.5, vjust=1)) +
                    theme( axis.text.x = element_text(family = "sans",size = 18, angle=45, face='plain', colour="#353535",   hjust=1, vjust=1) ) +
                    theme( axis.text.y = element_text(family = "sans",size = 18, face='plain', colour="#353535",  vjust=0.5) ) +
                    theme(axis.line.x = element_line(color="black", size = 0.5),
                          axis.line.y = element_line(color="black", size = 0.5)) +
                    theme(legend.background = element_rect()) + 
                    theme(legend.position="top") +
                    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),  panel.background = element_blank()) +
                    scale_x_discrete(labels = c("day0" = "D0", "day1" = "D1", "day2" = "D2", 
                                            "day3" = "D3", "day4" = "D4", "day5" = "D5", 
                                            "day6" = "D6", "day8" = "D8", "day15" = "D15", 
                                            "day30" = "D30", "day60" = "D60"))

BoxPlotSideBySide  

ggsave("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/figures_paper/Figure1/March22_2020_figure1/WithLegendSideBySideStatCompareDataWilcoxTestRefDay0Plot.tiff")

