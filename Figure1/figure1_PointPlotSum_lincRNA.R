#point plot of expression (at least 20 counts) of lincRNA across days 
library(dplyr)
library(tidyr)
library(pheatmap)
library(reshape)
library(data.table)
library(DESeq2)
library(ggplot2)
library(viridis)

###########
#Load Data 
###########

#updated lnc lincRNA and funcLincRNA annotations
load(file = paste0( "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/Practice 18.12.16/data/Tx.linc_lncRNA_FunctlncRNA_release92TSL.RData"))


#results of norm_counts 
load(paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/vsd.norm.counts.RData"))

#read in design information 
samples <- read.table(paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/metadata_file.csv"), sep = ",", header = TRUE)

samples_name <- as.list(paste0(samples$day, samples$replicate))

##################
#Data Manipulation
##################
colnames(mat2) <- samples_name

mat2 <- tibble::rownames_to_column(as.data.frame(mat2), "id")

mat2_long <- mat2 %>% tidyr::gather(sample, count, 2:ncol(mat2))

mat2_long$sample = gsub("Rep_1|Rep_2|Rep_3", "", mat2_long$sample)

long.df <- mat2_long %>% group_by(id, sample) %>%  summarize(count = mean(count))

tidy.wide <- spread(long.df, key = sample, value = count)

clean.counts <-  as.data.frame(tibble::column_to_rownames(tidy.wide, var = "id"))

#order the days in the chronological order
clean.counts <- clean.counts[c("day0", "day1", "day2", "day3", "day4", "day5", "day6",
                               "day8", "day15", "day30", "day60")]

########################
##Get lincRNAs for plot
########################

#get subset desired 
lincRNACounts <- clean.counts[rownames(clean.counts) %in% Tx.lincRNA$gene_id, ]


lincRNACounts <- tibble::rownames_to_column(as.data.frame(lincRNACounts), "id")
#turn into long format by sample and count
long_lincRNACounts <- lincRNACounts %>% tidyr::gather(sample, count, 2:ncol(lincRNACounts))

#look this up 
#this 20 counts was already filtered in import_transcripts.R when creating the dds obj.

#use filter of 20 counts i.e. log2(4.3)
long_lincRNACounts <- long_lincRNACounts[long_lincRNACounts$count >= 4.3, ]

#for the plot turn this into a table 
lincRNAdf <- setNames(as.data.frame(table(long_lincRNACounts$sample)), c("Var1", "Freq_lincRNA"))

#correct order of samples again
sample_names <- c("day0", "day1", "day2", "day3", "day4", "day5", "day6",
                  "day8", "day15", "day30", "day60")

#re-level the factor before using ggplot
lincRNAdf$Var1 <- factor(lincRNAdf$Var1, levels = sample_names)


#plot 
PointPlot <- ggplot(lincRNAdf, aes(x = Var1, color = Var1)) + 
            geom_point(aes(y = Freq_lincRNA), shape = 17, size = 5) +
            xlab(" ") + ylab("Number of lincRNA Expressed") +
  
            #######
            #legend
            ########
            theme(legend.title = element_blank()) +
            theme(legend.position = "none") +
            theme(legend.text = element_text(colour="black", size = 16, face = "plain")) +
            theme( axis.title.x = element_text(family="sans",size = 14, face="bold", hjust=0.5, vjust=-0.5),
                   axis.title.y = element_text(family="sans",size = 14, angle=90, face="bold", hjust=0.5, vjust=1)) +
            theme( axis.text.x = element_text(family = "sans",size = 14, angle=45, face='plain', colour="#353535",   hjust=1, vjust=1) ) +
            theme( axis.text.y = element_text(family = "sans",size = 14, face='plain', colour="#353535",  vjust=0.5) ) +
            theme(axis.line.x = element_line(color="black", size = 0.5),
                  axis.line.y = element_line(color="black", size = 0.5)) +
            theme(legend.background = element_rect()) + 
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),  panel.background = element_blank()) +
            #scale_y_continuous(sec.axis = sec_axis(~.*1.8 , name = "Number of TF Expressed")) +
            #scale_y_continuous(sec.axis = sec_axis(~., name = "Number of TF Expressed")) +
            scale_shape_manual(values= c("day0", "day1", "day2", "day3", "day4", "day5", "day6",
                                         "day8", "day15", "day30", "day60" = 2))
          
#this is for custom colors and custom legend shape because both shapes use the same legend. 
#this is need specific order for the x axis
PointPlot + scale_color_manual(values=c("#f9bdbd", "#f74a4a", "#ffaa56", "#f9e67a", "#a6dd6a", "#05aa60", 
                                            "#b0e5df",  "#4bb7e5","#648cea", "#caa9e5","#9a25d1" )) +
            guides(colour = guide_legend(override.aes = list(shape = 16))) +
            #adjust the x label to reflect more concise option
            scale_x_discrete(labels = c("day0" = "D0", "day1" = "D1", "day2" = "D2", 
                                        "day3" = "D3", "day4" = "D4", "day5" = "D5", 
                                        "day6" = "D6", "day8" = "D8", "day15" = "D15", 
                                        "day30" = "D30", "day60" = "D60"))
