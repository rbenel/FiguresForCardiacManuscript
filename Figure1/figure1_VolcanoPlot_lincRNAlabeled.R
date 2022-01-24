library(ggplot2)
library(dplyr)
library(ggrepel)

###########
#Load data
##########
#results from DESeq2 for Day 0 v Day 60 -- Wald
load(file = "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/resDay0Day60.RData")

#updated lnc lincRNA and funcLincRNA annotations
load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/Practice 18.12.16/data/Tx.linc_lncRNA_FunctlncRNA_release92TSL.RData"))

###############
#Results of DE
################
#turn into a data frame and check highly sig padj and can see fold change. 
res_df <- as.data.frame(resWald)

res_df$rownames <- rownames(res_df)

#######################
#Add names of lincRNAs
#######################
#merge this info now for the plot: so can use external gene names 
#res_df.merge <- merge(res_df, TFs_annot, by.x = "rownames", by.y = "ensembl_gene_id", all.x = TRUE)
res_df.merge <- merge(res_df, Tx.lincRNA, by.x = "rownames", by.y = "gene_id", all.x = TRUE)

#write output 
#write.csv(res_df.merge, file = "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Diff_exp/Wald_results/DEres_funclncRNALabels.csv",
#          row.names = F)

#lincRNAs we are interested in labeleing 
lincRNAmarkers <- c("TSIX", "TERC", "FENDRR", "RMST", "MALAT1", "CARMN", 
                    "NEAT1", "DANCR", "MIAT", "PVT1")

#this allows us to label just markers on our plot 
#res_df.merge$markers <- as.factor(ifelse(res_df.merge$external_gene_name %in% TFmarkers, paste(res_df.merge$external_gene_name) , NA))

##############
#Volcano Plot
##############
## need to merge the res_df first 
VolcanoPlot <- ggplot(res_df.merge, aes(log2FoldChange, -log10(pvalue))) + 
        geom_point(color = "grey") +
        geom_point(data=subset(res_df.merge, abs(log2FoldChange) > 0.58 & padj < 0.05 & (rownames %in% Tx.lincRNA$gene_id)), #label only what has gene id i.e. the lincRNAs
                   aes(log2FoldChange, -log10(pvalue)), color = "#00bfc4") + 
        scale_fill_manual(values = c("grey")) + 
        scale_y_continuous(limits = c(0, 65)) +
        scale_x_continuous(limits = c(-10, 10)) + 
        ############
        #Label Plot
        ############
        #use this instead of there are many labels 
        geom_label(data=subset(res_df.merge, gene_name %in% lincRNAmarkers), #label only what has external name or use "marker" column
                   aes(log2FoldChange, -log10(pvalue) ,label=gene_name), vjust = 0.5,
                   hjust = "right", nudge_x = TRUE, nudge_y = TRUE, fill = "#00bfc4", fontface = 'bold') +
        #can use this line if there are few lables
        # geom_label_repel(data=subset(res_df.merge, gene_name %in% lincRNAmarkers), #label only what has gene id or abs(log2FoldChange) > 5 & padj < 0.00001 &
        #            aes(log2FoldChange, -log10(pvalue), label = gene_name), 
        #            box.padding = unit(0.1, "lines"), force = 2, fill = "#00bfc4",
        #                              segment.color = NA,  fontface = 'bold') + 
        #######
        #legend
        #######
        theme(legend.title = element_blank()) +
        theme(legend.text = element_text(colour="black", size = 14, face = "plain")) +
        theme( axis.title.x = element_text(family="sans",size = 18, face="bold", hjust=0.5, vjust=1),
               axis.title.y = element_text(family="sans",size = 18, angle=90, face="bold", hjust=0.5, vjust=1)) +
        theme( axis.text.x = element_text(family = "sans",size = 14, angle=0, face='plain', colour="#353535",   hjust=1, vjust=1) ) +
        theme( axis.text.y = element_text(family = "sans",size = 14, face='plain', colour="#353535",  vjust=0.5) ) +
        theme(axis.line.x = element_line(color="black", size = 0.5),
              axis.line.y = element_line(color="black", size = 0.5)) +
        theme(legend.background = element_rect()) + 
        theme(legend.position="top") +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),  panel.background = element_blank())
      
VolcanoPlot
