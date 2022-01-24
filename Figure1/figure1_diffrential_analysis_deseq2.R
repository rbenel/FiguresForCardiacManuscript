# Note that the two transformations offered by DESeq2 are provided for applications other than 
# differential testing. For differential testing we recommend the DESeq function applied to raw counts

#see tutorial here for time-series RNAseq analysis 
#https://www.bioconductor.org/help/course-materials/2015/LearnBioconductorFeb2015/B02.1.1_RNASeqLab.html
#RNAseq workflow 
#http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

library(DESeq2)
library(apeglm)
library(ashr)
library(Rmosek)
library(REBayes)
library(vsn)
library(pheatmap)
library("RColorBrewer")
library(dplyr)
library(DEGreport)
library(tidyr)
library(tibble)
library(limma)
library("ggbeeswarm")
library(ggplot2)

local_path <- getwd()

#this is the ensembl version 
load(paste0( "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/ddsTxi_obj.RData"))
#this is the GENCODE version 
#load(paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/GENCODEddsTxi_obj.RData"))
#TF annotations
TFs_annot = read.csv(paste0( "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/Practice 18.12.16/data/tf_human_annotations2018.txt"),
                     header = TRUE, sep = ",", stringsAsFactors = F)

#TF_interactome <- read.csv(paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/Practice 18.12.16/data/TF_listInteractomes.csv"), header = T)

#updated March 24th 2020
#updated lnc lincRNA and funcLincRNA annotations
load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/Practice 18.12.16/data/Tx.linc_lncRNA_FunctlncRNA_release92TSL.RData"))


##############
####WALD TEST
##############

ddsTxi_wald <- DESeq(ddsTxi)

counts <- as.data.frame(counts(ddsTxi_wald, normalized = TRUE))
counts$rownames <- rownames(counts)

resultsNames(ddsTxi_wald)

resWald <- results(ddsTxi_wald, name="day_day60_vs_day0")

mcols(resWald, use.names = TRUE)
#baseMean, is a just the average of the normalized count values, divided by the size factors, 
#taken over all samples in the DESeqDataSet.
#lfcSE, the standard error estimate for the log2 fold change estimate.

# The purpose of a test for differential expression is to test whether the data provides sufficient 
# evidence to conclude that this value is really different from zero. DESeq2 performs for each gene a 
# hypothesis test to see whether evidence is sufficient to decide against the null hypothesis that there 
# is zero effect of the treatment on the gene and that the observed difference between treatment and control 
# was merely caused by experimental variability

summary(resWald)
#as we can see the FDR rate is 0.1, if we want to be stricter we should re-reun res 

res.05 <- results(ddsTxi_wald, alpha = 0.05)
table(res.05$padj < 0.05)

#or we can raise the log2FC threshold lfcthreshold = 1 means must double in counts
resLFC1 <- results(ddsTxi_wald, lfcThreshold=1)
table(resLFC1$padj < 0.1)

# we can't rely on pvalues, need padj
sum(resWald$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(resWald$pvalue))

#by definition the pvalye would have up to 5% below 0.05, so we use p.adj
# this method calculates for each gene an adjusted p value that answers the following question: 
# if one called significant all genes with an adjusted p value less than or equal to this gene’s adjusted 
# p value threshold, what would be the fraction of false positives (the false discovery rate, FDR) among them

sum(resWald$padj < 0.1, na.rm=TRUE)

resSig <- subset(resWald, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ]) #up-regulation 

####################
##PLOTTING RESULTS
###################

topGene <- rownames(resWald)[which.min(resWald$padj)]
plotCounts(ddsTxi_wald, gene = topGene, intgroup=c("day"))

geneCounts <- plotCounts(ddsTxi_wald, gene = topGene, intgroup = c("day","replicate"),
                         returnData = TRUE)

ggplot(geneCounts, aes(x = day, y = count, color = day, shape = replicate)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

ggplot(geneCounts, aes(x = day, y = count, color = day, shape = replicate, group = replicate)) +
  scale_y_log10() + geom_point(size =3) + geom_line()


#can also plot top genes
# We can furthermore cluster significant genes by their profiles. 
# We extract a matrix of the shrunken log2 fold changes using the coef function:
betas <- coef(ddsTxi_wald)
colnames(betas)

topGenes <- head(order(resWald$padj),100)

mat <- betas[topGenes, c(4:13)] #1:3 are the names and the intercept and in this case the replicate 1 v. 2. 3

thr <- 5
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr


#colors <- colorRampPalette( rev(brewer.pal(9, "RdBlu")) )(255)

pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101), show_rownames = F, trace="none", dendrogram="row",
         Colv=FALSE, mar=c(12,8))

#############
#MA-Plot
############
# On the y-axis, the “M” stands for “minus” – subtraction of log values is equivalent to the log of the ratio 
# – and on the x-axis, the “A” stands for “average”. You may hear this plot also referred to as a 
# mean-difference plot, or a Bland-Altman plot.

#SHRINKAGE
#the lfc use a normal prior distribution, centered on zero and with a scale 
#that is fit to the data. The shrunken log fold changes are useful for ranking and visualization

#three types of shrinkage 
#"ashr" the adaptive shrinkage estimator from the ashr package (Stephens 2016). Here DESeq2 uses 
#the ashr option to fit a mixture of Normal distributions to form the prior, with method="shrinkage".

#The DESeq2 package uses a Bayesian procedure to moderate (or “shrink”) log2 fold changes from genes 
#with very low counts and highly variable counts, as can be seen by the narrowing of the vertical spread 
#of points on the left side of the MA-plot


resLFC <- lfcShrink(ddsTxi_wald, coef="day_day60_vs_day0", type="ashr") #add threshold if need to be stricter
DESeq2::plotMA(resLFC, cex = 0.8) #, ylim=c(-2,2)) #must specifcy if using loading limma library!

abline(h=c(-1,1), col="dodgerblue", lwd=2)

# one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. 
# One can then recover the gene identifiers by saving the resulting indices:
idx <- identify(resLFC$baseMean, resLFC$log2FoldChange)
rownames(resLFC)[idx]

#view outliers 
par(mar=c(8,5,2,2))
boxplot(log10(assays(ddsTxi_wald)[["cooks"]]), range=0, las=2)

#plot dispersions of dds object 
#Plotting the dispersion estimates is a useful diagnostic
#with the final estimates shrunk from the gene-wise estimates towards the fitted estimates
plotDispEsts(ddsTxi_wald)

#check the hist of the pvalues 
hist(resWald$pvalue[resWald$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

###############################
##RESULTS 
##############################
#make a df of the pairwise contrasts using Wald..
#the problem with this is that there are no pvalues... so we dont know how trustworthy these 
#FCs are... 
FC_df <- data.frame()

days = paste0("day", c("1", "2", "3", "4", "5", "6", "8", "15", "30", "60"))

for (i in days){
  res_df <- results(ddsTxi_wald, contrast=c("day", i,  "day0"), test="Wald")
  res_df <- as.data.frame(res_df)
  ##if the padj is more than 0.05 input 0 in FC column
  #this solves the issue of FC and padj value 
  res_df$log2FoldChange <- ifelse(res_df$padj > 0.05, 0, res_df$log2FoldChange) 
  ##this line is for the FC, if it is less than abs (0.58) we will turn it into a zero
  #but we don't need this if we use another script for norm_FCanalysis
  #res_df$log2FoldChange <- ifelse(abs(res_df$log2FoldChange) < 0.58, 0, res_df$log2FoldChange)
  
  colnames(res_df) <- paste0(colnames(res_df), "_", i)
  
  if(ncol(FC_df) == 0) {
    FC_df = res_df
  } else {
    FC_df = cbind(FC_df, res_df)
    
  }
}

#in order to subset the rows of the diffrentially expressed genes in the df
FC <-c("log2FoldChange")
subset.FC.df <- dplyr::select(FC_df, matches(FC))



# for(i in 1:ncol(subset.FC_df)){
#   for(j in 2:ncol(subset.FC_df)){
#     tmp.df <- subset.FC_df[, 1:2]
#     test <- tmp.df %>% dplyr::filter(tmp.df[,1] >= log2(4) & (tmp.df[,2] < 0.05))
# 
#   }
# }


#write.csv(subset.FC.df, paste0(, "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/padj_OnlyFCresWald.txt"))

#write.csv(FC_df, paste0(, "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/allresWald.txt"))

###############
###VOLCANO PLOT
################
#save the results from the resWald which are 60 v. 0 so we can draw the updated volcano plot in a new script 
#save(resWald, file = "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/resDay0Day60.RData")

#turn into a data frame and check highly sig padj and can see fold change. 
res_df <- as.data.frame(res0X60)

res_df$rownames <- rownames(res_df)

############################################################################
###ALL THIS IS FOR GGPLOT -- NAMES OF TFS
#merge this info for later: so can use external gene names 
res_df.merge <- merge(res_df, TFs_annot, by.x = "rownames", by.y = "ensembl_gene_id", all.x = TRUE)

TFmarkers <- c("POU5F1", "T/Bra", "SOX17", "MESP1", "ISL1", "NKX2-5")

#this allows us to label just markers on our plot 
res_df.merge$markers <- as.factor(ifelse(res_df.merge$external_gene_name %in% TFmarkers, paste(res_df.merge$external_gene_name) , NA))

##########################################################################
#check the distribution of pvalues in diffrential expression 
hist(res_df$pvalue,60)

# Make a basic volcano plot to see all results
plot(res_df$log2FoldChange,-log10(res_df$pvalue), col='grey')
plot(res_df$log2FoldChange,-log10(res_df$pvalue), col='grey', xlim=c(-10,10), ylim=c(0, 50),  xlab = "log2 FoldChange", ylab="-log10 (pvalue)") #change x lim or ylim

# Add colored points: 
#points that overcome our qualifications 
with(subset(res_df, (padj < 0.05 & abs(log2FoldChange) > 0.58) & (rownames %in% TFs_annot$ensembl_gene_id)), points(log2FoldChange, -log10(pvalue), pch=20, col= '#F8766D')) #pink
#with(subset(res_df, (rownames %in% TFs_annot$ensembl_gene_id)), points(log2FoldChange, -log10(pvalue), pch=20, col='#F8766D')) #pink

with(subset(res_df, (padj < 0.05 & abs(log2FoldChange) > 0.58) & (rownames %in% Tx.lincRNA$gene_id)), points(log2FoldChange, -log10(pvalue), pch=20, col= '#00bfc4' )) #blue
#with(subset(res_df, (rownames %in% Tx.lincRNA$gene_id)), points(log2FoldChange, -log10(pvalue), pch=20, col='#00bfc4')) #blue

abline(v= c(-0.58, 0.58), col = c("black"), lty=c(2) , lwd = (3))
#with(subset(res_df, (rownames %in% TF_interactome$ensembl_gene_id)), points(log2FoldChange, -log10(pvalue), pch=20, col='#8c00c4')) #purple

#####ggplot that has labels of the TF!!! 
## need to merge the res_df first 
p <- ggplot(res_df.merge, aes(log2FoldChange, -log10(pvalue)))
p + geom_point(color = "grey") +
  geom_point(data=subset(res_df.merge, abs(log2FoldChange) > 0.58 & padj < 0.05 & (rownames %in% TFs_annot$ensembl_gene_id)), #label only what has external name i.e. the TFs
                        aes(log2FoldChange, -log10(pvalue)), color = "#F8766D") + 
  geom_point(data=subset(res_df.merge, abs(log2FoldChange) > 0.58 & padj < 0.05 & (rownames %in% Tx.lincRNA$gene_id)), #label only what has external name i.e. the TFs
                        aes(log2FoldChange, -log10(pvalue)), color = "#00bfc4") + 
  scale_fill_manual(values = c("grey")) + scale_y_continuous(limits = c(0, 50)) +
  scale_x_continuous(limits = c(-10, 10)) + 
  geom_label(data=subset(res_df.merge, abs(log2FoldChange) > 0.58 & padj < 0.05), #label only what has external name or use "marker" column
            aes(log2FoldChange, -log10(pvalue) ,label=markers), vjust = 0.5, 
            hjust = "right", nudge_x = TRUE, nudge_y = TRUE, fill = "#F8766D") +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size = 14, face = "plain")) +
  theme( axis.title.x = element_text(family="sans",size = 20, face="bold", hjust=0.5, vjust=-0.5),
         axis.title.y = element_text(family="sans",size = 20, angle=90, face="bold", hjust=0.5, vjust=1)) +
  theme( axis.text.x = element_text(family = "sans",size = 14, angle=0, face='plain', colour="#353535",   hjust=1, vjust=1) ) +
  theme( axis.text.y = element_text(family = "sans",size = 14, face='plain', colour="#353535",  vjust=0.5) ) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  theme(legend.background = element_rect()) + 
  theme(legend.position="top") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  panel.background = element_blank())

