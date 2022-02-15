library(tximport)
library(readr)
library(GenomicFeatures)
library(DRIMSeq)
library(DEXSeq)
# A disadvantage over the exon-level analysis is that we must know in advance all of the possible isoforms 
# that can be generated from a gene locus, all of which are assumed to be contained in the annotation files 
# (FASTA and GTF).
library(stageR)
library(ggplot2)
library(AnnotationDbi)
library(ensembldb)


##############
#IMPORT COUNTS
##############
#workflow for differential transcript usage, there is a tutorial and a paper. 
#https://f1000research.com/articles/7-952
#https://f1000researchdata.s3.amazonaws.com/manuscripts/17982/6c31429f-792e-409b-91cf-3b596f165af8_15398_-_michael_love_v3.pdf?doi=10.12688/f1000research.15398.3&numberOfBrowsableCollections=21&numberOfBrowsableInstitutionalCollections=5&numberOfBrowsableGateways=24

#gene level data can be found here
dir <- paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/quant_data/gencode10.1_duplicatesDir")

list.files(dir)

#read in design information 
samples <- read.table(paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/metadata_file.csv"), sep = ",", header = TRUE)
samples

colnames(samples) <- c("sample_id", "type", "condition", "replicate")

#change levels of the factors 
samples$condition <-factor(samples$condition, levels=c("day0", "day1", "day2", "day3", "day4", "day5", "day6",
                                           "day8", "day15", "day30", "day60"))


#get full path of the quant files 
files <- file.path(dir, samples$sample_id, "quant.sf")
names(files) <- samples$sample_id

all(file.exists(files))
############
#TXIMPORT
############

#By using scaledTPM counts, the estimated proportions fit by DRIMSeq, which are generated from counts, 
#will be equivalent to proportions of the abundance of the isoforms.

txi <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = "scaledTPM") #, tx2gene = txdf) #F gives gene-level summarization 


cts <- txi$counts
cts <- round(cts, digits = 2)
filter <- rowSums(cts >= 50) >= 1 #at least 10 reads in more than 1 sample

#see mikes answer on mcols(obj)$fulBetaConv
cts <- cts[filter,]

load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/transcript_results/cts50_oneSample.RData"))

################################
#make TXDB object from gtf file
################################
#after creating the TXDB object dont need to do it everytime, just load file 
#gtf <- paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/gencode/gencode.v28.primary_assembly.annotation.gtf")
txdb.filename <- (paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/gencode/gencode.v28.annotation.sqlite"))
#txdb <- makeTxDbFromGFF(gtf)
#saveDb(txdb, txdb.filename)


##Once the TxDb database has been generated and saved, it can be quickly reloaded:
txdb <- loadDb(txdb.filename)
txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID") #"EXONCHROM", "EXONEND", "EXONSTART"
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))] #ntx number of transcripts 

#############################
#STATISTICAL ANALYSIS FOR DTU
#############################

range(colSums(cts)/1e6)

all(rownames(cts) %in% txdf$TXNAME)

txdf <- txdf[match(rownames(cts),txdf$TXNAME),]
all(rownames(cts) == txdf$TXNAME)


counts <- data.frame(gene_id=txdf$GENEID,
                     feature_id = txdf$TXNAME,
                     cts)

d <- dmDSdata(counts = counts, samples = samples)

d

methods(class=class(d))
counts(d[1,])[,1:10]

#n = number of samples, but perhaps we want to allow for a few samples to have very low counts 
#min_samps_geen_exp = factor multiple < 1 of n if we think n is too high  
#n.small = transcript does not make up more than 10% of the gene’s expression for at least n.small samples, 
#it will be removed, so can remove this feature or lower the filter 


n = 16 #number of a half of the samples 
n.small = 10
d <- dmFilter(d, min_samps_feature_expr = n.small, min_feature_expr = 10,
              min_samps_feature_prop = n.small, min_feature_prop = 0.1,
              min_samps_gene_exp = n, min_gene_exp = 10)
d

table(table(counts(d)$gene_id))

design_full <- model.matrix(~replicate + condition, data=DRIMSeq::samples(d))
colnames(design_full)


########
#DEXseq
########
#DEXSeq with the NB distribution only can capture correlations among transcript counts within 
#a gene when the DTU is across sample groups defined by covariates in the design matrix. 
#DEXSeq may do a better job modeling such heterogeneity in the biological variability of transcript 
#expression, as it estimates a separate dispersion parameter for each transcript.

# DEXSeq will test – after accounting for total gene expression for this sample and for the proportion 
#of this transcript relative to the others – whether there is a condition-specific difference in the 
#transcript proportion relative to the others.

count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
rownames(samples) <- samples$sample_id
dxd <- DEXSeqDataSet(countData=count.data,
                     sampleData=samples,
                     design=~sample + exon + replicate:exon + condition:exon,#exon and transcript at the same thing
                     featureID=counts(d)$feature_id,
                     groupID=counts(d)$gene_id)

colData(dxd)
head( counts(dxd), 5 ) 
split( seq_len(ncol(dxd)), colData(dxd)$exon ) #first cols correspond to the number of reads the last to the sum of the counts mapped to the rest of the exons from the same gene 
head( featureCounts(dxd), 5 )

#In both cases, the rows are labelled with gene IDs, followed by a colon and the counting bin number.
#(As a counting bin corresponds to an exon or part of an exon, this ID is called the feature ID or exon 
#ID within DEXSeq.

sampleAnnotation(dxd)

system.time({
  dxd <- estimateSizeFactors(dxd) #NORMALIZATION
  dxd <- estimateDispersions(dxd, quiet=FALSE)
  #estimateDispersions is doing..Briefly, per-exon dispersions are calculated using
  # a Cox-Reid adjusted profile likelihood estimation, then a dispersion-mean relation is fitted
  # to this individual dispersion values and finally, the fitted values are taken as a prior in order
  # to shrink the per-exon estimates towards the fitted values.
  dxd <- testForDEU(dxd, reducedModel=~sample + exon + replicate:exon) #we add the replicate data term here
  #so that any dependence of exon usage on replicate will be absorbed bec we are interested in day/condition changes
})


#save(dxd, file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/transcript_results/dxd20_LoweredFilterN16.RData"))

#save time and load the dxd object 
load(file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/transcript_results/dxd50counts.RData"))
plotDispEsts( dxd ) #plot dispersions will show the fitted mean-dispersion values 

#this is the results table.
dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)

table( dxr$padj < 0.1 )
table (tapply(dxr$padj < 0.1, dxr$groupID, any ) ) #We may also ask how many genes are affected

#plotMA(dxr, cex=0.8 )

#compute a per-gene adjusted p-value with uses multiple testing with a gene to a single p-value 
qval <- perGeneQValue(dxr)
dxr.g <- data.frame(gene=names(qval),qval)

columns <- c("featureID","groupID","pvalue")
dxr <- as.data.frame(dxr[,columns])
head(dxr)

# head( geneIDs(dxd) )
# head( exonIDs(dxd) )

#this final table below represent only genes that passed the filter are included in the table, 
#so the table already represents screened genes

##########
#stageR
##########
strp <- function(x) substr(x,1,15) #function to get rid of version numbers 
pConfirmation <- matrix(dxr$pvalue,ncol=1) #pvalue matrix 
dimnames(pConfirmation) <- list(strp(dxr$featureID),"transcript") #add names to the pvalues
pScreen <- qval
names(pScreen) <- strp(names(pScreen))
tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i]) #df with transcript and then gene identifiers 

# stageR class and requires a (preferably named) vector of p-values for the screening hypothesis pScreen 
# and a (preferably named) matrix of p-values for the confirmation stage pConfirmation with columns 
# corresponding to the different contrasts of interest. Note that the rows in  pConfirmation correspond 
# to features (genes) and the features should be identically sorted in  pScreen and pConfirmation.

#Finally, the pScreenAdjusted argument specifies whether the screening p-values have already been 
#adjusted according to FDR control.

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=TRUE, tx2gene=tx2gene) 

#The function stageWiseAdjustment then adjusts the p-values according to a stage-wise analysis.
#The alpha argument specifies the target OFDR level that is used for controlling the fraction of 
#false positive genes across all rejected genes over the entire stage-wise testing procedure.

#The adjusted p-values for genes that did not pass the screening stage are by default set to NA
#method = "dtu" indicates the adapted Holm-Shaffer FWER correction that was specifically tailored for DTU analysis
#This method is used from the thirs transcript and onwards, takes the two most significant p-values and tests
#the adjusted sig level from the screening stage / by the # of transcripts -2 for each gene. 
#The method will return NA p-values for genes with only **one** transcript if the stage-wise testing method equals "dtu"

stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05) #overall False Discovery Rate (OFDR)

#accessory functions...
sigGenes <- (getSignificantGenes(stageRObj))

sigTx <- (getSignificantTx(stageRObj))


#the stage-wise adj p-values are returned using the the function below. These values are adjusted according to the 
#BH FDR criterion and the confirmation hypothesis were adjusted for FWER

#these values can be used directly 

##NOTE: these are the values that can be used only at alpha of 0.05! ONLY. If want to find more sig, need to run again. 
#Only genes that passed the filter are included in the table, so the table already represents screened 
#genes. The transcripts with values in the column, transcript, less than 0.05 pass the **confirmation stage**
#on a target 5% overall false discovery rate, or OFDR
suppressWarnings({
  dex.padj <- getAdjustedPValues(stageRObj, order=TRUE,
                                 onlySignificantGenes=TRUE) #maybe need to change to FALSE to get all genes
})


head(dex.padj)

dex.padj$geneID <- as.character(dex.padj$geneID)
dex.padj$txID <- as.character(dex.padj$txID)


#the transcipts need to overcome 0.05 also
sigDex.padj <- subset(dex.padj, dex.padj$transcript <= 0.05)
#save(dex.padj, file = paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/transcript_results/dex50_padjLoweredFilterN16.RData"))
