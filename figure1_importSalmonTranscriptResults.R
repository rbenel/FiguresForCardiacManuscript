#load libraries 
library(tximport)
library(ensembldb)
library(AnnotationHub)
library("EnsDb.Hsapiens.v92") #the ensembl package in this analysis 
library(readr)
library(rjson)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggfortify)
library(RColorBrewer)

#####################
#Read in Salmon Files
#####################
#get directory path of where quant files are
dir <- paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/quant_data/salmon10.1_quants_bootstrap")
list.files(dir)

#read in design information 
samples <- read.table(paste0("/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/metadata_file.csv"), sep = ",", header = TRUE)
rownames(samples) <- samples$sample_id
samples

#the factors are in alphabetical order, if we want to change that do manually now
#must to this at the start of the analysis
samples$day <-factor(samples$day, levels=c("day0", "day1", "day2", "day3", "day4", "day5", "day6",
                                                     "day8", "day15", "day30", "day60"))

#get full path of the quant files 
files <- file.path(dir, samples$sample_id, "quant.sf")
names(files) <- samples$sample_id

#check that files exist 
all(file.exists(files))

#############
#GET DATABASE
#############
#give name to the ensdb database i am using
edb <- EnsDb.Hsapiens.v92

#get transcripts out of the ensdb package in grange obj form 
#A tutorial for the package can be found here. 
#http://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#10_getting_or_building_ensdb_databasespackages

#Note that we are changing here the return.type to DataFrame, so the method will return 
#a DataFrame with the results instead of the default GRanges.

#listTables(edb)
Tx <- transcripts(edb,
                  columns = c("tx_id", "gene_id", "gene_name", "tx_biotype", "tx_support_level"),
                  return.type = "DataFrame")

#if we want to remove duplicate enteries for gene level analysis
#this will take the *first* entry for each gene 
Tx.genes = Tx[!duplicated(Tx$gene_id),]

#######################
#IMPORT SALMON RESULTS  
#######################
#use tximport function to import the files and use "ignoreTxVersion" to ignore version numbers 
#tutorial can be found here: https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
txi <- tximport(files, type = "salmon", tx2gene = Tx, ignoreTxVersion = TRUE) #F gives gene-level summarization 

#The function takes the transcript level information and summarizes to the gene-level. The “length” matrix can be used to 
#generate an offset matrix for downstream gene-level differential analysis of count matrices

#this is then used as input for deseq which corrects for changes to the average transcript length across samples. 
#When using tximport combined with this DESeqDataSetFromTximport function from Deseq


names(txi)
head(txi$counts)

#####################
#CREATE DESeq OBJECT  
#####################
#Remember to put the factor of highest interest last in the design equation for the design model
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ replicate + day)

nrow(ddsTxi)
#remove rows with sum of less than 20
ddsTxi <- ddsTxi[ rowSums(counts(ddsTxi)) > 20, ] 

nrow(ddsTxi)

#save(ddsTxi, file = paste0(local_path, "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/ddsTxi_obj.RData"))

####This is only for use if one needs the annotations with the counts bec they are found in this file
#vsd <- vst(ddsTxi, blind = FALSE)
#mat <- assay(vsd)
#mat2 <- limma::removeBatchEffect(mat, vsd$replicate)

#ddsTxi_counts_annot <- merge(mat2, Tx.genes, by.x = "row.names", by.y = "gene_id")
#save(mat2, ddsTxi_counts_annot, file = paste0(local_path, "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/data/deseq2/ddsCounts_annot_vsd.RData"))
