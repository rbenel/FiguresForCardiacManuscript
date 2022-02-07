
# FiguresForManuscript

This repo will contain all scripts used to produce the figures found in the manuscript

# Figure 1

## PCA plot: 
Step1: `figure1_importSalmonTranscriptResults.R` 
Step2: `figure1_exploratory_visualization_PCA_heatmap.R`

## boxplot lincRNA genes V. PCG:
Step1:figure1_boxplot_comparisonPCG_lncRNA.R

## Violin plot lincRNA gene expression:
Step1:`figure1_violinplot_GloballincRNAexpression.R`

## Point plot sum "daily" expression:
Step1:`figure1_PointPlotSum_lincRNA.R`

## Volcano plot lincRNA day0 v. day60:
step1:`figure1_VolcanoPlot_lincRNAlabeled.R`

## Heatmap lincRNA expression:
Step1:`figure1_importSalmonTranscriptResults.R` 
Step2:`figure1_LRT_timeseries_deseq2.R`
Step3:`figure1_normalizedCounts_heatmap.R`

# Figure 2

## Heatmap of all DE lincRNA
Step 1: figure2_SpearmanCorr_DElincRNA.R

## Network of DE lincRNA with correlation score of 0.9 with one another 
Step 1: figure2_networkFor_DElincRNA.R

## Topological Attributes for the DE lincRNA Network 
Step 1: Functions_TopolgicalFeatures.R -- list of functions needed for Step #2
Step 2: figure2_NetworkTopologicalCharacteristics_lincRNA_lincRNA.R
