###############################################

if(!require("maptools")) install.packages("maptools")
library("maptools")
if(!require(colorRamps)) install.packages("colorRamps")
library('colorRamps')
if(!require(RColorBrewer)) install.packages("RColorBrewer")
library(RColorBrewer)
if(!require(gplots)) install.packages("gplots")
library(gplots)
if(!require(plyr)) install.packages("plyr")
library(plyr)
if(!require(openxlsx)) install.packages("openxlsx")
library("openxlsx")
if(!require(reshape2)) install.packages("reshape2")
library("reshape2")
if(!require(ggplot2)) install.packages("ggplot2")
library("ggplot2")

library("edgeR")
library("limma")
library("biomaRt")
library("GenomicFeatures")

### Analyze percent of cytokines, growth factors and proteases included

# import list of cytokines, growth factors and proteases
cytokine_fp = "/Volumes/Portable Drive/Documents/Graduate/Data/GeneLists/GO_term_summary_cytokines.xlsx"
cytokine_df = as.data.frame(read.xlsx(cytokine_fp, colNames = TRUE))
cytokine_vector = toupper(as.vector(cytokine_df$Symbol))

growthfactor_fp = "/Volumes/Portable Drive/Documents/Graduate/Data/GeneLists/GO_term_summary_growthfactors.xlsx"
growthfactor_df = as.data.frame(read.xlsx(growthfactor_fp,colNames = TRUE))
growthfactor_vector = toupper(as.vector(growthfactor_df$Symbol))

protease_fp = "/Volumes/Portable Drive/Documents/Graduate/Data/GeneLists/Protease_Human.0111819.s004.xlsx"
protease_df = as.data.frame(read.xlsx(protease_fp, colNames = TRUE))
protease_vector = toupper(as.vector(protease_df$GeneSymbol))

# Identify datapoints and trim low expression genes
plot_matrix = tmmcounts_hgnc_pcexmtgenes_normedwithall_df

expr_min = 0.50

plot_matrix[,"trimmed_mean"] = apply(plot_matrix, 1, function(x) mean(x, trim = 0.1))
plot_matrix = plot_matrix[rownames(head(plot_matrix[order(plot_matrix["trimmed_mean"], decreasing = T),], expr_min*nrow(plot_matrix))),]
plot_matrix = subset(plot_matrix, select=-c(trimmed_mean))
dim(plot_matrix)

# check for percentage of cytokines, growth factors and proteases included
cytokine_pctinc = 1-length(setdiff(cytokine_vector, rownames(plot_matrix)))/length(cytokine_vector)
print(cytokine_pctinc)
growthfactor_pctinc = 1-length(setdiff(growthfactor_vector, rownames(plot_matrix)))/length(growthfactor_vector)
print(growthfactor_pctinc)
protease_pctinc = 1-length(setdiff(protease_vector, rownames(plot_matrix)))/length(protease_vector)
print(protease_pctinc)
