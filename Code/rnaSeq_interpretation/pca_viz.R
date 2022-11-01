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
if(!require(irlba)) install.packages("irlba")

library("irlba")
library("sva")
library("edgeR")
library("limma")
library("biomaRt")
library("GenomicFeatures")
library("stringr")
library("edgeR")
library("limma")
library("biomaRt")
library("GenomicFeatures")


### Perform PCA Analysis
# Identify datapoints and trim low expression genes
set.seed(61)
plot_matrix = tmmcounts_hgnc_pcexmtgenes_normedwithall_df
expr_min = 0.50

plot_matrix[,"mean"] = apply(plot_matrix, 1, function(x) mean(x))
plot_matrix = plot_matrix[rownames(head(plot_matrix[order(plot_matrix["mean"], decreasing = T),], expr_min*nrow(plot_matrix))),]
plot_matrix = subset(plot_matrix, select=-c(mean))
dim(plot_matrix)

plot_matrix = apply(plot_matrix, 1:2, function(x) log2(x+1))

# create labels
id = colnames(plot_matrix)
condition = substring(id,1,2)
sample = substring(id,4,7)
treatment = substring(id,9,13)
batch = substring(id,15)

# PCA and extract info
set.seed(58)
plot_PCA = prcomp_irlba(plot_matrix, n=3)
summary(plot_PCA)

a = plot_PCA$rotation[,1]
b = plot_PCA$rotation[,2]
c = plot_PCA$rotation[,3]

var_pc1 = paste(round(summary(plot_PCA)[[7]][2,1], 4) * 100, "%", sep = "")
var_pc2 = paste(round(summary(plot_PCA)[[7]][2,2], 4) * 100, "%", sep = "")
var_pc3 = paste(round(summary(plot_PCA)[[7]][2,3], 4) * 100, "%", sep = "")

xname = paste("PC1(", var_pc1, ")", sep = "")
yname = paste("PC2(", var_pc2, ")", sep = "")
zname = paste("PC3(", var_pc3, ")", sep = "")

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

setwd("/Users/bartonprideaux/OneDrive/Documents/Education/Graduate_School/UCSD/Research/Projects/RA-Disagg/Analysis/Round3/Figures/PCA")

# Plot chart #1
plotColor_driver = condition

uniqueColors_needed = length(unique(factor(plotColor_driver)))
plotColors = sample(col_vector, uniqueColors_needed)
names(plotColors) = levels(factor(plotColor_driver))
colors = unname(plotColors[plotColor_driver])

shapes = rep(19, length(plotColor_driver))

dd = as.character(id)
names(colors) = unique(dd)
names(shapes) = unique(dd)

png("pc1pc2_tmm_log2pseduo_PCexMTgenes-normedwithall_condition.png", width = 600, height = 600, units = "px")
plot(a, b, col = colors[dd], pch = shapes[dd], xlab = xname, ylab = yname, cex = 1.2)
pointLabel(a, b, labels = sample, cex = 1.2)
legend("topright", legend = c("OA", "RA"), col = plotColors, border = "black", bty = "n", pch = shapes, cex = 1.2)
dev.off()

png("pc2pc3_tmm_log2pseduo_PCexMTgenes-normedwithall_condition.png", width = 600, height = 600, units = "px")
plot(b, c, col = colors[dd], pch = shapes[dd], xlab = yname, ylab = zname, cex = 1.2)
pointLabel(b, c, labels = sample, cex = 1.2)
legend("bottomright", legend = c("OA", "RA"), col = plotColors, border = "black", bty = "n", pch = shapes, cex = 1.2)
dev.off()

# Plot chart #2
plotColor_driver = treatment

uniqueColors_needed = length(unique(factor(plotColor_driver)))
plotColors = sample(col_vector, uniqueColors_needed)
names(plotColors) = levels(factor(plotColor_driver))
colors = unname(plotColors[plotColor_driver])

shapes = rep(19, length(plotColor_driver))

dd = as.character(id)
names(colors) = unique(dd)
names(shapes) = unique(dd)

png("pc1pc2_tmm_log2pseduo_PCexMTgenes-normedwithall_treatment.png", width = 600, height = 600, units = "px")
plot(a, b, col = colors[dd], pch = shapes[dd], xlab = xname, ylab = yname, cex = 1.2)
pointLabel(a, b, labels = sample, cex = 1.2)
legend("topright", legend = c("Intact", "Lib Flavo", "Liberase", "RNA Stat", "Sub A"), col = plotColors, border = "black", bty = "n", pch = shapes, cex = 1.2)
dev.off()

png("pc2pc3_tmm_log2pseduo_PCexMTgenes-normedwithall_treatment.png", width = 600, height = 600, units = "px")
plot(b, c, col = colors[dd], pch = shapes[dd], xlab = yname, ylab = zname, cex = 1.2)
pointLabel(b, c, labels = sample, cex = 1.2)
legend("bottomright", legend = c("Intact", "Lib Flavo", "Liberase", "RNA Stat", "Sub A"), col = plotColors, border = "black", bty = "n", pch = shapes, cex = 1.2)
dev.off()

# Plot chart #3
plotColor_driver = batch

uniqueColors_needed = length(unique(factor(plotColor_driver)))
plotColors = sample(col_vector, uniqueColors_needed)
names(plotColors) = levels(factor(plotColor_driver))
colors = unname(plotColors[plotColor_driver])

shapes = rep(19, length(plotColor_driver))

dd = as.character(id)
names(colors) = unique(dd)
names(shapes) = unique(dd)

png("pc1pc2_tmm_log2pseduo_PCexMTgenes-normedwithall_batch.png", width = 600, height = 600, units = "px")
plot(a, b, col = colors[dd], pch = shapes[dd], xlab = xname, ylab = yname, cex = 1.2)
pointLabel(a, b, labels = sample, cex = 1.2)
legend("topright", legend = c("Batch 1", "Batch 2", "Batch 3"), col = plotColors, border = "black", bty = "n", pch = shapes, cex = 1.2)
dev.off()

png("pc2pc3_tmm_log2pseduo_PCexMTgenes-normedwithall_batch.png", width = 600, height = 600, units = "px")
plot(b, c, col = colors[dd], pch = shapes[dd], xlab = yname, ylab = zname, cex = 1.2)
pointLabel(b, c, labels = sample, cex = 1.2)
legend("bottomright", legend = c("Batch 1", "Batch 2", "Batch 3"), col = plotColors, border = "black", bty = "n", pch = shapes, cex = 1.2)
dev.off()