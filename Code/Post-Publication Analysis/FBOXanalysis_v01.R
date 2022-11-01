###############################################
library("biomaRt")
library("GenomicFeatures")
library("edgeR")
library("limma")
if(!require("plyr")) install.packages("plyr")
library("plyr")
if(!require("openxlsx")) install.packages("openxlsx")
library("openxlsx")
if(!require("qvalue")) install.packages("qvalue")
library("qvalue")
library("EnhancedVolcano")
library(gplots)
library("RColorBrewer")
library(ggplot2)

### FBOX family expression analysis - Disagg Data ###
# import and select genes for analysis
countsforanalysis_df = tmmcounts_hgnc_pcexmtgenes_normedwithall_df

analysisGenes_vector = c("FBXO2", "FBXO6",  "FBXO10",  "FBXO16", "FBXO17", "FBH1", "FBXL19", "FBXO27", 
                         "FBXO33", "FBXO39", "FBXO41", "FBXO42","FBXO44", "FBXO45",
                         "SKP2", "FBXL3", "FBXL6", "FBXL7", "FBXL8", "KDM2B", "KDM2A", "FBXL12", "FBXL14", 
                         "FBXL16", "FBXL17", "FBXL22",
                         "BTRC", "FBXW2", "FBXW5", "FBXW7", "FBXW10", "FBXW11")

# create dataframe of selected genes
FBOXanalysis_df = countsforanalysis_df[analysisGenes_vector,
                                       grep("Intac", colnames(countsforanalysis_df), value = TRUE)]

RAsample_vector = grep("RA", colnames(FBOXanalysis_df), value = TRUE)
OAsample_vector = grep("OA", colnames(FBOXanalysis_df), value = TRUE)

## visualize with heatmap
setwd("/Users/bartonprideaux/Library/CloudStorage/OneDrive-Personal/Documents/Education/Graduate_School/UCSD/Research/Projects/Miscellaneous/2022.09_FBOXanalysis/heatmaps")

# by disease
new_colOrder=c(RAsample_vector, OAsample_vector)

heatmap_df = FBOXanalysis_df[,new_colOrder]

color_gradient=7
hm_cols = colorRampPalette(rev(brewer.pal(color_gradient, "RdYlBu")), bias=0.6)

png("RAvsOA_FBOX_disaggIntac_hm.png", width = 1000, height = 750, units = "px")
heatmap.2(as.matrix(heatmap_df), 
          Rowv=F, Colv=F, dendrogram = "none",
          trace="none", scale = "row", density.info='histogram', denscol="black",
          labRow = rownames(heatmap_df), col=hm_cols, 
          symkey=F, 
          cexRow=1.0, cexCol=1.0, 
          margins=c(15,10)
)
dev.off()

## visualize with bar chart
setwd("/Users/bartonprideaux/Library/CloudStorage/OneDrive-Personal/Documents/Education/Graduate_School/UCSD/Research/Projects/Miscellaneous/2022.09_FBOXanalysis/barCharts")

# by region
gene_vector=rep(rownames(FBOXanalysis_df), 2)
disease_vector=c(rep("RA", 32), rep("OA", 32))
value_vector=c(rowSums(FBOXanalysis_df[,RAsample_vector]), rowSums(FBOXanalysis_df[,OAsample_vector]))
barChart_df=data.frame(gene_vector, disease_vector, value_vector)

png("RAvsOA_FBOX_disaggIntac_bc.png", width = 1250, height = 750, units = "px")
ggplot(barChart_df, aes(fill=disease_vector, y=value_vector, x=gene_vector)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle=90, hjust=1)) +
  xlab("Gene") + ylab("Gene Expression Value (geTMM)") + guides(fill=guide_legend(title="Disease"))
dev.off()




### FBOX family expression analysis - FLS DATA ###
# import and select genes for analysis
flsWorkbook_fp = "/Users/bartonprideaux/Library/CloudStorage/OneDrive-Personal/Documents/Education/Graduate_School/UCSD/Research/Projects/Miscellaneous/2021.04_RNAseqAnalysis/Data/Output/RNAseq_RAOAFLS_geTMMpcexmtgenes.xlsx"
countsforanalysis_df = as.data.frame(read.xlsx(flsWorkbook_fp, colNames = TRUE, rowNames = TRUE))

analysisGenes_vector = c("FBXO2", "FBXO6",  "FBXO10",  "FBXO16", "FBXO17", "FBH1", "FBXL19", "FBXO27", 
                         "FBXO33", "FBXO39", "FBXO41", "FBXO42","FBXO44", "FBXO45",
                         "SKP2", "FBXL3", "FBXL6", "FBXL7", "FBXL8", "KDM2B", "KDM2A", "FBXL12", "FBXL14", 
                         "FBXL16", "FBXL17", "FBXL22",
                         "BTRC", "FBXW2", "FBXW5", "FBXW7", "FBXW10", "FBXW11")

# create dataframe of selected genes
FBOXanalysis_df = countsforanalysis_df[analysisGenes_vector,]

RAsample_vector = grep("RA", colnames(FBOXanalysis_df), value = TRUE)
OAsample_vector = grep("OA", colnames(FBOXanalysis_df), value = TRUE)

## visualize with heatmap
setwd("/Users/bartonprideaux/Library/CloudStorage/OneDrive-Personal/Documents/Education/Graduate_School/UCSD/Research/Projects/Miscellaneous/2022.09_FBOXanalysis/heatmaps")

# by disease
new_colOrder=c(RAsample_vector, OAsample_vector)

heatmap_df = FBOXanalysis_df[,new_colOrder]

color_gradient=7
hm_cols = colorRampPalette(rev(brewer.pal(color_gradient, "RdYlBu")), bias=0.6)

png("RAvsOA_FBOX_FLS_hm.png", width = 1000, height = 750, units = "px")
heatmap.2(as.matrix(heatmap_df), 
          Rowv=F, Colv=F, dendrogram = "none",
          trace="none", scale = "row", density.info='histogram', denscol="black",
          labRow = rownames(heatmap_df), col=hm_cols, 
          symkey=F, 
          cexRow=1.0, cexCol=1.0, 
          margins=c(15,10)
)
dev.off()

## visualize with bar chart
setwd("/Users/bartonprideaux/Library/CloudStorage/OneDrive-Personal/Documents/Education/Graduate_School/UCSD/Research/Projects/Miscellaneous/2022.09_FBOXanalysis/barCharts")

# by region
gene_vector=rep(rownames(FBOXanalysis_df), 2)
disease_vector=c(rep("RA", 32), rep("OA", 32))
value_vector=c(rowSums(FBOXanalysis_df[,RAsample_vector]), rowSums(FBOXanalysis_df[,OAsample_vector]))
barChart_df=data.frame(gene_vector, disease_vector, value_vector)

png("RAvsOA_FBOX_FLS_bc.png", width = 1250, height = 750, units = "px")
ggplot(barChart_df, aes(fill=disease_vector, y=value_vector, x=gene_vector)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle=90, hjust=1)) +
  xlab("Gene") + ylab("Gene Expression Value (geTMM)") + guides(fill=guide_legend(title="Disease"))
dev.off()

