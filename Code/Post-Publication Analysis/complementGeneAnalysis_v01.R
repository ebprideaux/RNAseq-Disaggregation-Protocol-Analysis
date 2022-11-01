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

analysisGenes_vector = c("C2", "C3", "C4A", "C4B", "C5", "C6")

complementGenes_df = countsforanalysis_df[analysisGenes_vector,]

# define sample groups
intactSample_vector = grep("Intac", colnames(countsforanalysis_df), value = TRUE)
liberSample_vector = grep("Liber", colnames(countsforanalysis_df), value = TRUE)
libfSample_vector = grep("Lib.F", colnames(countsforanalysis_df), value = TRUE)
subaSample_vector = grep("Sub.A", colnames(countsforanalysis_df), value = TRUE)

# create sample group DFs
cg_intact_df = complementGenes_df[,intactSample_vector]
cg_liber_df = complementGenes_df[,liberSample_vector]
cg_libf_df = complementGenes_df[,libfSample_vector]
cg_suba_df = complementGenes_df[,subaSample_vector]

## visualize with bar chart
setwd("/Users/bartonprideaux/Library/CloudStorage/OneDrive-Personal/Documents/Education/Graduate_School/UCSD/Research/Projects/RA-Disagg/Analysis/PostProject")

# by treatment
gene_vector=rep(rownames(complementGenes_df), 4)
treatment_vector=c(rep("Intact", length(rownames(complementGenes_df))), 
                   rep("Liberase", length(rownames(complementGenes_df))),
                   rep("LibFlavo", length(rownames(complementGenes_df))),
                   rep("SubA", length(rownames(complementGenes_df))))

value_vector=c(rowMeans(cg_intact_df), rowMeans(cg_liber_df), rowMeans(cg_libf_df), rowMeans(cg_suba_df))
barChart_df=data.frame(gene_vector, treatment_vector, value_vector)

png("compGenes_byTreatment_bc.png", width = 1250, height = 750, units = "px")
ggplot(barChart_df, aes(fill=treatment_vector, y=value_vector, x=gene_vector)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle=90, hjust=1)) +
  xlab("Gene") + ylab("Gene Expression Value (geTMM)") + guides(fill=guide_legend(title="Treatment"))
dev.off()