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
if(!require(Hotelling)) install.packages("Hotelling")
library("Hotelling")

library("edgeR")
library("limma")
library("biomaRt")
library("GenomicFeatures")

### Perform T-Test analysis
countsforanalysis_df = tmmcounts_hgnc_pcexmtgenes_normedwithall_df
expr_min = 0.50

countsforanalysis_df[,"trimmed_mean"] = apply(countsforanalysis_df, 1, function(x) mean(x, trim = 0.1))
countsforanalysis_df = countsforanalysis_df[rownames(head(countsforanalysis_df[order(countsforanalysis_df["trimmed_mean"], decreasing = T),], expr_min*nrow(countsforanalysis_df))),]
countsforanalysis_df = subset(countsforanalysis_df, select=-c(trimmed_mean))
colnames(countsforanalysis_df)

combinedsample_foldchange_df = as.data.frame(countsforanalysis_df)
combinedsample_foldchange_df$RA_Intact = rowMedians(cbind(countsforanalysis_df[,"RA3408_Intact"], countsforanalysis_df[,"RA3514_Intact"], countsforanalysis_df[,"RA3520_Intact"]))
combinedsample_foldchange_df$RA_Lib = rowMedians(cbind(countsforanalysis_df[,"RA3408_Lib"], countsforanalysis_df[,"RA3514_Lib"], countsforanalysis_df[,"RA3520_Lib"]))
combinedsample_foldchange_df$RA_LibFlavo = rowMedians(cbind(countsforanalysis_df[,"OA3491_LibFlavo"], countsforanalysis_df[,"RA3514_LibFlavo"], countsforanalysis_df[,"RA3520_LibFlavo"]))
combinedsample_foldchange_df$RA_SubA = rowMedians(cbind(countsforanalysis_df[,"RA3408_SubA"], countsforanalysis_df[,"RA3514_SubA"], countsforanalysis_df[,"RA3520_SubA"]))
combinedsample_foldchange_df$OA_Intact = rowMedians(cbind(countsforanalysis_df[,"OA3491_Intact"], countsforanalysis_df[,"OA3494_Intact"], countsforanalysis_df[,"OA3496_Intact"]))
combinedsample_foldchange_df$OA_Lib = rowMedians(cbind(countsforanalysis_df[,"OA3491_Lib"], countsforanalysis_df[,"OA3494_Lib"], countsforanalysis_df[,"OA3496_Lib"]))
combinedsample_foldchange_df$OA_LibFlavo = rowMedians(cbind(countsforanalysis_df[,"OA3491_LibFlavo"], countsforanalysis_df[,"OA3494_LibFlavo"], countsforanalysis_df[,"OA3496_LibFlavo"]))
combinedsample_foldchange_df$OA_SubA = rowMedians(cbind(countsforanalysis_df[,"OA3491_SubA"], countsforanalysis_df[,"OA3494_SubA"], countsforanalysis_df[,"OA3496_SubA"]))

RA_Intact_sample_df = countsforanalysis_df[,Reduce(intersect, list(grep("RA", colnames(countsforanalysis_df), value = TRUE), grep("Intact", colnames(countsforanalysis_df), value = TRUE)))]
OA_Intact_sample_df = countsforanalysis_df[,Reduce(intersect, list(grep("OA", colnames(countsforanalysis_df), value = TRUE), grep("Intact", colnames(countsforanalysis_df), value = TRUE)))]
RA_Lib_sample_df = countsforanalysis_df[,Reduce(intersect, list(grep("RA", colnames(countsforanalysis_df), value = TRUE), grep("Lib$", colnames(countsforanalysis_df), value = TRUE)))]
OA_Lib_sample_df = countsforanalysis_df[,Reduce(intersect, list(grep("OA", colnames(countsforanalysis_df), value = TRUE), grep("Lib$", colnames(countsforanalysis_df), value = TRUE)))]
RA_LibFlavo_sample_df = countsforanalysis_df[,Reduce(intersect, list(grep("RA", colnames(countsforanalysis_df), value = TRUE), grep("LibFlavo", colnames(countsforanalysis_df), value = TRUE)))]
OA_LibFlavo_sample_df = countsforanalysis_df[,Reduce(intersect, list(grep("OA", colnames(countsforanalysis_df), value = TRUE), grep("LibFlavo", colnames(countsforanalysis_df), value = TRUE)))]
RA_SubA_sample_df = countsforanalysis_df[,Reduce(intersect, list(grep("RA", colnames(countsforanalysis_df), value = TRUE), grep("SubA", colnames(countsforanalysis_df), value = TRUE)))]
OA_SubA_sample_df = countsforanalysis_df[,Reduce(intersect, list(grep("OA", colnames(countsforanalysis_df), value = TRUE), grep("SubA", colnames(countsforanalysis_df), value = TRUE)))]

RAOA_Intact_hotelling = hotelling.test(RA_Intact_sample_df, OA_Intact_sample_df, shrinkage = FALSE)
RAOA_Intact_hotelling$pval
RAOA_Lib_hotelling = hotelling.test(RA_Lib_sample_df, OA_Lib_sample_df, shrinkage = FALSE)
RAOA_Lib_hotelling$pval
RAOA_LibFlavo_hotelling = hotelling.test(RA_LibFlavo_sample_df, OA_LibFlavo_sample_df, shrinkage = FALSE)
RAOA_LibFlavo_hotelling$pval
RAOA_SubA_hotelling = hotelling.test(RA_SubA_sample_df, OA_SubA_sample_df, shrinkage = FALSE)
RAOA_SubA_hotelling$pval

RA_IntactLib_hotelling = hotelling.test(RA_Intact_sample_df, RA_Lib_sample_df, shrinkage = FALSE)
RA_IntactLib_hotelling$pval
RA_IntactLibFlavo_hotelling = hotelling.test(RA_Intact_sample_df, RA_LibFlavo_sample_df, shrinkage = FALSE)
RA_IntactLibFlavo_hotelling$pval
RA_IntactSubA_hotelling = hotelling.test(RA_Intact_sample_df, RA_SubA_sample_df, shrinkage = FALSE)
RA_IntactSubA_hotelling$pval

RA_LibLibFlavo_hotelling = hotelling.test(RA_Lib_sample_df, RA_LibFlavo_sample_df, shrinkage = FALSE)
RA_LibLibFlavo_hotelling$pval
RA_LibFlavoSubA_hotelling = hotelling.test(RA_LibFlavo_sample_df, RA_SubA_sample_df, shrinkage = FALSE)
RA_LibFlavoSubA_hotelling$pval

OA_IntactLib_hotelling = hotelling.test(OA_Intact_sample_df, OA_Lib_sample_df, shrinkage = FALSE)
OA_IntactLib_hotelling$pval
OA_IntactLibFlavo_hotelling = hotelling.test(OA_Intact_sample_df, OA_LibFlavo_sample_df, shrinkage = FALSE)
OA_IntactLibFlavo_hotelling$pval
OA_IntactSubA_hotelling = hotelling.test(OA_Intact_sample_df, OA_SubA_sample_df, shrinkage = FALSE)
OA_IntactSubA_hotelling$pval

OA_LibLibFlavo_hotelling = hotelling.test(OA_Lib_sample_df, OA_LibFlavo_sample_df, shrinkage = FALSE)
OA_LibLibFlavo_hotelling$pval
OA_LibFlavoSubA_hotelling = hotelling.test(OA_LibFlavo_sample_df, OA_SubA_sample_df, shrinkage = FALSE)
OA_LibFlavoSubA_hotelling$pval

