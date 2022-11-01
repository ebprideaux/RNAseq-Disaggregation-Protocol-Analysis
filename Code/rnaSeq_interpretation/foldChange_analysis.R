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
library("reshape2")

### Perform Fold Change Analysis
countsforanalysis_df = tmmcounts_hgnc_pcexmtgenes_normedwithall_df
expr_min = 0.50

countsforanalysis_df[,"trimmed_mean"] = apply(countsforanalysis_df, 1, function(x) mean(x, trim = 0.1))
countsforanalysis_df = countsforanalysis_df[rownames(head(countsforanalysis_df[order(countsforanalysis_df["trimmed_mean"], decreasing = T),], expr_min*nrow(countsforanalysis_df))),]
countsforanalysis_df = subset(countsforanalysis_df, select=-c(trimmed_mean))
colnames(countsforanalysis_df)

combinedsample_foldchange_df = data.frame(matrix(ncol = 0, nrow = nrow(countsforanalysis_df)))
rownames(combinedsample_foldchange_df) = row.names(countsforanalysis_df)

# select needed vectors
RA_intact_vector = Reduce(intersect, list(grep("RA", colnames(countsforanalysis_df), value = TRUE), grep("Intac", colnames(countsforanalysis_df),  value = TRUE)))
OA_intact_vector = Reduce(intersect, list(grep("OA", colnames(countsforanalysis_df), value = TRUE), grep("Intac", colnames(countsforanalysis_df),  value = TRUE)))
RA_libera_vector = Reduce(intersect, list(grep("RA", colnames(countsforanalysis_df), value = TRUE), grep("Liber", colnames(countsforanalysis_df),  value = TRUE)))
OA_libera_vector = Reduce(intersect, list(grep("OA", colnames(countsforanalysis_df), value = TRUE), grep("Liber", colnames(countsforanalysis_df),  value = TRUE)))
RA_libfla_vector = Reduce(intersect, list(grep("RA", colnames(countsforanalysis_df), value = TRUE), grep("Lib.F", colnames(countsforanalysis_df),  value = TRUE)))
OA_libfla_vector = Reduce(intersect, list(grep("OA", colnames(countsforanalysis_df), value = TRUE), grep("Lib.F", colnames(countsforanalysis_df),  value = TRUE)))
RA_suba_vector = Reduce(intersect, list(grep("RA", colnames(countsforanalysis_df), value = TRUE), grep("Sub.A", colnames(countsforanalysis_df),  value = TRUE)))
OA_suba_vector = Reduce(intersect, list(grep("OA", colnames(countsforanalysis_df), value = TRUE), grep("Sub.A", colnames(countsforanalysis_df),  value = TRUE)))

# RA vs OA, for each of intact vs intact, lib vs lib, libflavo vs libflavo, suba vs suba
combinedsample_foldchange_df$RAIntactVsOAIntact = log2(rowMedians(as.matrix(countsforanalysis_df[,RA_intact_vector]))/(rowMedians(as.matrix(countsforanalysis_df[,OA_intact_vector]))))
combinedsample_foldchange_df$RALibVsOALib = log2(rowMedians(as.matrix(countsforanalysis_df[,RA_libera_vector]))/(rowMedians(as.matrix(countsforanalysis_df[,OA_libera_vector]))))
combinedsample_foldchange_df$RALibFlavoVsOALibFlavo = log2(rowMedians(as.matrix(countsforanalysis_df[,RA_libfla_vector]))/(rowMedians(as.matrix(countsforanalysis_df[,OA_libfla_vector]))))
combinedsample_foldchange_df$RASubAVsOASubA = log2(rowMedians(as.matrix(countsforanalysis_df[,RA_suba_vector]))/(rowMedians(as.matrix(countsforanalysis_df[,OA_suba_vector]))))

# OA treatment vs intact
combinedsample_foldchange_df$OALibVsOAIntact = log2(rowMedians(as.matrix(countsforanalysis_df[,OA_libera_vector]))/(rowMedians(as.matrix(countsforanalysis_df[,OA_intact_vector]))))
combinedsample_foldchange_df$OALibFlavoVsOAIntact = log2(rowMedians(as.matrix(countsforanalysis_df[,OA_libfla_vector]))/(rowMedians(as.matrix(countsforanalysis_df[,OA_intact_vector]))))
combinedsample_foldchange_df$OASubAVsOAIntact = log2(rowMedians(as.matrix(countsforanalysis_df[,OA_suba_vector]))/(rowMedians(as.matrix(countsforanalysis_df[,OA_intact_vector]))))

# RA treatment vs intact
combinedsample_foldchange_df$RALibVsRAIntact = log2(rowMedians(as.matrix(countsforanalysis_df[,RA_libera_vector]))/(rowMedians(as.matrix(countsforanalysis_df[,RA_intact_vector]))))
combinedsample_foldchange_df$RALibFlavoVsRAIntact = log2(rowMedians(as.matrix(countsforanalysis_df[,RA_libfla_vector]))/(rowMedians(as.matrix(countsforanalysis_df[,RA_intact_vector]))))
combinedsample_foldchange_df$RASubAVsRAIntact = log2(rowMedians(as.matrix(countsforanalysis_df[,RA_suba_vector]))/(rowMedians(as.matrix(countsforanalysis_df[,RA_intact_vector]))))

rownames(combinedsample_foldchange_df[order(combinedsample_foldchange_df$RALibVsRAIntact, decreasing = T),])[1:10]
rownames(combinedsample_foldchange_df[order(combinedsample_foldchange_df$RALibFlavoVsRAIntact, decreasing = T),])[1:10]
rownames(combinedsample_foldchange_df[order(combinedsample_foldchange_df$RASubAVsRAIntact, decreasing = T),])[1:10]


# Plot charts
setwd("/Users/bartonprideaux/OneDrive/Documents/Education/Graduate_School/UCSD/Research/Projects/RA-Disagg/Analysis/Round3/Figures/FoldChange")

png("Log2(FC)_OALib_OALibFlavo_OAIntact.png", width = 385, height = 261, units = "px")
plot(combinedsample_foldchange_df$OALibVsOAIntact, combinedsample_foldchange_df$OALibFlavoVsOAIntact, 
     col = "grey", main = "OA", xlab = "Log2(Liberase / Intact)", ylab = "Log2(LibFlavo / Intact)", cex.lab = 1.2, cex = 1.0,
     xlim=c(-5,5), ylim=c(-5,5), yaxs = "i", xaxs = "i")
dev.off()

png("Log2(FC)_OALib_OASubA_OAIntact.png", width = 385, height = 261, units = "px")
plot(combinedsample_foldchange_df$OALibVsOAIntact, combinedsample_foldchange_df$OASubAVsOAIntact, 
     col = "grey", main = "OA", xlab = "Log2(Liberase / Intact)", ylab = "Log2(SubA / Intact)", cex.lab = 1.2, cex = 1.0,
     xlim=c(-5,5), ylim=c(-5,5), yaxs = "i", xaxs = "i")
dev.off()

png("Log2(FC)_OALibFlavo_OASubA_OAIntact.png", width = 385, height = 261, units = "px")
plot(combinedsample_foldchange_df$OALibFlavoVsOAIntact, combinedsample_foldchange_df$OASubAVsOAIntact, 
     col = "grey", main = "OA", xlab = "Log2(LibFlavo / Intact)", ylab = "Log2(SubA / Intact)", cex.lab = 1.2, cex = 1.0,
     xlim=c(-5,5), ylim=c(-5,5), yaxs = "i", xaxs = "i")
dev.off()

png("Log2(FC)_RALib_RALibFlavo_RAIntact.png", width = 385, height = 261, units = "px")
plot(combinedsample_foldchange_df$RALibVsRAIntact, combinedsample_foldchange_df$RALibFlavoVsRAIntact, 
     col = "grey", main = "RA", xlab = "Log2(Liberase / Intact)", ylab = "Log2(LibFlavo / Intact)", cex.lab = 1.2, cex = 1.0,
     xlim=c(-5,5), ylim=c(-5,5), yaxs = "i", xaxs = "i")
dev.off()

png("Log2(FC)_RALib_RASubA_RAIntact.png", width = 385, height = 261, units = "px")
plot(combinedsample_foldchange_df$RALibVsRAIntact, combinedsample_foldchange_df$RASubAVsRAIntact, 
     col = "grey", main = "RA", xlab = "Log2(Liberase / Intact)", ylab = "Log2(SubA / Intact)", cex.lab = 1.2, cex = 1.0,
     xlim=c(-5,5), ylim=c(-5,5), yaxs = "i", xaxs = "i")
dev.off()

png("Log2(FC)_RALibFlavo_RASubA_RAIntact.png", width = 385, height = 261, units = "px")
plot(combinedsample_foldchange_df$RALibFlavoVsRAIntact, combinedsample_foldchange_df$RASubAVsRAIntact, 
     col = "grey", main = "RA", xlab = "Log2(LibFlavo / Intact)", ylab = "Log2(SubA / Intact)", cex.lab = 1.2, cex = 1.0,
     xlim=c(-5,5), ylim=c(-5,5), yaxs = "i", xaxs = "i")
dev.off()

#######

png("Log2(FC)_RA-LibIntactVsOA-LibIntact.png", width = 900, height = 450, units = "px")
plot(combinedsample_foldchange_df$RALibVsRAIntact, combinedsample_foldchange_df$OALibVsOAIntact, 
     col = "grey", xlab = "RA-Lib/Intact", ylab = "OA-Lib/Intact", cex = 1.0,
     xlim=c(-5,5), ylim=c(-5,5), yaxs = "i", xaxs = "i")
dev.off()

png("Log2(FC)_RA-LibFlavoIntactVsOA-LibFlavoIntact.png", width = 900, height = 450, units = "px")
plot(combinedsample_foldchange_df$RALibFlavoVsRAIntact, combinedsample_foldchange_df$OALibFlavoVsOAIntact, 
     col = "grey", xlab = "RA-LibFlavo/Intact", ylab = "OA-LibFlavo/Intact", cex = 1.0,
     xlim=c(-5,5), ylim=c(-5,5), yaxs = "i", xaxs = "i")
dev.off()

png("Log2(FC)_RA-SubAIntactVsOA-SubAIntact.png", width = 900, height = 450, units = "px")
plot(combinedsample_foldchange_df$RASubAVsRAIntact, combinedsample_foldchange_df$OASubAVsOAIntact, 
     col = "grey", xlab = "RA-SubA/Intact", ylab = "OA-SubA/Intact", cex = 1.0,
     xlim=c(-5,5), ylim=c(-5,5), yaxs = "i", xaxs = "i")
dev.off()


# histograms
setwd("/Users/bartonprideaux/OneDrive/Documents/Education/Graduate_School/UCSD/Research/Projects/RA-Disagg/Analysis/Round3/Figures/Histograms")

RAvsOAcomparison_df = combinedsample_foldchange_df[,c("RAIntactVsOAIntact", "RALibVsOALib", "RALibFlavoVsOALibFlavo", "RASubAVsOASubA")]
gg = melt(RAvsOAcomparison_df)
png("RAvsOA_Log2(FC)_hist.png", width = 875, height = 550, units = "px")
ggplot(gg, aes(x=value, fill=variable)) + geom_histogram(binwidth=0.02) + facet_grid(variable~.) + 
  scale_x_continuous(breaks = seq(-4, 4, 1), lim = c(-4, 4)) + scale_y_continuous(breaks = seq(0, 400, 100), lim = c(0, 400))
dev.off()

OAprotcomparison_df = combinedsample_foldchange_df[,c("OALibVsOAIntact", "OALibFlavoVsOAIntact", "OASubAVsOAIntact")]
colnames(OAprotcomparison_df) = c("Liberase", "LibFlavo", "SubA")
gg = melt(OAprotcomparison_df)
png("OAprot_Log2(FC)_hist.png", width = 525, height = 195, units = "px")
ggplot(gg, aes(x=value, fill=variable)) + 
  geom_histogram(binwidth=0.02) + 
  facet_grid(variable~.) + 
  scale_x_continuous(breaks = seq(-4, 4, 1), lim = c(-4, 4)) + 
  scale_y_continuous(breaks = seq(0, 300, 100), lim = c(0, 300)) + 
  ggtitle("OA - Log2(Treatment / Intact)")
dev.off()

RAprotcomparison_df = combinedsample_foldchange_df[,c("RALibVsRAIntact", "RALibFlavoVsRAIntact", "RASubAVsRAIntact")]
colnames(RAprotcomparison_df) = c("Liberase", "LibFlavo", "SubA")
gg = melt(RAprotcomparison_df)
png("RAprot_Log2(FC)_hist.png", width = 525, height = 195, units = "px")
ggplot(gg, aes(x=value, fill=variable)) + 
  geom_histogram(binwidth=0.02) + 
  facet_grid(variable~.) + 
  scale_x_continuous(breaks = seq(-4, 4, 1), lim = c(-4, 4)) + 
  scale_y_continuous(breaks = seq(0, 300, 100), lim = c(0, 300)) + 
  ggtitle("RA - Log2(Treatment / Intact)")
dev.off()

# calculate statistics for each dist
OA_LibvIntact_var = var(combinedsample_foldchange_df$OALibVsOAIntact,na.rm=TRUE)
OA_LibFlavovIntact_var = var(combinedsample_foldchange_df$OALibFlavoVsOAIntact,na.rm=TRUE)
OA_SubavIntact_var = var(combinedsample_foldchange_df$OASubAVsOAIntact,na.rm=TRUE)

RA_LibvIntact_var = var(combinedsample_foldchange_df$RALibVsRAIntact[abs(combinedsample_foldchange_df$RALibVsRAIntact) < 9999999],na.rm=TRUE)
RA_LibFlavovIntact_var = var(combinedsample_foldchange_df$RALibFlavoVsRAIntact[abs(combinedsample_foldchange_df$RALibFlavoVsRAIntact) < 9999999],na.rm=TRUE)
RA_SubavIntact_var = var(combinedsample_foldchange_df$RASubAVsRAIntact[abs(combinedsample_foldchange_df$RASubAVsRAIntact) < 9999999],na.rm=TRUE)


# sig test for histogram distribution (Wilcox Test)
OALib_OALibFlavoVsOAIntact_FC_MWresult = wilcox.test(combinedsample_foldchange_df$OALibVsOAIntact, combinedsample_foldchange_df$OALibFlavoVsOAIntact, alternative = "greater")
print(OALib_OALibFlavoVsOAIntact_FC_MWresult$p.value)
OALib_OASubAVsOAIntact_FC_MWresult = wilcox.test(combinedsample_foldchange_df$OALibVsOAIntact, combinedsample_foldchange_df$OASubAVsOAIntact, alternative = "greater")
print(OALib_OASubAVsOAIntact_FC_MWresult$p.value)
OALibFlavo_OASubAVsOAIntact_FC_MWresult = wilcox.test(combinedsample_foldchange_df$OALibFlavoVsOAIntact, combinedsample_foldchange_df$OASubAVsOAIntact, alternative = "greater")
print(OALibFlavo_OASubAVsOAIntact_FC_MWresult$p.value)

# get gene list for FC change greater than X for each
fc_min = log2(1.5)

OALib_highFCgenes_vector = as.vector(rownames(combinedsample_foldchange_df[combinedsample_foldchange_df$OALibVsOAIntact >= fc_min, , drop = FALSE]))
length(OALib_highFCgenes_vector)
clipr::write_clip(OALib_highFCgenes_vector)

OALibFlavo_highFCgenes_vector = as.vector(rownames(combinedsample_foldchange_df[combinedsample_foldchange_df$OALibFlavoVsOAIntact >= fc_min, , drop = FALSE]))
length(OALibFlavo_highFCgenes_vector)
clipr::write_clip(OALibFlavo_highFCgenes_vector)

OASubA_highFCgenes_vector = as.vector(rownames(combinedsample_foldchange_df[combinedsample_foldchange_df$OASubAVsOAIntact >= fc_min, , drop = FALSE]))
length(OASubA_highFCgenes_vector)
clipr::write_clip(OASubA_highFCgenes_vector)

RALib_highFCgenes_vector = as.vector(rownames(combinedsample_foldchange_df[combinedsample_foldchange_df$RALibVsRAIntact >= fc_min, , drop = FALSE]))
length(RALib_highFCgenes_vector)
clipr::write_clip(RALib_highFCgenes_vector)

RALibFlavo_highFCgenes_vector = as.vector(rownames(combinedsample_foldchange_df[combinedsample_foldchange_df$RALibFlavoVsRAIntact >= fc_min, , drop = FALSE]))
length(RALibFlavo_highFCgenes_vector)
clipr::write_clip(RALibFlavo_highFCgenes_vector)

RASubA_highFCgenes_vector = as.vector(rownames(combinedsample_foldchange_df[combinedsample_foldchange_df$RASubAVsRAIntact >= fc_min, , drop = FALSE]))
length(RASubA_highFCgenes_vector)
clipr::write_clip(RASubA_highFCgenes_vector)

### RA vs OA
# Intact
RAIntactVsOAIntact_DEG_vector = as.vector(rownames(combinedsample_foldchange_df[combinedsample_foldchange_df$RAIntactVsOAIntact >= fc_min, , drop = FALSE]))
length(RAIntactVsOAIntact_DEG_vector)
clipr::write_clip(sort(RAIntactVsOAIntact_DEG_vector))
# Lib
RALibVsOALib_DEG_vector = as.vector(rownames(combinedsample_foldchange_df[combinedsample_foldchange_df$RALibVsOALib >= fc_min, , drop = FALSE]))
length(RALibVsOALib_DEG_vector)
clipr::write_clip(sort(RALibVsOALib_DEG_vector))
# Lib Flavo
RALibFlavoVsOALibFlavo_DEG_vector = as.vector(rownames(combinedsample_foldchange_df[combinedsample_foldchange_df$RALibFlavoVsOALibFlavo >= fc_min, , drop = FALSE]))
length(RALibFlavoVsOALibFlavo_DEG_vector)
clipr::write_clip(sort(RALibFlavoVsOALibFlavo_DEG_vector))
# Sub A
RASubAVsOASubA_DEG_vector = as.vector(rownames(combinedsample_foldchange_df[combinedsample_foldchange_df$RASubAVsOASubA >= fc_min, , drop = FALSE]))
length(RASubAVsOASubA_DEG_vector)
clipr::write_clip(sort(RASubAVsOASubA_DEG_vector))


### Identify top FC genes
numbergenes=1000

# Setup dataframes
topdiffgenes_df = combinedsample_foldchange_df
topdiffgenes_df = combinedsample_foldchange_df[combinedsample_foldchange_df$RAIntactVsOAIntact != "Inf", , drop = FALSE]
RAIntactVsOAIntact_FC_topgenes_vector=as.vector(rownames(topdiffgenes_df[order(topdiffgenes_df$RAIntactVsOAIntact, decreasing = T),])[1:numbergenes])
RAIntactVsOAIntact_FC_topgenes_df = countsforanalysis_df[RAIntactVsOAIntact_FC_topgenes_vector, grep("Intact", colnames(countsforanalysis_df))]
clipr::write_clip(RAIntactVsOAIntact_FC_topgenes_vector)

topdiffgenes_df = combinedsample_foldchange_df
topdiffgenes_df = combinedsample_foldchange_df[combinedsample_foldchange_df$RALibVsOALib != "Inf", , drop = FALSE]
RALibVsOALib_FC_topgenes_vector=as.vector(rownames(topdiffgenes_df[order(topdiffgenes_df$RALibVsOALib, decreasing = T),])[1:numbergenes])
RALibVsOALib_FC_topgenes_df = countsforanalysis_df[RALibVsOALib_FC_topgenes_vector, grep("Lib$", colnames(countsforanalysis_df))]
clipr::write_clip(RALibVsOALib_FC_topgenes_vector)

topdiffgenes_df = combinedsample_foldchange_df
topdiffgenes_df = combinedsample_foldchange_df[combinedsample_foldchange_df$RALibFlavoVsOALibFlavo != "Inf", , drop = FALSE]
RALibFlavoVsOALibFlavo_FC_topgenes_vector=as.vector(rownames(topdiffgenes_df[order(topdiffgenes_df$RALibFlavoVsOALibFlavo, decreasing = T),])[1:numbergenes])
RALibFlavoVsOALibFlavo_FC_topgenes_df = countsforanalysis_df[RALibFlavoVsOALibFlavo_FC_topgenes_vector, grep("LibFlavo", colnames(countsforanalysis_df))]
clipr::write_clip(RALibFlavoVsOALibFlavo_FC_topgenes_vector)

topdiffgenes_df = combinedsample_foldchange_df
topdiffgenes_df = combinedsample_foldchange_df[combinedsample_foldchange_df$RASubAVsOASubA != "Inf", , drop = FALSE]
RASubAVsOASubA_FC_topgenes_vector=as.vector(rownames(topdiffgenes_df[order(topdiffgenes_df$RASubAVsOASubA, decreasing = T),])[1:numbergenes])
RASubAVsOASubA_FC_topgenes_df = countsforanalysis_df[RASubAVsOASubA_FC_topgenes_vector, grep("SubA", colnames(countsforanalysis_df))]
clipr::write_clip(RASubAVsOASubA_FC_topgenes_vector)