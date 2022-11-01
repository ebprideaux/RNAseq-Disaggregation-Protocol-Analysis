# Load Key Libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ReactomePA")
BiocManager::install("clusterProfiler")
if(!require("venn")) install.packages("venn")
if(!require("ggVennDiagram")) install.packages("ggVennDiagram")
if(!require(ggplot2)) install.packages("ggplot2")
BiocManager::install("org.Hs.eg.db")
if(!require(org.Hs.eg.db)) install.packages("org.Hs.eg.db")


library("ggplot2")
library("ggVennDiagram")
library("ReactomePA")
library("clusterProfiler")
library("openxlsx")
library("plyr")

######################################################################## Venn Diagram ########################################################################
setwd("/Users/bartonprideaux/OneDrive/Documents/Education/Graduate_School/UCSD/Research/Projects/RA-Disagg/Analysis/Round3/Figures/VennDiagram")

### RA vs OA under each treatment combined ###
RAvsOA_allTreatments_DEGs_list = list(Intact=RAIntactVsOAIntact_DEG_vector, Liberase=RALibVsOALib_DEG_vector, 
                                      LibFlavo=RALibFlavoVsOALibFlavo_DEG_vector, SubA=RASubAVsOASubA_DEG_vector)

RAvsOA_allTreatments_DEGs_chart = ggVennDiagram(RAvsOA_allTreatments_DEGs_list, color = "black", lwd = 0.8, lty = 1) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

png("RAvsOA_allTreatments_DEGs_VennDiagram.png", width = 575, height = 280, units = "px")
plot(RAvsOA_allTreatments_DEGs_chart)
dev.off()

RAvsOA_allTreatments_DEGs_tibl = process_region_data(Venn(RAvsOA_allTreatments_DEGs_list))

a_vector = sort(as.vector(RAvsOA_allTreatments_DEGs_tibl$item[[6]]))
length(a_vector)
clipr::write_clip(sort(a_vector))

b_vector = sort(as.vector(RAvsOA_allTreatments_DEGs_tibl$item[[13]]))
length(b_vector)
clipr::write_clip(sort(b_vector))

c_vector = sort(as.vector(RAvsOA_allTreatments_DEGs_tibl$item[[7]]))
length(c_vector)
clipr::write_clip(sort(c_vector))

### Create dataframe of DEGs
RAvsOA_allTreatments_DEGs_df = t(plyr::ldply(RAvsOA_allTreatments_DEGs_list, rbind))
colnames(RAvsOA_allTreatments_DEGs_df) = RAvsOA_allTreatments_DEGs_df[".id",]
RAvsOA_allTreatments_DEGs_df = RAvsOA_allTreatments_DEGs_df[-c(1),]

setwd("/Users/bartonprideaux/Library/CloudStorage/OneDrive-Personal/Documents/Education/Graduate_School/UCSD/Research/Projects/RA-Disagg/Analysis/Round3/xls_output")

DEG_wb=createWorkbook("DEGs.xlsx")
addWorksheet(DEG_wb, "DEG_list")
writeData(DEG_wb, "DEG_list", RAvsOA_allTreatments_DEGs_df, rowNames = FALSE, colNames = TRUE)
saveWorkbook(DEG_wb, file = "DEGs.xlsx", overwrite = TRUE)

### RA vs OA under each treatment separated ###
setwd("/Users/bartonprideaux/Library/CloudStorage/OneDrive-Personal/Documents/Education/Graduate_School/UCSD/Research/Projects/RA-Disagg/Analysis/Round3/Figures/VennDiagram")

## Intact vs Liberase ##
# Create chart
RAvsOA_intactLib_DEGs_list = list(Intact=RAIntactVsOAIntact_DEG_vector, Liberase=RALibVsOALib_DEG_vector)
RAvsOA_intactLib_DEGs_chart = ggVennDiagram(RAvsOA_intactLib_DEGs_list, color = "black", lwd = 0.8, lty = 1) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

png("RAvsOA_intactLib_DEGs_VennDiagram.png", width = 430, height = 210, units = "px")
plot(RAvsOA_intactLib_DEGs_chart)
dev.off()

# Define missed genes
RAvsOA_intactLib_DEGs_tibl = process_region_data(Venn(RAvsOA_intactLib_DEGs_list))
RAvsOA_intactLib_missedDEGs_vector = sort(as.vector(RAvsOA_intactLib_DEGs_tibl$item[[1]]))
length(RAvsOA_intactLib_missedDEGs_vector)
clipr::write_clip(sort(RAvsOA_intactLib_missedDEGs_vector))


## Intact vs LibFlavo ##
# Create chart
RAvsOA_intactLibFlavo_DEGs_list = list(Intact=RAIntactVsOAIntact_DEG_vector, LibFlavo=RALibFlavoVsOALibFlavo_DEG_vector)
RAvsOA_intactLibFlavo_DEGs_chart = ggVennDiagram(RAvsOA_intactLibFlavo_DEGs_list, color = "black", lwd = 0.8, lty = 1) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

png("RAvsOA_intactLibFlavo_DEGs_VennDiagram.png", width = 430, height = 210, units = "px")
plot(RAvsOA_intactLibFlavo_DEGs_chart)
dev.off()

# Define missed genes
RAvsOA_intactLibFlavo_DEGs_tibl = process_region_data(Venn(RAvsOA_intactLibFlavo_DEGs_list))
RAvsOA_intactLibFlavo_missedDEGs_vector = sort(as.vector(RAvsOA_intactLibFlavo_DEGs_tibl$item[[1]]))
length(RAvsOA_intactLibFlavo_missedDEGs_vector)
clipr::write_clip(sort(RAvsOA_intactLibFlavo_missedDEGs_vector))


## Intact vs Sub A ##
# Create chart
RAvsOA_intactSubA_DEGs_list = list(Intact=RAIntactVsOAIntact_DEG_vector, SubA=RASubAVsOASubA_DEG_vector)
RAvsOA_intactSubA_DEGs_chart = ggVennDiagram(RAvsOA_intactSubA_DEGs_list, color = "black", lwd = 0.8, lty = 1) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

png("RAvsOA_intactSubA_DEGs_VennDiagram.png", width = 430, height = 210, units = "px")
plot(RAvsOA_intactSubA_DEGs_chart)
dev.off()

# Define missed genes
RAvsOA_intactSubA_DEGs_tibl = process_region_data(Venn(RAvsOA_intactSubA_DEGs_list))
RAvsOA_intactSubA_missedDEGs_vector = sort(as.vector(RAvsOA_intactSubA_DEGs_tibl$item[[1]]))
length(RAvsOA_intactSubA_missedDEGs_vector)
clipr::write_clip(sort(RAvsOA_intactSubA_missedDEGs_vector))


### Create missed pathways list for RA vs OA under each treatment relative to Intact ###
adjPval_max = 0.01
## Create vectors of pathways for each treatment
# Intact
RAIntactVsOAIntact_DEGid_vector = bitr(RAIntactVsOAIntact_DEG_vector, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)
RAvsOA_intact_pathways_obj = enrichPathway(gene=RAIntactVsOAIntact_DEGid_vector$ENTREZID, pAdjustMethod = "BH")

RAvsOA_intact_pathways_description_vector = as.vector(RAvsOA_intact_pathways_obj@result[RAvsOA_intact_pathways_obj@result$p.adjust < adjPval_max,]$Description)
RAvsOA_intact_pathways_geneRatio_vector = as.vector(RAvsOA_intact_pathways_obj@result[RAvsOA_intact_pathways_obj@result$p.adjust < adjPval_max,]$GeneRatio)
RAvsOA_intact_pathways_pvalue_vector = as.vector(RAvsOA_intact_pathways_obj@result[RAvsOA_intact_pathways_obj@result$p.adjust < adjPval_max,]$pvalue)
RAvsOA_intact_pathways_adjPvalue_vector = as.vector(RAvsOA_intact_pathways_obj@result[RAvsOA_intact_pathways_obj@result$p.adjust < adjPval_max,]$p.adjust)
RAvsOA_intact_pathways_df = data.frame("Description" = RAvsOA_intact_pathways_description_vector, "Gene Ratio" = RAvsOA_intact_pathways_geneRatio_vector, "P-Value" = RAvsOA_intact_pathways_pvalue_vector, "Adj. P-Value" = RAvsOA_intact_pathways_adjPvalue_vector)

# Liberase
RALibVsOALib_DEGid_vector = bitr(RALibVsOALib_DEG_vector, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)
RAvsOA_lib_pathways_obj = enrichPathway(gene=RALibVsOALib_DEGid_vector$ENTREZID, pAdjustMethod = "BH")

RAvsOA_lib_pathways_description_vector = as.vector(RAvsOA_lib_pathways_obj@result[RAvsOA_lib_pathways_obj@result$p.adjust < adjPval_max,]$Description)
RAvsOA_lib_pathways_geneRatio_vector = as.vector(RAvsOA_lib_pathways_obj@result[RAvsOA_lib_pathways_obj@result$p.adjust < adjPval_max,]$GeneRatio)
RAvsOA_lib_pathways_pvalue_vector = as.vector(RAvsOA_lib_pathways_obj@result[RAvsOA_lib_pathways_obj@result$p.adjust < adjPval_max,]$pvalue)
RAvsOA_lib_pathways_adjPvalue_vector = as.vector(RAvsOA_lib_pathways_obj@result[RAvsOA_lib_pathways_obj@result$p.adjust < adjPval_max,]$p.adjust)
RAvsOA_lib_pathways_df = data.frame("Description" = RAvsOA_lib_pathways_description_vector, "Gene Ratio" = RAvsOA_lib_pathways_geneRatio_vector, "P-Value" = RAvsOA_lib_pathways_pvalue_vector, "Adj. P-Value" = RAvsOA_lib_pathways_adjPvalue_vector)

# LibFlavo
RALibFlavoVsOALibFlavo_DEGid_vector = bitr(RALibFlavoVsOALibFlavo_DEG_vector, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)
RAvsOA_libFlavo_pathways_obj = enrichPathway(gene=RALibFlavoVsOALibFlavo_DEGid_vector$ENTREZID, pAdjustMethod = "BH")

RAvsOA_libFlavo_pathways_description_vector = as.vector(RAvsOA_libFlavo_pathways_obj@result[RAvsOA_libFlavo_pathways_obj@result$p.adjust < adjPval_max,]$Description)
RAvsOA_libFlavo_pathways_geneRatio_vector = as.vector(RAvsOA_libFlavo_pathways_obj@result[RAvsOA_libFlavo_pathways_obj@result$p.adjust < adjPval_max,]$GeneRatio)
RAvsOA_libFlavo_pathways_pvalue_vector = as.vector(RAvsOA_libFlavo_pathways_obj@result[RAvsOA_libFlavo_pathways_obj@result$p.adjust < adjPval_max,]$pvalue)
RAvsOA_libFlavo_pathways_adjPvalue_vector = as.vector(RAvsOA_libFlavo_pathways_obj@result[RAvsOA_libFlavo_pathways_obj@result$p.adjust < adjPval_max,]$p.adjust)
RAvsOA_libFlavo_pathways_df = data.frame("Description" = RAvsOA_libFlavo_pathways_description_vector, "Gene Ratio" = RAvsOA_libFlavo_pathways_geneRatio_vector, "P-Value" = RAvsOA_libFlavo_pathways_pvalue_vector, "Adj. P-Value" = RAvsOA_libFlavo_pathways_adjPvalue_vector)

# SubA
RASubAVsOASubA_DEGid_vector = bitr(RASubAVsOASubA_DEG_vector, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)
RAvsOA_subA_pathways_obj = enrichPathway(gene=RASubAVsOASubA_DEGid_vector$ENTREZID, pAdjustMethod = "BH")

RAvsOA_subA_pathways_description_vector = as.vector(RAvsOA_subA_pathways_obj@result[RAvsOA_subA_pathways_obj@result$p.adjust < adjPval_max,]$Description)
RAvsOA_subA_pathways_geneRatio_vector = as.vector(RAvsOA_subA_pathways_obj@result[RAvsOA_subA_pathways_obj@result$p.adjust < adjPval_max,]$GeneRatio)
RAvsOA_subA_pathways_pvalue_vector = as.vector(RAvsOA_subA_pathways_obj@result[RAvsOA_subA_pathways_obj@result$p.adjust < adjPval_max,]$pvalue)
RAvsOA_subA_pathways_adjPvalue_vector = as.vector(RAvsOA_subA_pathways_obj@result[RAvsOA_subA_pathways_obj@result$p.adjust < adjPval_max,]$p.adjust)
RAvsOA_subA_pathways_df = data.frame("Description" = RAvsOA_subA_pathways_description_vector, "Gene Ratio" = RAvsOA_subA_pathways_geneRatio_vector, "P-Value" = RAvsOA_subA_pathways_pvalue_vector, "Adj. P-Value" = RAvsOA_subA_pathways_adjPvalue_vector)

## Create vectors of missed pathways
liberase_missedPaths_description_vector = setdiff(RAvsOA_intact_pathways_description_vector, RAvsOA_lib_pathways_description_vector)
liberase_extraPaths_description_vector = setdiff(RAvsOA_lib_pathways_description_vector, RAvsOA_intact_pathways_description_vector)
liberase_differentPaths_maxLength = max(length(liberase_missedPaths_description_vector), length(liberase_extraPaths_description_vector))
liberase_missedPaths_description_vector = c(liberase_missedPaths_description_vector, rep("", liberase_differentPaths_maxLength - length(liberase_missedPaths_description_vector)))
liberase_extraPaths_description_vector = c(liberase_extraPaths_description_vector, rep("", liberase_differentPaths_maxLength - length(liberase_extraPaths_description_vector)))
liberase_differentPaths_df = data.frame("Missed Paths" = liberase_missedPaths_description_vector, 
                                        "Extra Paths Found" = liberase_extraPaths_description_vector)
clipr::write_clip(liberase_missedPaths_description_vector)
clipr::write_clip(liberase_extraPaths_description_vector)

libFlavo_missedPaths_description_vector = setdiff(RAvsOA_intact_pathways_description_vector, RAvsOA_libFlavo_pathways_description_vector)
libFlavo_extraPaths_description_vector = setdiff(RAvsOA_libFlavo_pathways_description_vector, RAvsOA_intact_pathways_description_vector)
libFlavo_differentPaths_maxLength = max(length(libFlavo_missedPaths_description_vector), length(libFlavo_extraPaths_description_vector))
libFlavo_missedPaths_description_vector = c(libFlavo_missedPaths_description_vector, rep("", libFlavo_differentPaths_maxLength - length(libFlavo_missedPaths_description_vector)))
libFlavo_extraPaths_description_vector = c(libFlavo_extraPaths_description_vector, rep("", libFlavo_differentPaths_maxLength - length(libFlavo_extraPaths_description_vector)))
libFlavo_differentPaths_df = data.frame("Missed Paths" = libFlavo_missedPaths_description_vector, 
                                        "Extra Paths Found" = libFlavo_extraPaths_description_vector)
clipr::write_clip(libFlavo_missedPaths_description_vector)
clipr::write_clip(libFlavo_extraPaths_description_vector)

subA_missedPaths_description_vector = setdiff(RAvsOA_intact_pathways_description_vector, RAvsOA_subA_pathways_description_vector)
subA_extraPaths_description_vector = setdiff(RAvsOA_subA_pathways_description_vector, RAvsOA_intact_pathways_description_vector)
subA_differentPaths_maxLength = max(length(subA_missedPaths_description_vector), length(subA_extraPaths_description_vector))
subA_missedPaths_description_vector = c(subA_missedPaths_description_vector, rep("", subA_differentPaths_maxLength - length(subA_missedPaths_description_vector)))
subA_extraPaths_description_vector = c(subA_extraPaths_description_vector, rep("", subA_differentPaths_maxLength - length(subA_extraPaths_description_vector)))
subA_differentPaths_df = data.frame("Missed Paths" = subA_missedPaths_description_vector, 
                                        "Extra Paths Found" = subA_extraPaths_description_vector)
clipr::write_clip(subA_missedPaths_description_vector)
clipr::write_clip(subA_extraPaths_description_vector)

### Export pathways for external analysis
setwd("/Users/bartonprideaux/Library/CloudStorage/OneDrive-Personal/Documents/Education/Graduate_School/UCSD/Research/Projects/RA-Disagg/Analysis/Round3/Figures/ReactomePathways")

pathwayEnrichment_wb=createWorkbook("enrichedPathways.xlsx")
addWorksheet(pathwayEnrichment_wb, "RAOA_intact")
writeData(pathwayEnrichment_wb, "RAOA_intact", RAvsOA_intact_pathways_df, rowNames = TRUE)
addWorksheet(pathwayEnrichment_wb, "RAOA_liberase")
writeData(pathwayEnrichment_wb, "RAOA_liberase", RAvsOA_lib_pathways_df, rowNames = TRUE)
addWorksheet(pathwayEnrichment_wb, "RAOA_libflavo")
writeData(pathwayEnrichment_wb, "RAOA_libflavo", RAvsOA_libFlavo_pathways_df, rowNames = TRUE)
addWorksheet(pathwayEnrichment_wb, "RAOA_subA")
writeData(pathwayEnrichment_wb, "RAOA_subA", RAvsOA_subA_pathways_df, rowNames = TRUE)
addWorksheet(pathwayEnrichment_wb, "liberase_diffPathways")
writeData(pathwayEnrichment_wb, "liberase_diffPathways", liberase_differentPaths_df, rowNames = TRUE)
addWorksheet(pathwayEnrichment_wb, "libFlavo_diffPathways")
writeData(pathwayEnrichment_wb, "libFlavo_diffPathways", libFlavo_differentPaths_df, rowNames = TRUE)
addWorksheet(pathwayEnrichment_wb, "subA_diffPathways")
writeData(pathwayEnrichment_wb, "subA_diffPathways", subA_differentPaths_df, rowNames = TRUE)
saveWorkbook(pathwayEnrichment_wb, file = "enrichedPathways.xlsx", overwrite = TRUE)


### Missed pathways Venn Diagram ###
missedPathways_list = list(Liberase=liberase_missedPaths_description_vector, 
                           LibFlavo=libFlavo_missedPaths_description_vector, 
                           SubA=subA_missedPaths_description_vector)

missedPathways_list_chart = ggVennDiagram(missedPathways_list, color = "black", lwd = 0.8, lty = 1) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

missedPathways_tibl = process_region_data(Venn(missedPathways_list))

specificMissedPaths_lib_vector = sort(as.vector(missedPathways_tibl$item[[1]]))
clipr::write_clip(specificMissedPaths_lib_vector)