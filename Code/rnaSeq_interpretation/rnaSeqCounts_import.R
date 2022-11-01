###############################################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ddply")

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
if(!require(qdapRegex)) install.packages("qdapRegex")
library(qdapRegex)
if(!require(stringr)) install.packages("stringr")
library(stringr)

library("edgeR")
library("limma")
library("biomaRt")
library("GenomicFeatures")

### Import RNAseq Files

### Create gene x sample matrix
combinedSample_list = list()
sample_dir = "/Volumes/Portable Drive/Documents/Graduate/Data/Projects/2021.Disagg_RA_Analysis/"
sample_files = list.files(path = sample_dir, pattern = ".counts$", recursive = TRUE)

for (fp in sample_files){
  sampleTable_fp = paste0(sample_dir, fp)
  sampleTable_df = as.data.frame(read.delim(sampleTable_fp, row.names = 1, header = FALSE))
  condition_name = substring(basename(fp), 1, 2)
  sample_number = substring(basename(fp), 3, 6)
  treatment_name = substring(basename(fp), 8, 12)
  treatment_name = str_replace(treatment_name, "_F$", ".F")
  treatment_name = str_replace(treatment_name, "Lib_S$", "Liber")
  treatment_name = str_replace(treatment_name, "-S$", "St")
  treatment_name = str_replace(treatment_name, "_A$", ".A")
  experiment_name = substring(fp, 8, 8)
  combined_name = paste0(condition_name, "_", sample_number, "_", treatment_name, "_", experiment_name)
  print(combined_name)
  colnames(sampleTable_df) = c(combined_name)
  combinedSample_list[[combined_name]] = sampleTable_df
}

combinedSample_df = Reduce(function(x,y) transform(merge(x, y, by="row.names", all = TRUE), row.names=Row.names, Row.names = NULL), combinedSample_list)
samplecombined_matrix = as.matrix(combinedSample_df)

## Make Combined amd QC Matrices
QCcombined_df = combinedSample_df

# remove samples with failed QC
failed_samples = c("")
samplecombined_matrix = combinedSample_df[6:(nrow(samplecombined_matrix)), !colnames(combinedSample_df) %in% failed_samples]
samplecombined_df = as.data.frame(samplecombined_matrix)

###  Convert ENSEMBL ID to HGNC Symbol for counts

# Create HGNC gene ID symbol matrix
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes = rownames(samplecombined_df)
genes = gsub("\\..*", "",genes)
idname_df = getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart=mart)

# Match HGNC symbol to ensembl ID
samplecombined_hgnc_df = samplecombined_df
samplecombined_hgnc_df$ensembl_gene_id = gsub("\\..*", "", rownames(samplecombined_df))
samplecombined_hgnc_df = merge(samplecombined_hgnc_df, idname_df, by="ensembl_gene_id")

# Remove non-mapped and duplicate HGNC IDs (take max count from duplicated HGNC IDs)
id_null = which(nchar(samplecombined_hgnc_df$hgnc_symbol)==0)
samplecombined_hgnc_df = samplecombined_hgnc_df[-id_null,]

samplecombined_hgnc_df = ddply(samplecombined_hgnc_df,"hgnc_symbol",numcolwise(max))
rownames(samplecombined_hgnc_df) = samplecombined_hgnc_df$hgnc_symbol
samplecombined_hgnc_df = subset(samplecombined_hgnc_df, select=-c(hgnc_symbol))

# Filter for protein-coding only
all_genes = rownames(samplecombined_hgnc_df)
genefunction_df = getBM(attributes= c("hgnc_symbol", "transcript_biotype"), filters= "hgnc_symbol", values=all_genes, mart=mart)
pc_genes = genefunction_df$hgnc_symbol[genefunction_df$transcript_biotype == "protein_coding"]

samplecombined_hgnc_pcgenes_df = samplecombined_hgnc_df[pc_genes,]

# Check percent of MT encoded genes and remove MT genes
mtgenes_index = Reduce(union, list(grep("^MT-", row.names(samplecombined_hgnc_pcgenes_df)),
                                   grep("^MTATP", row.names(samplecombined_hgnc_pcgenes_df)),
                                   grep("^MTND", row.names(samplecombined_hgnc_pcgenes_df)),
                                   grep("^MTCO", row.names(samplecombined_hgnc_pcgenes_df)),
                                   grep("^MTCY", row.names(samplecombined_hgnc_pcgenes_df))))


mtgenes_df = samplecombined_hgnc_pcgenes_df[mtgenes_index,]
mtgene_percent_df = as.data.frame(colSums(mtgenes_df)/colSums(samplecombined_hgnc_pcgenes_df))

samplecombined_hgnc_exmtgenes_df = samplecombined_hgnc_pcgenes_df[-mtgenes_index,]

### Normalize counts matrix to geTMM

# create gtf data; # make TxDb from GTF file and get gene information
gtf_fp <- "/Volumes/Portable Drive/Documents/Graduate/Data/Annotation/GTF/gencode.v37.annotation.gtf"
txdb <- makeTxDbFromGFF(gtf_fp)
all.genes <- genes(txdb)
all.genes.lengths_df <- as.data.frame(all.genes@ranges@width)
rownames(all.genes.lengths_df) <- all.genes@ranges@NAMES

##  FOR GENE LENGTH FILE
# Match HGNC symbol to ensembl ID - FOR GENE LENGTH FILE
all.genes.lengths_hgnc_df = all.genes.lengths_df
all.genes.lengths_hgnc_df$ensembl_gene_id = gsub("\\..*", "", rownames(all.genes.lengths_df))
all.genes.lengths_hgnc_df = merge(all.genes.lengths_hgnc_df, idname_df, by="ensembl_gene_id")

# Remove non-mapped and duplicate HGNC IDs (take max count from duplicated HGNC IDs)  - FOR GENE LENGTH FILE
id_null = which(nchar(all.genes.lengths_hgnc_df$hgnc_symbol)==0)
all.genes.lengths_hgnc_df = all.genes.lengths_hgnc_df[-id_null,]

all.genes.lengths_hgnc_df = ddply(all.genes.lengths_hgnc_df,"hgnc_symbol",numcolwise(max))
rownames(all.genes.lengths_hgnc_df) = all.genes.lengths_hgnc_df$hgnc_symbol
all.genes.lengths_hgnc_df = subset(all.genes.lengths_hgnc_df, select=-c(hgnc_symbol))

## FOR BOTH ALL GENES DF and EX-MT GENES DFs
## import list of gene names
my.genes = c(rownames(samplecombined_hgnc_df))
my.notmtgenes = c(rownames(samplecombined_hgnc_exmtgenes_df))
my.pcgenes = c(rownames(samplecombined_hgnc_pcgenes_df))

## get the length of each of those genes
my.genes.lengths_df = subset(all.genes.lengths_hgnc_df, rownames(all.genes.lengths_hgnc_df) %in% my.genes)
my.notmtgenes.lengths_df = subset(all.genes.lengths_hgnc_df, rownames(all.genes.lengths_hgnc_df) %in% my.notmtgenes)
my.pcgenes.lengths_df = subset(all.genes.lengths_hgnc_df, rownames(all.genes.lengths_hgnc_df) %in% my.pcgenes)

## RPK normalize
rpkcounts_hgnc_df = samplecombined_hgnc_df
genelength_kb=c(my.genes.lengths_df[,1])
genelength_kb=genelength_kb/1000

for(y in colnames(samplecombined_hgnc_df)){
  rawsample=c(samplecombined_hgnc_df[,y])
  rpksample=rawsample/genelength_kb
  rpkcounts_hgnc_df[,y]=rpksample
}

rpkcounts_hgnc_exmtgenes_df = samplecombined_hgnc_exmtgenes_df
genelength_kb=c(my.notmtgenes.lengths_df[,1])
genelength_kb=genelength_kb/1000

for(y in colnames(samplecombined_hgnc_exmtgenes_df)){
  rawsample=c(samplecombined_hgnc_exmtgenes_df[,y])
  rpksample=rawsample/genelength_kb
  rpkcounts_hgnc_exmtgenes_df[,y]=rpksample
}

rpkcounts_hgnc_pcgenes_df = samplecombined_hgnc_pcgenes_df
genelength_kb=c(my.pcgenes.lengths_df[,1])
genelength_kb=genelength_kb/1000

for(y in colnames(samplecombined_hgnc_pcgenes_df)){
  rawsample=c(samplecombined_hgnc_pcgenes_df[,y])
  rpksample=rawsample/genelength_kb
  rpkcounts_hgnc_pcgenes_df[,y]=rpksample
}

# TMM normalize
rpkcounts_hgnc_dge = DGEList(counts = rpkcounts_hgnc_df)
calcNormFactors(rpkcounts_hgnc_dge, method = "TMM")
tmmcounts_hgnc_matrix = cpm(rpkcounts_hgnc_dge)
tmmcounts_hgnc_df = as.data.frame(tmmcounts_hgnc_matrix)

rpkcounts_hgnc_exmtgenes_dge = DGEList(counts = rpkcounts_hgnc_exmtgenes_df)
calcNormFactors(rpkcounts_hgnc_exmtgenes_dge, method = "TMM")
tmmcounts_hgnc_exmtgenes_matrix = cpm(rpkcounts_hgnc_exmtgenes_dge)
tmmcounts_hgnc_exmtgenes_df = as.data.frame(tmmcounts_hgnc_exmtgenes_matrix)

rpkcounts_hgnc_pcgenes_dge = DGEList(counts = rpkcounts_hgnc_pcgenes_df)
calcNormFactors(rpkcounts_hgnc_pcgenes_dge, method = "TMM")
tmmcounts_hgnc_pcgenes_matrix = cpm(rpkcounts_hgnc_pcgenes_dge)
tmmcounts_hgnc_pcgenes_df = as.data.frame(tmmcounts_hgnc_pcgenes_matrix)

# Create protein-coding and ex-MT gene gene dataframe with normalization done on entire set

tmmcounts_hgnc_pcexmtgenes_normedwithall_df = tmmcounts_hgnc_df[pc_genes,]

mtgenes_index = Reduce(union, list(grep("^MT-", row.names(tmmcounts_hgnc_pcexmtgenes_normedwithall_df)),
                                   grep("^MTATP", row.names(tmmcounts_hgnc_pcexmtgenes_normedwithall_df)),
                                   grep("^MTND", row.names(tmmcounts_hgnc_pcexmtgenes_normedwithall_df)),
                                   grep("^MTCO", row.names(tmmcounts_hgnc_pcexmtgenes_normedwithall_df)),
                                   grep("^MTCY", row.names(tmmcounts_hgnc_pcexmtgenes_normedwithall_df))))

tmmcounts_hgnc_pcexmtgenes_normedwithall_df = tmmcounts_hgnc_pcexmtgenes_normedwithall_df[-mtgenes_index,]

### Export workbook for external analysis
setwd("/Users/bartonprideaux/OneDrive/Documents/Education/Graduate_School/UCSD/Research/Projects/RA-Disagg/Analysis/Round3/xls_output")
RAcounts_pcexmtgenes_wb=createWorkbook("RAOA_Disagg_counts_pcexmtgenes.xlsx")
addWorksheet(RAcounts_pcexmtgenes_wb, "TMMcounts")
writeData(RAcounts_pcexmtgenes_wb, "TMMcounts", data.frame(tmmcounts_hgnc_pcexmtgenes_normedwithall_df), rowNames = TRUE)
saveWorkbook(RAcounts_pcexmtgenes_wb, file = "RAOA_Disagg_counts_pcexmtgenes.xlsx", overwrite = TRUE)