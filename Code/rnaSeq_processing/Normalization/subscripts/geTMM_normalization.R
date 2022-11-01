#load libraries; install if necessary
library(limma)
library(edgeR)
library(biomaRt)
library(tidyverse)
library(plyr)
library(GenomicFeatures)

# import passed variables
args = commandArgs(trailingOnly=TRUE)
countfilelist_fp = args[1]
gtf_fp = args[2]
rpk_dir = args[3]
getmm_dir = args[4]

# create gtf data; # make TxDb from GTF file and get gene information
txdb = makeTxDbFromGFF(gtf_fp)
all.genes = genes(txdb)

all.genes.lengths_df = as.data.frame(all.genes@ranges@width)
ensembl_ids = all.genes@ranges@NAMES
all.genes.lengths_df$ensembl_gene_index_id = ensembl_ids
trimmedEnsembl_ids = sub('\\..*', '', ensembl_ids)
all.genes.lengths_df$ensembl_gene_id = trimmedEnsembl_ids

# Create gene ID symbol matrix
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes = trimmedEnsembl_ids
idname_df = getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart=mart)

all.genes.lengths_df = merge(all.genes.lengths_df, idname_df, by="ensembl_gene_id")
all.genes.lengths_df = all.genes.lengths_df[!duplicated(all.genes.lengths_df$hgnc_symbol),]

######### First, take files and convert to HGNC symbols #########

# Import raw counts
filelist = read.table(countfilelist_fp, sep = "\n", col.names = "files")
rawcounts_dge = readDGE(filelist, header=F)
normcounts_df = rawcounts_dge$counts
normcounts_df = head(normcounts_df, -5)
print("DGE import complete")

# Format columns
samples = c(rawcounts_dge$samples[,1])
cleanSamples = sub(".*/", "", samples)
colnames(normcounts_df) = cleanSamples

### Change gene id to gene symbol
# Import and format norm counts matrix
normcounts_df = as.data.frame(normcounts_df)
normcounts_df$ensembl_gene_index_id = rownames(normcounts_df)
#normcounts_df$ensembl_gene_id = sub('\\..*', '', normcounts_df$ensembl_gene_index_id)
rownames(normcounts_df) = NULL

# Match gene symbols to IDs
normcounts_df$ensembl_gene_index_id
all.genes.lengths_df$ensembl_gene_index_id
combined_df = merge(normcounts_df, all.genes.lengths_df, by="ensembl_gene_index_id")
#combined_df = merge(normcounts_df, all.genes.lengths_df, by="ensembl_gene_id")
id_null = which(nchar(combined_df$hgnc_symbol)==0)
combined_df = combined_df[-id_null,]

# Remove gene ID column and change gene symbol column to rowname
#clean_df = subset(combined_df, select=-c(ensembl_gene_id))

# Sum and remove duplicates where HGNC Gene Symbol in multiple locations
nondup_df = ddply(combined_df,"hgnc_symbol",numcolwise(sum))
rownames(nondup_df) = nondup_df[,"hgnc_symbol"]
nondup_df = subset(nondup_df, select = c(cleanSamples))

######### Second, convert raw counts to RPK counts #########
# import counts files to be transformed
rawcounts_df = as.data.frame(nondup_df)

## import list of gene names
my.genes = c(rownames(rawcounts_df))
my.genes = my.genes[1:(length(my.genes))]

## get the length of each of those genes
my.genes.lengths_df = all.genes.lengths_df[all.genes.lengths_df$hgnc_symbol %in% my.genes,]

## perform RPK transformation
rpkcounts_df = rawcounts_df
genelength_kb=c(my.genes.lengths_df$`all.genes@ranges@width`)
genelength_kb=genelength_kb/1000

for(y in colnames(rawcounts_df)){
  rawsample=c(rawcounts_df[,y])
  rpksample=rawsample/genelength_kb
  rpkcounts_df[,y]=rpksample
}

## separate into individual files
for(y in colnames(rpkcounts_df)){
  old_fn=y
  new_fn=paste0(old_fn,".RPK")
  new_fp=paste0(rpk_dir,"/",new_fn)
  names_df=rownames(rpkcounts_df)
  values_df=rpkcounts_df[,y]
  output_df=cbind.data.frame(names_df,values_df,stringsAsFactors=FALSE)
  write.table(output_df, file=new_fp, sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE)
}

######### Third, convert RPK counts to geTMM counts #########

# Import rpk counts
file_list = paste0(rpk_dir,"/",list.files(path = rpk_dir, pattern = ".RPK$", recursive = TRUE))
rawcounts_dge = readDGE(file_list, header=F)
rawcounts_dge = calcNormFactors(rawcounts_dge, method = "TMM")
normcounts_df = cpm(rawcounts_dge)

# Format columns
samples <- c(rawcounts_dge$samples[,1])
colnames(normcounts_df) <- samples

# Format columns
samples = c(rawcounts_dge$samples[,1])
cleanSamples = sub(".*/", "", samples)
colnames(normcounts_df) = cleanSamples

y=colnames(normcounts_df)[1]
## separate into individual files
for(y in colnames(normcounts_df)){
  old_fn=y
  new_fn=paste0(old_fn,".geTMM")
  new_fp=paste0(getmm_dir,"/",new_fn)
  names_df=rownames(normcounts_df)
  values_df=normcounts_df[,y]
  output_df=cbind.data.frame(names_df,values_df,stringsAsFactors=FALSE)
  write.table(output_df, file=new_fp, sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE)
}
