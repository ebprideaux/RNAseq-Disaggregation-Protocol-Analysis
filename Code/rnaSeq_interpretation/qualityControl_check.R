### Quality Control ###
QCcombined_df = as.data.frame(QCcombined_df)
head(QCcombined_df)

###  Convert ENSEMBL ID to HGNC Symbol for counts

# Create HGNC gene ID symbol matrix
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes = rownames(QCcombined_df)
genes = gsub("\\..*", "",genes)
idname_df = getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart=mart)

# Match HGNC symbol to ensembl ID
samplecombined_hgnc_df = QCcombined_df
samplecombined_hgnc_df$ensembl_gene_id = gsub("\\..*", "", rownames(QCcombined_df))
samplecombined_hgnc_df = merge(samplecombined_hgnc_df, idname_df, by="ensembl_gene_id")

# Remove non-mapped and duplicate HGNC IDs (take max count from duplicated HGNC IDs)
id_null = which(nchar(samplecombined_hgnc_df$hgnc_symbol)==0)
samplecombined_hgnc_df = samplecombined_hgnc_df[-id_null,]

samplecombined_hgnc_df = ddply(samplecombined_hgnc_df,"hgnc_symbol",numcolwise(max))
rownames(samplecombined_hgnc_df) = samplecombined_hgnc_df$hgnc_symbol
samplecombined_hgnc_df = subset(samplecombined_hgnc_df, select=-c(hgnc_symbol))

# MT-encoded genes
mtgenes_index = Reduce(union, list(grep("^MT-", row.names(samplecombined_hgnc_df)),
                                   grep("^MTATP", row.names(samplecombined_hgnc_df)),
                                   grep("^MTND", row.names(samplecombined_hgnc_df)),
                                   grep("^MTCO", row.names(samplecombined_hgnc_df)),
                                   grep("^MTCY", row.names(samplecombined_hgnc_df))))

mtgenes_df = samplecombined_hgnc_df[mtgenes_index,]
mtgene_percent_df = as.data.frame(colSums(mtgenes_df)/colSums(samplecombined_hgnc_df))

RAintact_MTavg = mean(mtgene_percent_df[Reduce(intersect, list(grep("RA", rownames(mtgene_percent_df), value = TRUE),grep("Intac", rownames(mtgene_percent_df), value = TRUE))),])
OAintact_MTavg = mean(mtgene_percent_df[Reduce(intersect, list(grep("OA", rownames(mtgene_percent_df), value = TRUE),grep("Intac", rownames(mtgene_percent_df), value = TRUE))),])
RALib_MTavg = mean(mtgene_percent_df[Reduce(intersect, list(grep("RA", rownames(mtgene_percent_df), value = TRUE),grep("Liber", rownames(mtgene_percent_df), value = TRUE))),])
OALib_MTavg = mean(mtgene_percent_df[Reduce(intersect, list(grep("OA", rownames(mtgene_percent_df), value = TRUE),grep("Liber", rownames(mtgene_percent_df), value = TRUE))),])
RALibFlavo_MTavg = mean(mtgene_percent_df[Reduce(intersect, list(grep("RA", rownames(mtgene_percent_df), value = TRUE),grep("Lib.F", rownames(mtgene_percent_df), value = TRUE))),])
OALibFlavo_MTavg = mean(mtgene_percent_df[Reduce(intersect, list(grep("OA", rownames(mtgene_percent_df), value = TRUE),grep("Lib.F", rownames(mtgene_percent_df), value = TRUE))),])
RASubA_MTavg = mean(mtgene_percent_df[Reduce(intersect, list(grep("RA", rownames(mtgene_percent_df), value = TRUE),grep("Sub.A", rownames(mtgene_percent_df), value = TRUE))),])
OASubA_MTavg = mean(mtgene_percent_df[Reduce(intersect, list(grep("OA", rownames(mtgene_percent_df), value = TRUE),grep("Sub.A", rownames(mtgene_percent_df), value = TRUE))),])


# rRNA
rRNAgenes_index = Reduce(union, list(grep("^RNA5S", row.names(samplecombined_hgnc_df)), 
                                     grep("^MT-RNR", row.names(samplecombined_hgnc_df))
))

rRNAgenes_df = samplecombined_hgnc_df[rRNAgenes_index,]
rRNA_percent_df = as.data.frame(colSums(rRNAgenes_df)/colSums(samplecombined_hgnc_df))

MTrRNA_QC_df = cbind.data.frame(mtgene_percent_df, rRNA_percent_df)
colnames(MTrRNA_QC_df) = c("%MT-Encoded", "rRNA")

# Standard QC
standardQC_df = QCcombined_df[nrow(QCcombined_df):1,]
standardQC_df["TotalReads",] = colSums(standardQC_df)
standardQC_df["#Genes(Counts>=1)",] = apply(standardQC_df, 2, function(x) sum(x >= 1))
standardQC_df["#Genes(Counts>=4)",] = apply(standardQC_df, 2, function(x) sum(x >= 4))
standardQC_df["#Genes(Counts>=10)",] = apply(standardQC_df, 2, function(x) sum(x >= 10))
standardQC_df["#Genes(Counts>=100)",] = apply(standardQC_df, 2, function(x) sum(x >= 100))
standardQC_df = tail(standardQC_df, n=10)
standardQC_df = as.data.frame(t(standardQC_df))

QC_df = cbind.data.frame(standardQC_df, MTrRNA_QC_df)

summary(QC_df)

clusters = kmeans(QC_df[,c("%MT-Encoded")], 3)
qc_cluster_vector = as.vector(clusters$cluster)
qc_cluster_df = as.data.frame(qc_cluster_vector)
qc_cluster_df$samples = rownames(QC_df)

mtqc_vector = as.vector(mtgene_percent_df[,1])

## Export QC workbook
setwd("/Users/bartonprideaux/OneDrive/Documents/Education/Graduate_School/UCSD/Research/Projects/RA-Disagg/Data/QualityControl/Round3")

QC_wb=createWorkbook("qc_metrics.xlsx")
addWorksheet(QC_wb, "qc_metrics")
writeData(QC_wb, "qc_metrics", data.frame(QC_df), rowNames = TRUE)
saveWorkbook(QC_wb, file = "qc_metrics.xlsx", overwrite = TRUE)