
library(Seurat)
library(ggplot2)
library(RColorBrewer)


## Load data
seuset.rna <- readRDS('/home/ajl1213/Projects/PD/SciAdv_Data/snRNAseq/02_SeuratProcess/SeuratObjects/PD.SN.snRNA.postAlign.Anno.rds')
dir.create('CountByCellType')

#==================================================================================================================================================
## Extract raw count
df <- NULL
for (cellType in levels(seuset.rna@meta.data$ident)){
    print(cellType)
    idx_list <- rownames(seuset.rna@meta.data)[which(seuset.rna@meta.data$ident == cellType)]
    cellTypeData <- as.matrix(seuset.rna@assays$RNA@counts)[, as.character(idx_list)]
    sumVec <- rowSums(cellTypeData)
    df <- cbind(df, sumVec)
}

colnames(df) <- levels(seuset.rna@meta.data$ident)
df.final <- data.frame(geneID=rownames(df), df)
write.table(df.final, file='CountByCellType/RawCount.byCellType.txt', row.names=F, col.names=T, sep='\t', quote=F)

##==================================================================================================================================================
## quantile normalize
library(preprocessCore)

values <- as.matrix(df)
qq_values <- normalize.quantiles.robust(values)
colnames(qq_values) <- colnames(df)
rownames(qq_values) <- rownames(df)

df.final <- data.frame(geneID=row.names(qq_values), qq_values)
write.table(df.final, file='CountByCellType/RawCount.byCellType.qq.txt', row.name=F, col.names=T, sep="\t", quote=F)

##==================================================================================================================================================
## TPM
valMat <- read.table('CountByCellType/RawCount.byCellType.qq.txt', header=T, row.names=1)
totalCounts <- colSums(valMat)
TPM_mat <- valMat / totalCount * 1000000

df.final <- data.frame(geneID=row.names(TPM_mat), TPM_mat)
write.table(df.final, file='CountByCellType/snRNA.TPM.byCellType.txt', row.name=F, col.names=T, sep='\t', quote=F)

##==================================================================================================================================================
## log2 transform
log2_valmat <- log2(qq_values+1)
df.final <- data.frame(geneID=row.names(log2_valmat), log2_valmat)
write.table(df.final, file='CountByCellType/RawCount.byCellType.qq.log2.txt', row.name=F, col.names=T, sep="\t", quote=F)

#==================================================================================================================================================
## z-transform
valMat <- log2_valmat
get_zscore <- function(x){(x - mean(x)) /sd(x)}
zmat <- apply(valMat, 1, get_zscore)
zmat <- t(zmat)
zmat[is.na(zmat)] <- 0
df.final <- data.frame(geneID=row.names(zmat), zmat)
write.table(df.final, file='CountByCellType/RawCount.byCellType.qq.log2.zscore.txt', row.name=F, col.names=T, sep="\t", quote=F)


##==================================================================================================================================================
## expression ratio
valMat <- qq_values

get_ratio <- function(x){x/sum(x)*100}
ratio_mat <- apply(valMat, 1, get_ratio) 
ratio_mat <- t(ratio_mat)
ratio_mat[is.na(ratio_mat)] <- 0

df.final <- data.frame(geneID=row.names(ratio_mat), ratio_mat)
write.table(df.final, file='CountByCellType/ExpRatio.byCellType.txt', row.name=F, col.names=T, sep="\t", quote=F)


