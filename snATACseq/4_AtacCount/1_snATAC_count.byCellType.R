
library(Seurat)
library(Signac)
library(ggplot2)
library(RColorBrewer)


## Load data
seuset.atac <- readRDS('/home/ajl1213/Projects/PD/SciAdv_Data/snATACseq/2_SignacProcess/SeuratObjects/PD.SN.snATAC.postAlign.geneExp.label.anno.rds')
dir.create('CountByCellType')

#==================================================================================================================================================
## Extract raw count
df <- NULL
for (cellType in levels(seuset.atac@meta.data$ident)){
    print(cellType)
    idx_list <- rownames(seuset.atac@meta.data)[which(seuset.atac@meta.data$ident == cellType)]
    cellTypeData <- as.matrix(seuset.atac@assays$RNA@counts)[, as.character(idx_list)]
    sumVec <- rowSums(cellTypeData)
    df <- cbind(df, sumVec)
}

colnames(df) <- levels(seuset.atac@meta.data$ident)
df.final <- data.frame(geneID=rownames(df), df)
write.table(df.final, file='CountByCellType/RawCount.snATAC.byCellType.txt', row.names=F, col.names=T, sep='\t', quote=F)

##==================================================================================================================================================
## quantile normalize
library(preprocessCore)

values <- as.matrix(df)
qq_values <- normalize.quantiles.robust(values)
colnames(qq_values) <- colnames(df)
rownames(qq_values) <- rownames(df)

df.final <- data.frame(geneID=row.names(qq_values), qq_values)
write.table(df.final, file='CountByCellType/RawCount.snATAC.byCellType.qq.txt', row.name=F, col.names=T, sep="\t", quote=F)

##==================================================================================================================================================
## TPM
valMat <- read.table('CountByCellType/RawCount.snATAC.byCellType.qq.txt', header=T, row.names=1)
totalCounts <- colSums(valMat)
TPM_mat <- valMat / totalCounts * 1000000

df.final <- data.frame(geneID=row.names(TPM_mat), TPM_mat)
write.table(df.final, file='CountByCellType/snATAC.TPM.byCellType.txt', row.name=F, col.names=T, sep='\t', quote=F)

##==================================================================================================================================================
## expression ratio
valMat <- qq_values

get_ratio <- function(x){x/sum(x)*100}
ratio_mat <- apply(valMat, 1, get_ratio) 
ratio_mat <- t(ratio_mat)
ratio_mat[is.na(ratio_mat)] <- 0

df.final <- data.frame(geneID=row.names(ratio_mat), ratio_mat)
write.table(df.final, file='CountByCellType/ExpRatio.snATAC.byCellType.txt', row.name=F, col.names=T, sep="\t", quote=F)


