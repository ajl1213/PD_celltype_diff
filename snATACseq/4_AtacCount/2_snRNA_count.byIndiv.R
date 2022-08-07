
library(Seurat)
library(Signac)
library(ggplot2)
library(RColorBrewer)


## Load data
seuset.atac <- readRDS('/home/ajl1213/Projects/PD/SciAdv_Data/snATACseq/2_SignacProcess/SeuratObjects/PD.SN.snATAC.postAlign.geneExp.label.anno.rds')
dir.create('CountByIndiv')

#==================================================================================================================================================
## Extract raw count
for (cellType in levels(seuset.atac@meta.data$ident)){
    print(cellType)
    df <- NULL
    for (sampleID in unique(seuset.atac@meta.data$orig.ident)){
	print(paste0(cellType,',',sampleID))
	idx_list <- rownames(seuset.atac@meta.data)[which(seuset.atac@meta.data$ident == cellType & seuset.atac@meta.data$orig.ident==sampleID)]

	if (length(idx_list) > 1){
	cellTypeData <- as.matrix(seuset.atac@assays$RNA@counts)[, as.character(idx_list)]
	sumVec <- rowSums(cellTypeData)
	}
	if (length(idx_list) == 1){
	cellTypeData <- as.matrix(seuset.atac@assays$RNA@counts)[, as.character(idx_list)]
	sumVec <- cellTypeData 
	}
	if (length(idx_list) == 0){
	sumVec <- rep(0, dim(seuset.atac@assays$RNA@counts)[1])
	names(sumVec) <- rownames(seuset.atac@assays$RNA@counts)
	}
	df <- cbind(df, sumVec)
    }

    colnames(df) <- unique(seuset.atac@meta.data$orig.ident) 

    print(colSums(df))
    df.final <- data.frame(geneID=rownames(df), df) 
    write.table(df.final, file=paste0('CountByIndiv/', cellType,'.snATAC.RawCount.txt'), row.names=F, col.names=T, sep='\t', quote=F)
}

##==================================================================================================================================================
## quantile normalize
library(preprocessCore)

for (cellType in levels(seuset.atac@meta.data$ident)){
    print(cellType)

    df <- read.table(paste0('CountByIndiv/', cellType,'.snATAC.RawCount.txt'), header=T, row.names=1)

    values <- matrix(as.numeric(unlist(df)), nrow=nrow(df))
    qq_values <- normalize.quantiles.robust(values)
    colnames(qq_values) <- colnames(df)
    rownames(qq_values) <- rownames(df)

    df.final <- data.frame(geneID=row.names(qq_values), qq_values)
    write.table(df.final, file=paste0('CountByIndiv/', cellType, '.snATAC.RawCount.qq.txt'), row.name=F, col.names=T, sep="\t", quote=F)

}


##==================================================================================================================================================
## TPM
for (cellType in levels(seuset.atac@meta.data$ident)){
    print(cellType)

    df <- read.table(paste0('CountByIndiv/', cellType,'.snATAC.RawCount.qq.txt'), header=T, row.names=1)
    countSum <- colSums(df)
    TPM_DF <- df / countSum * 1000000

    df.final <- data.frame(geneID=rownames(TPM_DF), TPM_DF)
    write.table(df.final, file=paste0('CountByIndiv/', cellType, '.snATAC.TPM.txt'), row.name=F, col.names=T, sep="\t", quote=F)
}


##==================================================================================================================================================
## z-transform

for (cellType in levels(seuset.atac@meta.data$ident)){
    print(cellType)

    df <- read.table(paste0('CountByIndiv/', cellType, '.snATAC.TPM.txt'), header=T, row.names=1)

    get_zscore <- function(x){(x - mean(x)) /sd(x)}
    zmat <- apply(df, 1, get_zscore)
    zmat <- t(zmat)
    zmat[is.na(zmat)] <- 0
    df.final <- data.frame(geneID=row.names(zmat), zmat)
    write.table(df.final, file=paste0('CountByIndiv/', cellType, '.snATAC.TPM.zscore.txt'), row.name=F, col.names=T, sep="\t", quote=F)

}





