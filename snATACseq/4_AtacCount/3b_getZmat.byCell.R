
library(Seurat)
library(Signac)



args <- commandArgs(TRUE)

snATAC_object <- args[1]
cellType <- args[2]
outFile <- args[3]

##
seuset.atac <- readRDS(snATAC_object)

subseuset <- subset(seuset.atac, idents=cellType)
valMat <- subseuset@assays$RNA@data

z_score <- function(x){(x - mean(x)) /sd(x)}
zmat <- apply(valMat, 1, z_score)
zmat <- t(zmat)
zmat[is.na(zmat)] <- 0

df <- data.frame(geneID=rownames(zmat), zmat)
write.table(df, file=paste0('ZmatByCells/', cellType, '.snATAC.zscore.txt'), row.names=F, col.names=T, sep="\t", quote=F)



