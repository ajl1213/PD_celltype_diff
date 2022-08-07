

#======================================================================================================================
# Quantile normalization
library(preprocessCore)
data <- read.table("AdjustedBulk/RNA.cellAdjusted.iter_3.txt",head=T,row.names=1)
valMat=as.matrix(data[, 1:ncol(data)])

qq_values <- normalize.quantiles.robust(valMat)
colnames(qq_values)=colnames(valMat)

#
df.final=data.frame(geneID=row.names(data), qq_values)
write.table(df.final, file='AdjustedBulk/RNA.cellAdjusted.qq.txt', row.name=F, col.names=T, sep="\t", quote=F)

#======================================================================================================================
# Batch correction
library(sva)

#
batch1 <- c(
'X4689NOSN',
'X5628_1NOSN',
'X4870NOSN',
'X5006_1NOSN',
'X5244_1PDSN',
'X4831_1PDSN',
'X5649PDSN'
)

batch2 <- c(
'X5130_1NOSN',
'X5130_2NOSN',
'X4996_1NOSN',
'X5006_2NOSN',
'X5215_1PDSN',
'X5627_1PDSN',
'X5778PDSN',
'X5742PDSN',
'X5591PDSN'
)

#
diagnosis <- c()
diagnosis[grep('NOSN', colnames(data))] <- 'NOSN'
diagnosis[grep('PDSN', colnames(data))] <- 'PDSN'

metadata <- data.frame(row.names = colnames(data), diagnosis)

metadata$batch <- rownames(metadata)
metadata$batch[metadata$batch %in% batch1] <- 'b1'
metadata$batch[metadata$batch %in% batch2] <- 'b2'

#
valMat <- qq_values[,colnames(qq_values) %in% rownames(metadata)]
#model_combat <- model.matrix(~1, data=metadata)
model_combat <- model.matrix(~as.factor(metadata$diagnosis), data=metadata)
combat_data <- ComBat(dat=valMat, batch=metadata$batch, mod=model_combat, par.prior=TRUE, prior.plots=FALSE)
combat_data[combat_data<0]=0

#
#df.final=data.frame(geneID=row.names(data), combat_data)
#write.table(df.final, file='AdjustedBulk/RNA.cellAdjusted.qq.combat.txt', row.name=F, col.names=T, sep="\t", quote=F)

#======================================================================================================================
# Quantile normalization
valMat <- combat_data

qq_values <- normalize.quantiles.robust(valMat)
colnames(qq_values)=colnames(data)
rownames(qq_values)=rownames(data)

#
df.final=data.frame(geneID=row.names(data), qq_values)
write.table(df.final, file='AdjustedBulk/RNA.cellAdjusted.combat.qq.txt', row.name=F, col.names=T, sep="\t", quote=F)

#======================================================================================================================
# Zscore
valMat <- qq_values
get_zscore <- function(x){(x-mean(x))/sd(x)}
zmat <- t(apply(valMat, 1, get_zscore))
zmat[is.na(zmat)] <- 0

#
df.final=data.frame(geneID=row.names(data), zmat)
write.table(df.final, file='AdjustedBulk/RNA.cellAdjusted.combat.qq.zscore.txt', row.name=F, col.names=T, sep="\t", quote=F)




