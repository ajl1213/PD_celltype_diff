library(preprocessCore)



#=========================================================================================================
## Quantile normalization

data <- read.table("CountMatrix/RNA.rawCount.txt",head=T,row.names=1)
valMat=as.matrix(data[, 2:ncol(data)])

qq_values <- normalize.quantiles.robust(valMat)
colnames(qq_values)=colnames(valMat)


df.final=data.frame(ensembleID=row.names(data), geneID=data$geneID, qq_values)
write.table(df.final, file='CountMatrix/RNA.rawCount.qq.txt', row.name=F, col.names=T, sep="\t", quote=F)

df.final=data.frame(ensembleID=row.names(data), geneID=data$geneID, log2(qq_values+1))
write.table(df.final, file='CountMatrix/RNA.rawCount.qq.log2.txt', row.name=F, col.names=T, sep="\t", quote=F)



