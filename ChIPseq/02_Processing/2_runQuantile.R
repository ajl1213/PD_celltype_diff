library(preprocessCore)



#=========================================================================================================
## Quantile normalization

# H3K27ac
data <- read.table("CountMatrix/H3K27ac.readCount.denoised.txt",head=T,row.names=1)
peakList <- read.table('/home/ajl1213/Projects/PD/data/snATAC/04_GetAtacCount/CellLabel/AtacPeak.CellLabel.txt', header=T)$peakID

data <- data[peakList,]
data=as.matrix(data)

qq_values <- normalize.quantiles.robust(data)
colnames(qq_values)=colnames(data)
qq_values=data.frame(peakID=row.names(data), qq_values)

write.table(qq_values, file='CountMatrix/H3K27ac.readCount.denoised.qq.txt', row.name=F, col.names=T, sep="\t", quote=F)





