

library(ggplot2)
library(pheatmap)
library(RColorBrewer)


dir.create('Plots')

#pval_thresh <- 0.05
adjPval_thresh <- 0.05
log2fc_thresh <- 1


##
valTable <- read.table(paste0('ValTable/RNA.valTable.txt'), header=T)
#valTable.filter <- valTable[which(valTable$adjPval < adjPval_thresh & abs(valTable$logFC) > log2fc_thresh),]
valTable.filter <- valTable[which(valTable$adj.P.Val < adjPval_thresh & abs(valTable$logFC) > log2fc_thresh),]

print(paste0('down: ', sum(valTable.filter$logFC < 0)))
print(paste0('up: ', sum(valTable.filter$logFC >= 0)))

label <- c() 
label[valTable.filter$logFC < 0] <- 'Down'
label[valTable.filter$logFC >= 0] <- 'Up'
valTable.filter$label <- label

##
valMat <- read.table(paste0('AdjustedBulk/RNA.cellAdjusted.combat.qq.txt'), header=T, row.names=1)
sampleList <- colnames(valMat)

get_zscore <- function(x){(x - mean(x)) /sd(x)}
zmat <- t(apply(valMat, 1, get_zscore))
zmat[is.na(zmat)] <- 0

zmat <- zmat[valTable.filter$geneID[valTable.filter$geneID %in% rownames(zmat)],]

minVal <- -2
maxVal <- 2
zmat[zmat < minVal]= minVal
zmat[zmat > maxVal]= maxVal

##
valColor=brewer.pal(11, 'RdBu')

pdf(paste0('Plots/RNA.diffPeaks.zscore.pdf'))
zmat <- zmat[, sampleList]
pheatmap(zmat,
color= valColor,
cluster_rows=T, cluster_cols=T,
show_colnames=T, show_rownames=F,
cutree_rows=2, cutree_cols=1)
dev.off()

write.table(valTable.filter, 'ValTable/RNA.valTable.filter.txt', col.names=T, row.names=F, sep='\t', quote=F)





