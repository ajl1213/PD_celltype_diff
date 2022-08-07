library(ggplot2)
library(pheatmap)
library(RColorBrewer)


dir.create('Plots')


#======================================================================================================================================
# zscore heatmap (diff peaks)
#pval_thresh <- 0.05
adjPval_thresh <- 0.05
log2fc_thresh <- 0


##
valTable <- read.table(paste0('ValTable/H3K27ac.valTable.txt'), header=T)
#valTable.filter <- valTable[which(valTable$adjPval < adjPval_thresh & abs(valTable$logFC) > log2fc_thresh),]
valTable.filter <- valTable[which(valTable$adj.P.Val < adjPval_thresh & abs(valTable$logFC) > log2fc_thresh),]

print(paste0('down: ', sum(valTable.filter$logFC < 0)))
print(paste0('up: ', sum(valTable.filter$logFC >= 0)))

label <- c()
label[valTable.filter$logFC < 0] <- 'Down'
label[valTable.filter$logFC >= 0] <- 'Up'
valTable.filter$label <- label

##
#valMat <- read.table(paste0('AdjustedBulk/H3K27ac.cellAdjusted.qq.combat.txt'), header=T, row.names=1)
valMat <- read.table(paste0('AdjustedBulk/H3K27ac.cellAdjusted.combat.qq.txt'), header=T, row.names=1)
sampleList <- colnames(valMat)

get_zscore <- function(x){(x - mean(x)) /sd(x)}
zmat <- t(apply(valMat, 1, get_zscore))
zmat[is.na(zmat)] <- 0

zmat <- zmat[valTable.filter$peakID[valTable.filter$peakID %in% rownames(zmat)],]

minVal <- -1.5
maxVal <- 1.5
zmat[zmat < minVal]= minVal
zmat[zmat > maxVal]= maxVal

##
valColor=brewer.pal(11, 'RdBu')

pdf(paste0('Plots/H3K27ac.diffPeaks.zscore.pdf'))
zmat <- zmat[, sampleList]
pheatmap(zmat,
color= valColor,
cluster_rows=T, cluster_cols=T,
show_colnames=T, show_rownames=F,
cutree_rows=2, cutree_cols=1)
dev.off()

write.table(valTable.filter, 'ValTable/H3K27ac.valTable.filter.txt', col.names=T, row.names=F, sep='\t', quote=F)


#======================================================================================================================================
# atac cor
library(preprocessCore)

#
h3k27acFile <- read.table('/home/ajl1213/Projects/PD/SciAdv_Data/ChIP/02_Processing/CountMatrix/H3K27ac.readCount.denoised.txt', header=T, row.names=1)
NOSN_h3k27acCount <- h3k27acFile[, grep('NOSN', colnames(h3k27acFile))]
PDSN_h3k27acCount <- h3k27acFile[, grep('PDSN', colnames(h3k27acFile))]
#
atacFile <- read.table('/home/ajl1213/Projects/PD/data/snATAC/04_GetAtacCount/PeakCount/atacCount.byCellType.txt', header=T, row.names=1)

#
mean_NOSN <- rowSums(NOSN_h3k27acCount)
mean_PDSN <- rowSums(PDSN_h3k27acCount)
atacCount <- rowSums(atacFile)
df <- data.frame(atacCount, h3k27ac_NOSN=mean_NOSN[names(atacCount)], h3k27ac_PDSN=mean_PDSN[names(atacCount)])

#
qq_values <- normalize.quantiles.robust(as.matrix(df))
colnames(qq_values)=colnames(df)
rownames(qq_values)=rownames(df)
qq_values <- data.frame(qq_values)

print(dim(qq_values))

# H3K27ac NOSN vs snATAC
data1 <- data.frame(atacCount=qq_values$atacCount, h3k27ac_NOSN=qq_values$h3k27ac_NOSN)
min_thresh <- quantile(as.matrix(data1), 0.001)
minVals <- apply(data1, 1, min)
data1 <- data1[minVals > min_thresh, ]

g <- ggplot(data1) +
    geom_hex(aes(log2(data1$h3k27ac_NOSN+1), log2(data1$atacCount+1)), bins=350) +
    scale_fill_gradientn('Density', colours=rev(c('yellow','darkorange','red','navy','black')), limits=c(0,50)) +
    labs(x='log2 NOSN H3K27ac ChIP-seq reads', y='log2 merged snATAC-seq reads') +
    theme_classic() +
    theme(legend.position='none')

pdf('Plots/H3K27acNOSN.snATAC.cor.pdf')
plot(g)
dev.off()

summary(lm(log2(data1$h3k27ac_NOSN+1)~log2(data1$atacCount+1)))$adj.r.squared
cor(log2(data1$h3k27ac_NOSN+1), log2(data1$atacCount+1), method='spearman')
cor.test(log2(data1$h3k27ac_NOSN+1), log2(data1$atacCount+1), method='spearman')

# H3K27ac PDSN vs snATAC
data1 <- data.frame(atacCount=qq_values$atacCount, h3k27ac_PDSN=qq_values$h3k27ac_PDSN)
min_thresh <- quantile(as.matrix(data1), 0.001)
minVals <- apply(data1, 1, min)
data1 <- data1[minVals > min_thresh, ]

g <- ggplot(data1) +
    geom_hex(aes(log2(data1$h3k27ac_PDSN+1), log2(data1$atacCount+1)), bins=350) +
    scale_fill_gradientn('Density', colours=rev(c('yellow','darkorange','red','navy','black')), limits=c(0,50)) +
    labs(x='log2 PDSN H3K27ac ChIP-seq reads', y='log2 merged snATAC-seq reads') +
    theme_classic() +
    theme(legend.position='none')

pdf('Plots/H3K27acPDSN.snATAC.cor.pdf')
plot(g)
dev.off()

summary(lm(log2(data1$h3k27ac_PDSN+1)~log2(data1$atacCount+1)))$adj.r.squared
cor(log2(data1$h3k27ac_PDSN+1), log2(data1$atacCount+1), method='spearman')
cor.test(log2(data1$h3k27ac_PDSN+1), log2(data1$atacCount+1), method='spearman')




