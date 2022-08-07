
library(pheatmap)
library(RColorBrewer)


dir.create('Plots')

cellTypeList <- c('DopaN','GabaN','Oligo','OPC','Ast','Micro','Endo','Peri')
maxPval <- 10


# DOWN
col_range <- brewer.pal(9, 'YlGnBu')
data1 <- read.table(paste0('PvalMat/Down.logPval.txt'), header=T)

dataLabel <- paste0(data1$cellType, '.', data1$termID)
valMat <- data1[, 4:ncol(data1)]
rownames(valMat) <- dataLabel


pdf('Plots/Down.combinedGO.pdf')
pheatmap(valMat,
color= col_range,
cluster_rows=F, cluster_cols=F,
show_colnames=T, show_rownames=T,
cutree_rows=1, cutree_cols=1,
breaks=seq(0, maxPval, length.out=9),
display_numbers = T, number_format='%.1f'
)
dev.off()


# UP
col_range <- brewer.pal(9, 'YlOrBr')
data1 <- read.table(paste0('PvalMat/Up.logPval.txt'), header=T)

dataLabel <- paste0(data1$cellType, '.', data1$termID)
valMat <- data1[, 4:ncol(data1)]
rownames(valMat) <- dataLabel


pdf('Plots/Up.combinedGO.pdf')
pheatmap(valMat,
color= col_range,
cluster_rows=F, cluster_cols=F,
show_colnames=T, show_rownames=T,
cutree_rows=1, cutree_cols=1,
breaks=seq(0, maxPval, length.out=9),
display_numbers = T, number_format='%.1f'
)
dev.off()






