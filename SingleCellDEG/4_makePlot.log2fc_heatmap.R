library(pheatmap)
library(RColorBrewer)


##===================================== make log2fc heatmap
data1 <- read.table('PDGenes.valMat.txt', header=T, row.names=1)


color_range <- colorRampPalette(c('blue','white','red'))(15)

pdf(paste0('Plots/PDGenes.log2fc.heatmap.pdf'))
pheatmap(data1, color= color_range, 
cluster_rows=T, cluster_cols=F,
show_colnames=T, show_rownames=T,
breaks=seq(-1, 1, length.out=15)
)
dev.off()



#
data2 <- read.table('PDGenes.valList.txt', header=T)

sigDEGs <- data2[which(abs(data2$Log2FC) > 0.2 & data2$AdjPval < 0.05 & data2$pct1 > 0.2 & max(data2$pct1, data2$pct2) > 0.2), ]


print(sigDEGs)




