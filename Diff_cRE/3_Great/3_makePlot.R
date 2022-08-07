
library(pheatmap)
library(RColorBrewer)


dir.create('Plots')

maxVal <- 15


#
col_range <- brewer.pal(10, 'YlOrBr')
downPlotDF <- read.table('Down.plotMat.txt', header=T, row.names=1)
downPlotDF[downPlotDF > maxVal] <- maxVal


pdf('Plots/Down.GO.pdf', width=10, height=8, pointsize=3)
pheatmap(downPlotDF,
color= col_range,
cluster_rows=F, cluster_cols=F,
show_colnames=T, show_rownames=T,
cutree_rows=1, cutree_cols=1
)

dev.off()

#
col_range <- brewer.pal(10, 'YlGnBu')

upPlotDF <- read.table('Up.plotMat.txt', header=T, row.names=1)
upPlotDF[upPlotDF > maxVal] <- maxVal

pdf('Plots/Up.GO.pdf', width=10, height=8, pointsize=3)
pheatmap(upPlotDF,
color= col_range,
cluster_rows=F, cluster_cols=F,
show_colnames=T, show_rownames=T,
cutree_rows=1, cutree_cols=1
)
dev.off()




