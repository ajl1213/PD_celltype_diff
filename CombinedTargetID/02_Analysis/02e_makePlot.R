

library(pheatmap)
library(RColorBrewer)


valColor <- brewer.pal(9, 'PuBu')
maxLogPval <- 5

data1 <- read.table('GOTABLESORTED/MP.celltype.txt', header=T, row.names=1, sep='\t')
data1[data1>maxLogPval] <- maxLogPval


pdf('Plots/MP.heatmap.pdf')
pheatmap(data1,
annotation_colors = setAnnoColor, color=valColor,
show_colnames=T, show_rownames=T,
cluster_rows = F, cluster_cols = F,
breaks=seq(0, maxLogPval, length.out=9)
)
dev.off()




