

library(pheatmap)
library(RColorBrewer)
library(viridis)


maxLogPval <- 5
data1 <- read.table('logPvalMat.txt', header=T, row.names=1)
data1[data1>maxLogPval] <- maxLogPval

print(summary(as.vector(as.matrix(data1))))

#valColor <- magma(15)
valColor <- brewer.pal(9, 'PuBu')


pdf('allPeak.LDSC.logPval.pdf')
pheatmap(data1, 
color=valColor,
cluster_rows=F, cluster_cols=F,
show_colnames=T, show_rownames=T, 
breaks=seq(0, maxLogPval, length.out=9)
)
dev.off()





