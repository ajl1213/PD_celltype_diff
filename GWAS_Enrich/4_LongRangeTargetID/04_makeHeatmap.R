
library(pheatmap)
library(RColorBrewer)


col_range <- colorRampPalette(c('white','yellow','darkorange','red','darkred'))(20)
celltype_colors <- brewer.pal(8, 'Dark2') 
maxVal <- 60

##
data1 <- read.table('AbcScoreMat.filtered.ordered.txt', header=T, row.names=1)


# row annotation
annoRow <- data.frame(maxCell=data1$maxCell, nearbyTarget=data1$NearbyTarget, longRangeTarget=data1$LongRangeTarget)
rownames(annoRow) <- rownames(data1)

setAnnoColor <- list(
maxCell=c(
DopaN=celltype_colors[1], GabaN=celltype_colors[2],
Oligo=celltype_colors[3], OPC=celltype_colors[4],
Ast=celltype_colors[5], Micro=celltype_colors[6],
Endo=celltype_colors[7], Peri=celltype_colors[8]
),
nearbyTarget=c(
Target='black', nonTarget='white'
),
longRangeTarget=c(
Target='black', nonTarget='white'
)
)


#
valMat <- data1[, 4:ncol(data1)]
summary(as.vector(as.matrix(valMat)))


pdf('AbcScoreMat.heatmap.pdf')
pheatmap(valMat,
color= col_range, 
annotation_row=annoRow, annotation_colors=setAnnoColor,
cluster_rows=F, cluster_cols=F,
show_colnames=T, show_rownames=T,
breaks=seq(0, maxVal, length.out=20)
)
dev.off()

