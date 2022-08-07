
library(pheatmap)
library(RColorBrewer)


data1 <- read.table('matchedGWAS.celltype_asso.sample.txt', header=T)


df <- NULL
for (sampleID in unique(data1$SampleID)){
    valList <- c()
    for (cellType in unique(data1$CellType)){
	
	fc <- data1$foldEnrich[which(data1$SampleID==sampleID & data1$CellType==cellType)]
	valList <- c(valList, fc)
    }
    df <- rbind(df, valList)

}

rownames(df) <- unique(data1$SampleID)
colnames(df) <- unique(data1$CellType)

print(min(df))
print(max(df))


# heatmap
valColor <- brewer.pal(11, 'RdBu')

pdf('Plots/matchedGWAS.celltype.sample.heatmap.pdf', width=6, height=6, pointsize=3)
pheatmap(df, color= valColor,
cluster_rows=F, cluster_cols=F,
show_colnames=T, show_rownames=T,
breaks=seq(0.2, 1.5, length.out=11)
)
dev.off()


#
print(data1[which(data1$BinomPval < 0.05),])







