library(ggplot2)


dir.create('Plots')


data1 <- read.table('GOTABLESORTED/Merged.MP.txt', header=T, sep='\t')
data1$termID <- factor(data1$termID, levels=data1$termID)



## barplot
g1 <- ggplot(data=data1, aes(x=data1$termID, y=data1$nGene))+ 
    geom_bar(stat='identity') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

pdf('Plots/MP.barplot.pdf', width=10, height=7, pointsize=3)
plot(g1)
dev.off() 



## heatmap
library(pheatmap)
library(RColorBrewer)
library(viridis)


valColor <- brewer.pal(9, 'PuBu')
maxLogPval <- 5
nMax <- 10

#data1 <- read.table('GOTABLESORTED/MP.celltype.txt', header=T, row.names=1, sep='\t')
data1 <- read.table('GOTABLESORTED/Merged.MP.txt', header=T, row.names=1, sep='\t')

data1_filtered <- data1[1:nMax,]
data1_filtered$logP <- -log10(data1_filtered$P.value)
data1_filtered$logP[data1_filtered$logP > maxLogPval] <- maxLogPval

df <- data.frame(logP = data1_filtered$logP)
rownames(df) <- rownames(data1_filtered)
print(df)

pdf('Plots/MergedMP.heatmap.pdf')
pheatmap(df, color=valColor,
show_colnames=T, show_rownames=T,
cluster_rows = F, cluster_cols = F,
breaks=seq(0, maxLogPval, length.out=9)
)
dev.off()


##
color_range <- colorRampPalette(c('dodgerblue3', 'lightblue'))(nMax)
pdf('Plots/MergedMP.barplot.pdf')
barplot(df$logP, names.arg = rownames(df),  col=color_range, las=2)
abline(h=-log10(0.05), col='darkgrey', size=2, lty=2)
dev.off()



