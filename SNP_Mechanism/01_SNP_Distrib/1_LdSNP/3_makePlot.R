

library(ggplot2)
library(RColorBrewer)


dir.create('Plots')


data1 <- read.table('PD_GWAS.match.txt', header=T)


##===========================================================================================================================
# Barplot
ctrlDF <- data1[grep('NOSN', data1$SampleID), ]
pdDF <- data1[grep('PDSN', data1$SampleID), ]

ctrlDF.order <- ctrlDF[order(ctrlDF$GWAS_frac, decreasing=T),]
pdDF.order <- pdDF[order(pdDF$GWAS_frac, decreasing=T),]

#
label <- c()
label[grep('NOSN', data1$SampleID)] <- 'NOSN'
label[grep('PDSN', data1$SampleID)] <- 'PDSN'

#
mergedDF <- rbind(ctrlDF.order, pdDF.order)
mergedDF$label <- label
mergedDF$SampleID <- factor(mergedDF$SampleID, levels=mergedDF$SampleID)
mergedDF$label <- factor(mergedDF$label, levels=c('NOSN','PDSN'))

#
p1 <- ggplot(data=mergedDF, aes(x=SampleID, y=GWAS_frac, fill=label)) +
    geom_bar(stat="identity", width=1, color='black') +
    scale_fill_manual(values=c( brewer.pal(9, 'PuBu')[7], 'darkorange')) +
    geom_hline(yintercept=c(mean(ctrlDF$GWAS_frac), mean(pdDF$GWAS_frac)), linetype='dashed', color='black')+
    ylim(0, 6) +
    labs(title=NULL, x=NULL, y='Matched with PD GWAS-SNPs / total variants called * 1000') +
    theme_classic() +
    theme(axis.text.x=element_text(angle=90, hjust=1))

pdf('Plots/PD_GWAS_match.bar.pdf')
plot(p1)
dev.off()

#
print(mean(ctrlDF$GWAS_frac))
print(mean(pdDF$GWAS_frac))

#
t.test(ctrlDF$GWAS_frac, pdDF$GWAS_frac)
t.test(ctrlDF$GWAS_frac, pdDF$GWAS_frac, alternative='less')

wilcox.test(ctrlDF$GWAS_frac, pdDF$GWAS_frac)
wilcox.test(ctrlDF$GWAS_frac, pdDF$GWAS_frac, alternative='less')

##===========================================================================================================================
# Boxplot
col_list <- c(brewer.pal(8, 'Set2')[1], brewer.pal(8, 'Set2')[2])

p1 <- ggplot() + 
    geom_boxplot(aes(mergedDF$label, mergedDF$GWAS_frac), outlier.colour=NA) + 
    geom_jitter(aes(mergedDF$label, mergedDF$GWAS_frac, col=mergedDF$label), shape=1, position=position_jitter(0.2), cex=2) +
    scale_color_manual(values=col_list) +
    labs(y='Matched with PD GWAS-SNPs / total variants called * 1000', x='') +
    theme_classic()

pdf('Plots/PD_GWAS_match.box.pdf')
plot(p1)
dev.off()





