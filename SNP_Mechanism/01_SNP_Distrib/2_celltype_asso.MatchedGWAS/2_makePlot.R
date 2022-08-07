
library(ggplot2)
library(RColorBrewer)


dir.create('Plots')

colList <- c(brewer.pal(9, 'PuBu')[7], 'darkorange')

#
data1 <- read.table('var_cRE_asso.v2.txt', header=T)

data1$PeakType <- factor(data1$PeakType, levels=c('Down','Up','nonDys'))
label <- c()
label[grep('NOSN', data1$SampleID)] <- 'NOSN'
label[grep('PDSN', data1$SampleID)] <- 'PDSN'
data1$label <- label
data1$label <- factor(data1$label, levels=c('NOSN','PDSN'))


#
p1 <- ggplot(data1, aes(x=data1$PeakType, y=data1$log2Enrich, fill=data1$label)) +
    geom_boxplot(outlier.shape=NA) +
    scale_fill_manual(values=colList) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(axis.title.x=element_blank())

pdf('Plots/var_cRE_asso.box.pdf')
plot(p1)
dev.off()



