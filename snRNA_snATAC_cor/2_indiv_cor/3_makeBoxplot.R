
library(ggplot2)
library(RColorBrewer)



cellTypeList <- c('DopaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri')

data1 <- read.table('Merged.cor.label.txt', header=T)
data1$label <- factor(data1$label, levels=c('match','unmatch'))
data1$cellType <- factor(data1$cellType, levels=cellTypeList)
data1$atacLibType <- factor(data1$atacLibType, levels=c('singleLib','pooledLib'))



## by celltype
p1 <- ggplot() +
    geom_boxplot(data=data1, aes(x=cellType, y=pcc, fill=label), outlier.colour=NA) +
    geom_hline(yintercept=0, linetype='dashed', color='black') +
    scale_fill_manual(values=c('darkorange','grey')) +
    labs(x=NULL, y='PCC with snRNA-seq and snATAC-seq samples') +
    theme_classic()

pdf('Plots/snRNA_snATAC.cor.byCellType.box.pdf')
plot(p1)
dev.off()

#
for (cellType in cellTypeList){
    print(cellType)

    tmpDF <- data1[which(data1$cellType==cellType), ]
    a <- wilcox.test(tmpDF$pcc[which(tmpDF$label=='match')], tmpDF$pcc[which(tmpDF$label=='unmatch')])
    print(a$p.value)
}


## merged
p1 <- ggplot() +
    geom_boxplot(data=data1, aes(x=label, y=pcc, fill=label), outlier.colour=NA) +
    geom_hline(yintercept=0, linetype='dashed', color='black') +
    scale_fill_manual(values=c('darkorange','grey')) +
    labs(x=NULL, y='PCC with snRNA-seq and snATAC-seq samples') +
    theme_classic()

pdf('Plots/snRNA_snATAC.cor.merged.box.pdf')
plot(p1)
dev.off()

#
a <- wilcox.test(data1$pcc[which(data1$label=='match')], data1$pcc[which(data1$label=='unmatch')])
print(a$p.value)



## Pooled Libs
p1 <- ggplot() +
    geom_boxplot(data=data1, aes(x=atacLibType, y=pcc, fill=label), outlier.colour=NA) +
    geom_hline(yintercept=0, linetype='dashed', color='black') +
    scale_fill_manual(values=c('darkorange','grey')) +
    labs(x=NULL, y='PCC with snRNA-seq and snATAC-seq samples') +
    theme_classic()

pdf('Plots/snRNA_snATAC.cor.byLibType.box.pdf')
plot(p1)
dev.off()

#
a <- wilcox.test(data1$pcc[which(data1$label=='match')], data1$pcc[which(data1$label=='unmatch')])
print(a$p.value)


for (atacLibType in c('singleLib','pooledLib')){
    print(atacLibType)

    tmpDF <- data1[which(data1$atacLibType==atacLibType), ]
    a <- wilcox.test(tmpDF$pcc[which(tmpDF$label=='match')], tmpDF$pcc[which(tmpDF$label=='unmatch')])
    print(a$p.value)
}




