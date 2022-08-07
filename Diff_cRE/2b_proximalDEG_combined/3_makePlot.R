
library(ggplot2)


dir.create('Plots')


## Down
# factor list
cellTypeList <- c('DopaN','GabaN','Oligo','OPC','Ast','Micro','Endo','Peri')
windowList <- c('100kb')

# random values
randDF <- NULL
for (cellType in cellTypeList){
    for (winDist in windowList){

	data1 <- read.table(paste0('RandPeakPermute/Down.', cellType, '.', winDist, '.RandPeakPermuteRecord.txt'), header=T)

	df <- data.frame(cellType=cellType, fraction=data1$fraction)
	randDF <- rbind(randDF, df)
    }
}

# observed values
data1 <- read.table('RESULT/EnrichStatFinal.txt', header=T)
observedDF <- data1[which(data1$dataType=='Down'),]

randDF$cellType <- factor(randDF$cellType, levels=cellTypeList)
observedDF$cellType <- factor(observedDF$cellType, levels=cellTypeList)

#
p <- ggplot() +
    geom_point(aes(x=observedDF$cellType, y=observedDF$peakDegFrac), color='navy', shape=1, size=3) +
    geom_violin(aes(x=randDF$cellType, y=randDF$fraction), fill='bisque', trim=F, scale='width', width=0.9, adjust=4) +
    labs(x=NULL, y='Proportion of DEGs within genomic window') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

pdf('Plots/Down.violin.pdf', width=10, height=4, pointsize=3)
plot(p)
dev.off()



## Up
# factor list
cellTypeList <- c('DopaN','GabaN','Oligo','OPC','Ast','Micro','Endo','Peri')
windowList <- c('100kb')

# random values
randDF <- NULL
for (cellType in cellTypeList){
    for (winDist in windowList){

	data1 <- read.table(paste0('RandPeakPermute/Up.', cellType, '.', winDist, '.RandPeakPermuteRecord.txt'), header=T)

	df <- data.frame(cellType=cellType, fraction=data1$fraction)
	randDF <- rbind(randDF, df)
    }
}

# observed values
data1 <- read.table('RESULT/EnrichStatFinal.txt', header=T)
observedDF <- data1[which(data1$dataType=='Up'),]

randDF$cellType <- factor(randDF$cellType, levels=cellTypeList)
observedDF$cellType <- factor(observedDF$cellType, levels=cellTypeList)

#
p <- ggplot() +
    geom_point(aes(x=observedDF$cellType, y=observedDF$peakDegFrac), color='darkorange', shape=1, size=3) +
    geom_violin(aes(x=randDF$cellType, y=randDF$fraction), fill='bisque', trim=F, scale='width', width=0.9, adjust=4) +
    labs(x=NULL, y='Proportion of DEGs within genomic window') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

pdf('Plots/Up.violin.pdf', width=10, height=4, pointsize=3)
plot(p)
dev.off()





