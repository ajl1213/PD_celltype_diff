
library(ggplot2)


dir.create('Plots')


## Down
# factor list
cellTypeList <- c('DopaN','GabaN','Oligo','OPC','Ast','Micro','Endo','Peri')
#windowList <- c('15kb','100kb','200kb')
windowList <- c('100kb')
dataTypeList <- c('snDEG','bulkDEG')
list_factors <- c()
for (cellType in cellTypeList){
    for (winDist in windowList){
	for (dataType in dataTypeList){
	    label <- paste0(cellType, '.', dataType,'.', winDist)
	    list_factors <- c(list_factors, label)
	}
    }
}

# random values
randDF <- NULL
for (cellType in cellTypeList){
    for (winDist in windowList){

	data1 <- read.table(paste0('RandDegPermute/Down.', cellType, '.', winDist, '.RandDegPermuteRecord.txt'), header=T)

	df <- data.frame(DEG_type=data1$DEG_type ,cellType=cellType, winDist=winDist, fraction=data1$fraction)
	randDF <- rbind(randDF, df)
    }
}

randDF$dataLabel <- paste0(randDF$cellType,'.',randDF$DEG_type,'.',randDF$winDist)
randDF$dataLabel <- factor(randDF$dataLabel, levels=list_factors)
randDF$DEG_type <- factor(randDF$DEG_type, levels=dataTypeList)

# observed values
data1 <- read.table('RESULT/EnrichStatFinal.txt', header=T)
observedDF <- data1[which(data1$dataType=='Down'),]
observedDF$dataLabel <- paste0(observedDF$cellType, '.', observedDF$DEG_type,'.',observedDF$winDist)
observedDF$dataLabel <- factor(observedDF$dataLabel, levels=list_factors)
observedDF$DEG_type <- factor(observedDF$DEG_type, levels=dataTypeList)

#
p <- ggplot() +
    geom_point(aes(x=observedDF$dataLabel, y=observedDF$peakDegFrac, color=observedDF$DEG_type), shape=1, size=3) +
    geom_violin(aes(x=randDF$dataLabel, y=randDF$fraction, fill=randDF$DEG_type), trim=F, scale='width', width=0.9, adjust=4) +
    scale_color_manual(values=c('darkorange','navy'))+
    scale_fill_manual(values=c('bisque','lightblue'))+
    labs(x=NULL, y='Proportion of DEGs within genomic window') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

pdf('Plots/Down.violin.pdf', width=10, height=4, pointsize=3)
plot(p)
dev.off()



## Up
# factor list
cellTypeList <- c('DopaN','GabaN','Oligo','OPC','Ast','Micro','Endo','Peri')
#windowList <- c('15kb','100kb','200kb')
windowList <- c('100kb')
dataTypeList <- c('snDEG','bulkDEG')
list_factors <- c()
for (cellType in cellTypeList){
    for (winDist in windowList){
	for (dataType in dataTypeList){
	    label <- paste0(cellType, '.', dataType,'.', winDist)
	    list_factors <- c(list_factors, label)
	}
    }
}

# random values
randDF <- NULL
for (cellType in cellTypeList){
    for (winDist in windowList){

	data1 <- read.table(paste0('RandDegPermute/Up.', cellType, '.', winDist, '.RandDegPermuteRecord.txt'), header=T)

	df <- data.frame(DEG_type=data1$DEG_type ,cellType=cellType, winDist=winDist, fraction=data1$fraction)
	randDF <- rbind(randDF, df)
    }
}

randDF$dataLabel <- paste0(randDF$cellType,'.',randDF$DEG_type,'.',randDF$winDist)
randDF$dataLabel <- factor(randDF$dataLabel, levels=list_factors)
randDF$DEG_type <- factor(randDF$DEG_type, levels=dataTypeList)

# observed values
data1 <- read.table('RESULT/EnrichStatFinal.txt', header=T)
observedDF <- data1[which(data1$dataType=='Up'),]
observedDF$dataLabel <- paste0(observedDF$cellType, '.', observedDF$DEG_type,'.',observedDF$winDist)
observedDF$dataLabel <- factor(observedDF$dataLabel, levels=list_factors)
observedDF$DEG_type <- factor(observedDF$DEG_type, levels=dataTypeList)

#
p <- ggplot() +
    geom_point(aes(x=observedDF$dataLabel, y=observedDF$peakDegFrac, color=observedDF$DEG_type), shape=1, size=3) +
    geom_violin(aes(x=randDF$dataLabel, y=randDF$fraction, fill=randDF$DEG_type),trim=F, scale='width', width=0.9, adjust=4) +
    scale_color_manual(values=c('darkorange','navy'))+
    scale_fill_manual(values=c('bisque','lightblue'))+
    labs(x=NULL, y='Proportion of DEGs within genomic window') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

pdf('Plots/Up.violin.pdf', width=10, height=4, pointsize=3)
plot(p)
dev.off()






