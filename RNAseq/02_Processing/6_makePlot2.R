library(ggplot2)
library(RColorBrewer)


nIter <- 3
sampleList <- c(
'X4870NOSN',
'X4689NOSN',
'X5628NOSN',
'X5130_1NOSN',
'X5130_2NOSN',
'X4996NOSN',
'X5006_1NOSN',
'X5006_2NOSN',
'X5244PDSN',
'X4831PDSN',
'X5215PDSN',
'X5627PDSN',
'X5778PDSN',
'X5742PDSN',
'X5649PDSN',
'X5591PDSN'
)

cellTypeList <- c('DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri')
celltype_colors <- c(brewer.pal(8, 'Dark2'))
#===========================================================================================================================================================
## Barplot

# initial composition
initialComp <- read.table('MarkerGeneRatio/marker_composition.txt', header=T)
initialComp$sampleID <- factor(initialComp$sampleID, levels=rev(sampleList))
initialComp$cellType <- factor(initialComp$cellType, levels=rev(cellTypeList))

g1 <- ggplot(data=initialComp, aes(x=sampleID, y=ratio, fill=cellType)) +
    geom_bar(stat='identity') + 
    scale_fill_manual(values=rev(celltype_colors)) + 
    coord_flip() + 
    labs(x=NULL, y='count') + 
    theme_classic()

pdf('Plots/initial_composition.RNA.pdf')
plot(g1)
dev.off()

# final composition
finalComp <- read.table(paste0('AdjustedBulk/marker_composition.iter_', nIter, '.txt'), header=T)
finalComp$sampleID <- factor(finalComp$sampleID, levels=rev(sampleList))
finalComp$cellType <- factor(finalComp$cellType, levels=rev(cellTypeList))

g1 <- ggplot(data=finalComp, aes(x=sampleID, y=ratio, fill=cellType)) +
    geom_bar(stat='identity') + 
    scale_fill_manual(values=rev(celltype_colors)) + 
    coord_flip() + 
    labs(x=NULL, y='count') + 
    theme_classic()

pdf('Plots/final_composition.RNA.pdf')
plot(g1)
dev.off()


#===========================================================================================================================================================
## Iteration plot
#
avgRatioDF <- read.table('MarkerGeneRatio/avg_marker_composition.txt', header=T)
# 
initialComp <- read.table('MarkerGeneRatio/marker_composition.txt', header=T)
initialComp$label <- 'initial'
#
df <- initialComp
for (i in 1:nIter){
    iterComp <- read.table(paste0('AdjustedBulk/marker_composition.iter_', i, '.txt'), header=T)
    iterComp$label <- paste0('iter_', i)
    df <- rbind(df, iterComp)
}

#
diagnosis <- c()
diagnosis[grep('NOSN', df$sampleID)] <- 'NOSN'
diagnosis[grep('PDSN', df$sampleID)] <- 'PDSN'
df$diagnosis <- diagnosis

df$dataType <- paste0(df$sampleID,':',df$cellType)

#
df$label <- factor(df$label, levels=c('initial','iter_1','iter_2','iter_3'))
df$diagnosis <- factor(df$diagnosis, levels=c('NOSN','PDSN'))
#
g1 <- ggplot(df, aes(x=label, y=ratio, group=dataType)) +
    geom_line(aes(color=diagnosis), size=0.5) +
    geom_point(aes(color=diagnosis), shape=1, size=1) +
    scale_color_manual(values=c('green','red')) +
    geom_hline(yintercept=avgRatioDF$avgRatio, color=celltype_colors, size=0.5, linetype='dashed') +
    theme_classic() +
    labs(x='Iteration',y='composition') +

pdf('Plots/IterativeComposition.pdf', height=4, width=4, pointsize=3)
plot(g1)
dev.off()


#===========================================================================================================================================================
## hclust

rawCountDF <- read.table('/home/ajl1213/Projects/PD/SciAdv_Data/RNA/01_Mapping/CountMatrix/RNA.rawCount.qq.txt', header=T, row.names=1)
rawCountDF <- rawCountDF[,2:ncol(rawCountDF)]
finalCountDF <- read.table('AdjustedBulk/RNA.cellAdjusted.combat.qq.txt', header=T, row.names=1)

#
get_zscore <- function(x){(x-mean(x))/sd(x)}
rawZmat <- as.matrix(t(apply(rawCountDF, 1, get_zscore)))
finalZmat <- as.matrix(t(apply(finalCountDF, 1, get_zscore)))
rawZmat[is.na(rawZmat)] <- 0
finalZmat[is.na(finalZmat)] <- 0

#
hr <- hclust(as.dist(1-cor(rawZmat, method='pearson')), method='complete')
pdf('Plots/Raw.hclust.pdf')
plot(hr)
dev.off()

#
hr <- hclust(as.dist(1-cor(finalZmat, method='pearson')), method='complete')
pdf('Plots/Final.hclust.pdf')
plot(hr)
dev.off()







