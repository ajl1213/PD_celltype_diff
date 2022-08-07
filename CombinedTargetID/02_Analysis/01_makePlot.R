library(ggplot2)
library(ggrepel)
library(RColorBrewer)


dir.create('Plots')

inputDIR <- '/home/ajl1213/Projects/PD/SciAdv_Anal/CombinedTargetID/01_FindTarget/RESULT'
cellTypeList <- c('DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri')

minAbcScore <- 10
maxAbcScore <- 50

for (cellType in cellTypeList){
    print(cellType)

    #
    data1 <- read.table(paste0(inputDIR, '/', cellType, '.selectList.txt'), header=T)

    #
    data1$abcScore[data1$abcScore > maxAbcScore] <- maxAbcScore
    data1$creType <- factor(data1$creType, levels=c('DownPeak','UpPeak','GwasPeak'))

    # 
    nonSigDF <- data1[which(data1$abcScore <= minAbcScore),]
    SigDF <- data1[which(data1$abcScore > minAbcScore),]
    print(table(SigDF$creType))
    print(table(SigDF$geneLabel, SigDF$creType))
    print('DownPeakTarget / DownDEG')
    print(unique(SigDF$geneID[which(SigDF$creType=='DownPeak' & SigDF$geneLabel=='downDEG')]))
    print('UpPeakTarget / UpDEG')
    print(unique(SigDF$geneID[which(SigDF$creType=='UpPeak' & SigDF$geneLabel=='upDEG')]))
    print('GwasPeakTarget / DownDEG')
    print(unique(SigDF$geneID[which(SigDF$creType=='GwasPeak' & SigDF$geneLabel=='downDEG')]))
    print('GwasPeakTarget / nonDEG')
    print(unique(SigDF$geneID[which(SigDF$creType=='GwasPeak' & SigDF$geneLabel=='nonDEG')]))
    print('GwasPeakTarget / UpDEG')
    print(unique(SigDF$geneID[which(SigDF$creType=='GwasPeak' & SigDF$geneLabel=='upDEG')]))

    # 
    labelDF <- SigDF[which(SigDF$abcScore > minAbcScore),]
    top10Genes <- labelDF$geneID[order(labelDF$abcScore, decreasing=T)][1:10]
    PdGenes <- labelDF$geneID[which(labelDF$creType=='GwasPeak' | labelDF$PdGeneLabel=='knownPdGene')]
    labelGeneList <- c(top10Genes, PdGenes)
    labelDF <- labelDF[labelDF$geneID %in% labelGeneList,]

    #
    g1 <- ggplot() + 
	geom_point(aes(x=nonSigDF$abcScore, y=log2(nonSigDF$contactVal), shape=nonSigDF$creType), color='grey', size=0.5) +
	geom_point(aes(x=SigDF$abcScore, y=log2(SigDF$contactVal), color=SigDF$creType, shape=SigDF$creType), size=1.5) +
	geom_text_repel(aes(x=labelDF$abcScore, y=log2(labelDF$contactVal), label=labelDF$geneID, color=labelDF$creType), size=2) +
	geom_vline(xintercept=minAbcScore, linetype='dashed') + 
	scale_color_manual(values=c('darkgreen','darkorange','darkblue')) +
	scale_shape_manual(values=c(1, 1, 4)) +
	xlim(0, maxAbcScore) +
	labs(title=cellType, x='ABC score by dysregulated cREs)', y='log2(Cumulative contact value (cRE X interaction) by dysregulated cREs)') +
	theme_classic()

    ##
    pdf(paste0('Plots/', cellType,'.ContactScorePlot.pdf'), width=7, height=4, pointsize=3)
    plot(g1)
    dev.off()
}




