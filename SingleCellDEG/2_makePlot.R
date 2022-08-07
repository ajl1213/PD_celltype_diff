library(pheatmap)
library(RColorBrewer)



dir.create('Plots')

knownGeneList <- c(
'SNCA', 'UCHL1', 'PINK1', 'PARK7','LRRK2',
'ATP13A2', 'FBXO7', 'VPS35', 'DNAJC6', 'SYNJ1',
'CHCHD2', 'VPS13C', 'GBA', 'MAPT', 'GAK',
'SMPD1', 'SCARB2', 'SLC17A5','ATP6V0A1','CTSB'
)

for (dataType in c('Down','Up')){
    print(dataType)

    #=============================== get order (label)
    orderFile <- read.table(paste0('DegOrdered.', dataType,'.txt'), header=T)
    orderList <- orderFile$geneID
    print(intersect(orderList, knownGeneList))

    #=============================== get order (label)
    data1 <- read.table(paste0('DegMatrix/DegMatrix.', dataType,'.txt'), header=T, row.names=1)
    labelMat.order <- as.matrix(data1[orderList,])

    # gene label
    gene_label <- NULL
    gene_label[rownames(labelMat.order) %in% intersect(orderList, knownGeneList)] <- rownames(labelMat.order)[rownames(labelMat.order) %in% intersect(orderList, knownGeneList)]
    gene_label[is.na(gene_label)] <- ''

    valColor <- c('white', brewer.pal(8, 'Dark2'))

    pdf(paste0('Plots/snRNA_DEG.label.', dataType,'.pdf'))
    pheatmap(labelMat.order, color= valColor, 
    cluster_rows=F, cluster_cols=F,
    show_colnames=T, show_rownames=T, labels_row=gene_label
    )
    dev.off()

    #=============================== make log2fc heatmap
    data2 <- read.table(paste0('DegMatrix/DegMatrix.log2fc.', dataType, '.txt'), header=T, row.names=1)
    valMat.order <- as.matrix(data2[orderList,])

    valColor=colorRampPalette(c('blue','white','red'))(15)

    minVal <- -1
    maxVal <- 1

    valMat.order[valMat.order > maxVal] <- maxVal
    valMat.order[valMat.order < minVal] <- minVal

    pdf(paste0('Plots/snRNA_DEG.log2fc.', dataType,'.pdf'))
    pheatmap(valMat.order, color= valColor, 
    cluster_rows=F, cluster_cols=F,
    show_colnames=T, show_rownames=F)
    dev.off()
}





