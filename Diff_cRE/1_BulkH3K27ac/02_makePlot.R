library(pheatmap)
library(RColorBrewer)


dir.create('Plots')

for (i in c('Down','Up')){
    print(i)

    #========================================================================================================================================================================
    # order
    data1 <- read.table(paste0('InputData/',i,'.funcCell.txt'), header=T, row.names=1)
    d <- dist(data1, method='minkowski')
    hr <- hclust(d, method='average')
    hClustOrder <- as.character(rownames(data1)[hr$order])

    #========================================================================================================================================================================
    # functional celltype
    valMat <- data1[hClustOrder,]
    print(colSums(valMat)/nrow(valMat)*100)

    valColor <- c('white', 'black')

    pdf(paste0('Plots/H3K27ac.',i,'.funcCell.heatmap.pdf'))
    pheatmap(valMat,
	color= valColor,
	cluster_rows=F, cluster_cols=F,
	show_colnames=T, show_rownames=F,# labels_row=gene_label,
	cutree_rows=1, cutree_cols=1
    )
    dev.off()

    #========================================================================================================================================================================
    # log2fc heatmap
    minVal <- -2
    maxVal <- 2
    countMat <- read.table('/home/ajl1213/Projects/PD/SciAdv_Data/ChIP/02_Processing/AdjustedBulk/H3K27ac.cellAdjusted.qq.combat.txt', header=T, row.names=1)
    get_log2fc <- function(x){log2((x+1)/(mean(x)+1))}
    log2fcMat <- t(apply(countMat, 1, get_log2fc))
    log2fcMat[log2fcMat < minVal] <- minVal
    log2fcMat[log2fcMat > maxVal] <- maxVal

    degLog2fcMat <- log2fcMat[hClustOrder, ]

    # metadata
    metadata <- read.table('/home/ajl1213/Projects/PD/SciAdv_Data/ChIP/01_Mapping/Metadata/metadata.final.txt', header=T, row.names=1)
    valColor <- colorRampPalette(c('blue','white','red'))(15)

    #
    annotation_col <- data.frame(row.names=rownames(metadata), Diagnosis=metadata$Diagnosis, Gender=metadata$Gender, AGE=metadata$AGE, BraakCombo=metadata$BraakCombo, MMSE=metadata$MMSE, DRS=metadata$DRS)

    anno_colors <- list(
    BraakCombo=c('white','red'),
    MMSE=c('red', 'white'),
    DRS=c('red','white'),
    Gender=c(Male=brewer.pal(8,'Accent')[5], Female=brewer.pal(8,'Accent')[3]),
    AGE=c('white','forestgreen'),
    Diagnosis=c(Normal='bisque',Parkinson='orangered')
    )

    pdf(paste0('Plots/H3K27ac.', i, '.log2fc.heatmap.pdf'))
    pheatmap(degLog2fcMat,
	color= valColor,
	annotation_col=annotation_col,
	annotation_colors=anno_colors,
	cluster_rows=F, cluster_cols=F,
	show_colnames=T, show_rownames=F,
	cutree_rows=1, cutree_cols=1
    )
    dev.off()

}



