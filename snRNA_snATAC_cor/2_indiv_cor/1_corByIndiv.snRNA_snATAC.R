

library(pheatmap)
library(RColorBrewer)


dir.create('Plots')
dir.create('CorMat')

cellTypeList <- c('DopaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri')

rna_DIR <- '/home/ajl1213/Projects/PD/SciAdv_Data/snRNAseq/03_GetRnaCount/CountByIndiv'
atac_DIR <- '/home/ajl1213/Projects/PD/SciAdv_Data/snATACseq/4_AtacCount/CountByIndiv'

for (cellType in cellTypeList){
    print(cellType)

    rna_file <- read.table(paste0(rna_DIR,'/', cellType, '.RawCount.qq.txt'), header=T, row.names=1)
    atac_file <- read.table(paste0(atac_DIR,'/', cellType, '.snATAC.RawCount.qq.txt'), header=T, row.names=1)

    rna_file <- rna_file[, grep('X', colnames(rna_file))]
    atac_file <- atac_file[, grep('X', colnames(atac_file))]

    ## filter non-expressed genes
    minVal <- apply(rna_file, 1, min)
    rnaVal <- rna_file[which(minVal > 0), ]
    minVal <- apply(atac_file, 1, min)
    atacVal <- atac_file[which(minVal > 0), ]

    print(dim(rnaVal))
    print(dim(atacVal))

    ## match two data sets
    rnaVal <- rnaVal[rownames(rnaVal) %in% rownames(atacVal),]
    atacVal <- atacVal[rownames(atacVal) %in% rownames(rnaVal),]

    colnames(rnaVal) <- paste0(colnames(rnaVal),'_rna')
    colnames(atacVal) <- paste0(colnames(atacVal),'_atac')

    atacVal <- atacVal[match(rownames(atacVal), rownames(rnaVal)),]

    print(dim(rnaVal))
    print(dim(atacVal))

    ## ztransform
    get_zscore <- function(x){(x - mean(x)) /sd(x)}
    rnaZval <- t(apply(rnaVal, 1, get_zscore))
    rnaZval[is.na(rnaZval)] <- 0
    atacZval <-t(apply(atacVal, 1, get_zscore))
    atacZval[is.na(atacZval)] <- 0

    ## compute correlation
    df <- NULL
    for (rna_cell in colnames(rnaZval)){
	cor_vec <- c() 
	for (atac_cell in colnames(atacZval)){
	    corVal <- cor(rnaZval[,rna_cell], atacZval[,atac_cell])
	    cor_vec <- c(cor_vec, corVal)
	}   
	names(cor_vec) <- colnames(atacZval)
	df <- rbind(df, cor_vec)
    }
    rownames(df) <- colnames(rnaZval)
    df

    # save the cor matrix
    corMat <- data.frame(sampleID=rownames(df ), df)
    rownames(corMat) <- NULL
    write.table(corMat, file=paste0('CorMat/', cellType,'.corMat.txt'), row.names=F, col.names=T, sep="\t", quote=F)

    ## plot
    #valColor <- brewer.pal(9, 'YlGnBu')
    valColor <- colorRampPalette(c('yellow','white', 'navy'))(12)
    pdf(paste0('Plots/', cellType,'.corHeatmap.byIndiv.snRNA_snATAC.pdf'))
    pheatmap(as.matrix(df), color= valColor, 
      cluster_rows=T, cluster_cols=T,
      show_colnames=T, show_rownames=T,
      breaks=seq(-0.15, 0.15, length.out=12)
      )
    dev.off()

}






