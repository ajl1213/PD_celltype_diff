
library(pheatmap)
library(RColorBrewer)


dir.create('Plots')


cellTypeList <- c('DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri')
rnaVal <- read.table('/home/ajl1213/Projects/PD/SciAdv_Data/snRNAseq/03_GetRnaCount/CountByCellType/RawCount.byCellType.qq.txt', header=T, row.names=1)
atacVal <- read.table('/home/ajl1213/Projects/PD/SciAdv_Data/snATACseq/4_AtacCount/CountByCellType/RawCount.snATAC.byCellType.qq.txt', header=T, row.names=1)

## filter non-expressed genes
minVal <- apply(rnaVal, 1, min)
rnaVal <- rnaVal[which(minVal > 0), ]
minVal <- apply(atacVal, 1, min)
atacVal <- atacVal[which(minVal > 0), ]

print(dim(rnaVal))
print(dim(atacVal))

## match two data sets
rnaVal <- rnaVal[rownames(rnaVal) %in% rownames(atacVal),]
atacVal <- atacVal[rownames(atacVal) %in% rownames(rnaVal),]

colnames(rnaVal) <- paste0(colnames(rnaVal),'_rna')
colnames(atacVal) <- paste0(colnames(atacVal),'_atac')

atacVal <- atacVal[match(rownames(atacVal), rownames(rnaVal)),]

dim(atacVal)

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

## plot
#valColor <- brewer.pal(9, 'YlGnBu')
valColor <- colorRampPalette(c('yellow','white', 'navy'))(12)
pdf('Plots/corHeatmap.snRNA_snATAC.pdf')
pheatmap(as.matrix(df), color= valColor, 
  cluster_rows=F, cluster_cols=F,
  show_colnames=T, show_rownames=T,
  breaks=seq(-0.5, 0.5, length.out=12)
  )
dev.off()




