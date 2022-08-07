
library(Signac)
library(GenomicRanges)
library(ggplot2)


seuset.atac <- readRDS('/home/ajl1213/Projects/PD/SciAdv_Data/snATACseq/2_SignacProcess/SeuratObjects/PD.SN.snATAC.postAlign.geneExp.label.anno.rds')
Macs2Path <- '/home/ajl1213/bin/macs2'

#================================ Peakcall
seuset.atac@active.assay <- 'peaks'
peakObj <- CallPeaks(object = seuset.atac, group.by = "ident", macs2.path = Macs2Path)
df <- data.frame(chrID=as.character(peakObj@seqnames), start=start(peakObj), end=end(peakObj), cellID=peakObj$peak_called_in)

dir.create('PeakList')
write.table(df, file='PeakList/snATAC.peaks.bed', row.names=F, col.names=F, sep="\t", quote=F)


#================================ Extract barcodes
dir.create('Barcodes')

# SampleID
for (sampleID in levels(seuset.atac@meta.data$orig.ident)){
    df=data.frame(barcodes=rownames(seuset.atac@meta.data)[which(seuset.atac@meta.data$orig.ident==sampleID)], cellType=seuset.atac@meta.data$ident[seuset.atac@meta.data$orig.ident==sampleID])
    write.table(df, file=paste0('Barcodes/',sampleID,'.AllCells.barcodes.txt'), row.names=F, col.names=T, sep="\t", quote=F)
}

# Cell-Type
for (sampleID in levels(seuset.atac@meta.data$orig.ident)){
    df=data.frame(barcodes=rownames(seuset.atac@meta.data)[which(seuset.atac@meta.data$orig.ident==sampleID)], cellType=seuset.atac@meta.data$ident[seuset.atac@meta.data$orig.ident==sampleID])

    for (cellType in levels(seuset.atac@meta.data$ident)){
        tmp=data.frame(barcodes=df$barcodes[which(df$cellType==cellType)], cellType=df$cellType[which(df$cellType==cellType)])
        write.table(tmp, file=paste0('Barcodes/', sampleID,'.',cellType,'.barcodes.txt'), row.names=F, col.names=T, sep='\t', quote=F)
        }   
}




