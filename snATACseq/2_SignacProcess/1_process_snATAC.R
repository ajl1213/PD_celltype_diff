library(Seurat)
library(Signac)
library(ggplot2)
library(MASS)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(patchwork)
library(rtracklayer)
library(viridis)


sampleList <- c(
'X4870NOSN',
'X4996NOSN',
'X5006NOSN',
'X5130NOSN',
'X5628NOSN',

'Public1NOSN',
'Public2NOSN',

'X4831PDSN',
'X5215PDSN',
'X5244PDSN',
'X5532PDSN',
'X5591PDSN',
'X5627PDSN',
'X5742PDSN',
'X5778PDSN'
)

##===================================================================================================================================================
## Load data
peaks <- Read10X_h5('/home/ajl1213/Projects/PD/SciAdv_Data/snATACseq/1_Process/MergedAll/outs/filtered_peak_bc_matrix.h5')
fragFile <- '/home/ajl1213/Projects/PD/SciAdv_Data/snATACseq/1_Process/MergedAll/outs/fragments.tsv.gz'

chrom_assay <- CreateChromatinAssay(counts=peaks, sep=c(":", "-"), genome='hg19', fragments=fragFile)

cell_meta <- read.csv('/home/ajl1213/Projects/PD/SciAdv_Data/snATACseq/1_Process/MergedAll/outs/singlecell.csv', header=T, row.names=1)
seuset.atac <-  CreateSeuratObject(counts=chrom_assay, assay="peaks", meta.data=cell_meta)

## gene annotation
#annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
#seqlevelsStyle(annotations) <- 'UCSC'
#genome(annotations) <- "hg19"
annotations <- import('~/genome.info/gencode/hg19.release38/GTF/gencode.v38lift37.annotation.gtf')
genome(annotations) <- 'hg19'
seqlevelsStyle(annotations) <- 'UCSC'
annotations$gene_biotype <- annotations$gene_type
Annotation(seuset.atac) <- annotations

## add sample metadata (Diagnosis, Age, Sex, etc)
metaTable <- read.table('/home/ajl1213/Projects/PD/SciAdv_Data/snATACseq/1_Process/DoubletStat/DoubletStatMerged.txt', header=T)
rownames(metaTable) <- metaTable$barcodeID; metaTable$barcodeID=NULL

metaTable <- metaTable[rownames(metaTable) %in% rownames(seuset.atac@meta.data),]
seuset.atac <- AddMetaData(seuset.atac, metadata = metaTable)

seuset.atac$orig.ident <- factor(seuset.atac$orig.ident, levels=sampleList)
table(seuset.atac@meta.data$orig.ident)

disease=c()
disease[grep('NOSN', seuset.atac@meta.data$orig.ident)]='NOSN'
disease[grep('PDSN', seuset.atac@meta.data$orig.ident)]='PDSN'
seuset.atac@meta.data$disease = disease
seuset.atac$disease=factor(disease, levels=c('NOSN','PDSN'))
table(seuset.atac@meta.data$disease)

##===================================================================================================================================================
## Quality control
seuset.atac$pct_reads_in_peaks <- seuset.atac$peak_region_fragments / seuset.atac$passed_filters * 100
seuset.atac$blacklist_ratio <- seuset.atac$blacklist_region_fragments / seuset.atac$peak_region_fragments
seuset.atac <- NucleosomeSignal(object=seuset.atac)
seuset.atac=TSSEnrichment(object=seuset.atac, fast=FALSE)

doublet_signal <- ifelse(seuset.atac$doubletStat == 'True', 'LikelyDoublet', 'Singlet')
table(doublet_signal)

summary(seuset.atac$peak_region_fragments)
lowfrag_signal <- ifelse(seuset.atac$peak_region_fragments < 2000, 'Lowfrag', 'pass')
table(lowfrag_signal)
highfrag_signal <- ifelse(seuset.atac$peak_region_fragments > 30000, 'highfrag', 'pass') 
table(highfrag_signal)

frip_signal <- ifelse(seuset.atac$pct_reads_in_peaks < 15, 'lowFRIP','pass')
table(frip_signal)

nuc_signal <- ifelse(seuset.atac$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
table(nuc_signal)
tss_signal <- ifelse(seuset.atac$TSS.enrichment > 2, 'High', 'Low')
table(tss_signal)

## QC plot
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

dir.create('QualityControl')
for (sample in levels(seuset.atac$orig.ident)){
  cur_data <- seuset.atac@meta.data[which(seuset.atac$orig.ident==sample),]
  densityVec <- get_density(log10(cur_data$nCount_peaks), cur_data$pct_reads_in_peaks, n=100)
  p1 <- ggplot(cur_data) + geom_point(aes(x=log10(cur_data$nCount_peaks), y=cur_data$pct_reads_in_peaks, color=densityVec)) +
  scale_color_viridis() +
  geom_vline(xintercept=log10(2000), color = "black") +
  labs(x='log10(nCount)', y='pct_reads_in_peaks') +
  theme_classic() + theme(legend.position='none')

  pdf(paste0('QualityControl/', sample,'.Frip.FragCount.pdf'), width=4, height=4, pointsize=3)
  plot(p1)
  dev.off()
}

seuset.atac <- subset(x=seuset.atac, subset = peak_region_fragments > 2000 & peak_region_fragments < 30000 & pct_reads_in_peaks > 15 & nucleosome_signal < 10 & TSS.enrichment > 2 & doubletStat =='False')
seuset.atac

##===================================================================================================================================================
## Save the work

outDF=seuset.atac@meta.data
outDF$barcodeID=rownames(outDF); rownames(outDF)=NULL
write.table(outDF, 'SeuratObjects/snATAC.metadata.RAW.txt', row.names=F, col.names=T, quote=F, sep='\t')

dir.create('SeuratObjects')
saveRDS(seuset.atac, 'SeuratObjects/PD.SN.snATAC.postAlign.rds')





