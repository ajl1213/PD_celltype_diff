library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)


args=commandArgs(TRUE)
sampleID=args[1]
print(sampleID)


### Load data
inputFile=paste(sampleID,'/outs/filtered_peak_bc_matrix.h5', sep='')
counts=Read10X_h5(filename=inputFile)
inputFile=paste(sampleID,'/outs/singlecell.csv', sep='')
metadata=read.csv(file=inputFile, header=TRUE, row.names=1)

fragFile=paste(sampleID,'/outs/fragments.tsv.gz', sep='')
chrom_assay=CreateChromatinAssay(counts = counts, sep = c(":", "-"), genome = 'hg19',
 fragments=fragFile, min.cells = 1)

seuset=CreateSeuratObject(counts=chrom_assay, project = sampleID, assay="peaks", meta.data=metadata)

### get gene annotations from EnsDb
annotations=GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

### change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations)='UCSC'
genome(annotations)="hg19"

### add the gene information to the object
Annotation(seuset)=annotations

### QC metrics
### nucleosome signal score per cell
seuset=NucleosomeSignal(object=seuset)

### TSS enrichment score per cell
seuset=TSSEnrichment(object=seuset, fast=FALSE)

### add blacklist ratio and fraction of reads in peaks
seuset$pct_reads_in_peaks=seuset$peak_region_fragments / seuset$passed_filters * 100
seuset$blacklist_ratio=seuset$blacklist_region_fragments / seuset$peak_region_fragments

seuset$high.tss <- ifelse(seuset$TSS.enrichment > 2, 'High', 'Low')
table(seuset$high.tss)
seuset$nucleosome_group <- ifelse(seuset$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
table(seuset$nucleosome_group)

outFile=paste('Plots/',sampleID,'.QC.pdf', sep='')
pdf(outFile, width=9, height=5, pointsize=2)
VlnPlot(object=seuset, features=c('pct_reads_in_peaks', 'peak_region_fragments', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'), pt.size = 0, ncol = 5)
dev.off()

seuset <- subset(x=seuset, subset = peak_region_fragments > 1000 & peak_region_fragments < 100000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 4 & TSS.enrichment > 2)
seuset


### Normalization and linear dimensional reduction
seuset <- RunTFIDF(seuset)
seuset <- FindTopFeatures(seuset, min.cutoff = 'q0')
seuset <- RunSVD(seuset)

### Non-linear dimension reduction and clustering
seuset <- RunUMAP(object = seuset, reduction = 'lsi', dims = 2:30)
seuset <- FindNeighbors(object = seuset, reduction = 'lsi', dims = 2:30)
seuset <- FindClusters(object = seuset, verbose = FALSE, algorithm = 3)

outFile=paste('Plots/',sampleID,'.UMAP.pdf', sep='')
pdf(outFile)
DimPlot(object = seuset, label = TRUE) + NoLegend()
dev.off()




