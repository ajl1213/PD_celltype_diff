
library(Seurat)
library(Signac)
library(ggplot2)
library(RColorBrewer)
library(viridis)


dir.create('Plots')

##===================================================================================================================================================
## theme:
my_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  legend.position="none",
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)

# sample list:
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

# color:
nCtrl <- 7
nDisease <- 8
sample_colors <- c(rev(brewer.pal(9, 'YlGnBu'))[1:nCtrl], rev(brewer.pal(9, 'RdPu'))[1:nDisease])
celltype_colors <- c(brewer.pal(8, 'Dark2'), brewer.pal(8, 'Accent')[5])
##===================================================================================================================================================
##
seuset.atac <- readRDS('SeuratObjects/PD.SN.snATAC.postAlign.geneExp.label.rds')

# clusters
png('Plots/snATAC.cluster.png', width=8, height=8, res=1000, units='in')
p1 <- DimPlot(seuset.atac, group.by = 'seurat_clusters', label = T, repel = T) + my_theme + ggtitle('')
plot(p1)
dev.off()

# prediction without label
png('Plots/snATAC.cellprediction.png', width=8, height=8, res=1000, units='in')
p1 <- DimPlot(seuset.atac, group.by = 'predicted.id', label = F, repel = T, cols = celltype_colors) + my_theme + ggtitle('')
plot(p1)
dev.off()

# prediction with label
png('Plots/snATAC.cellprediction.label.png', width=8, height=8, res=1000, units='in')
p1 <- DimPlot(seuset.atac, group.by = 'predicted.id', label = T, repel = T, cols = celltype_colors) + my_theme + ggtitle('')
plot(p1)
dev.off()

# prediction score
predictionscore_colors <- rev(magma(6))
png('Plots/snATAC.predictionscore.png', width=8, height=8, res=1000, units='in')
FeaturePlot(seuset.atac, reduction='umap', features='prediction.score.max', label=F, cols=predictionscore_colors) + my_theme + ggtitle('')
dev.off()

# prediction score barplot
pdf('Plots/snATAC.predictionscore.bar.pdf')
hist(seuset.atac$prediction.score.max, breaks=100, xlim=c(0, 1))
abline(v=0.6, lty=2, col='grey')
dev.off()

print(sum(seuset.atac$prediction.score.max > 0.6))
print(sum(seuset.atac$prediction.score.max > 0.6)/length(seuset.atac$prediction.score.max))


# doubletScore
doubletScore_colors <- rev(magma(6))
png('Plots/snATAC.doubletscore.png', width=8, height=8, res=1000, units='in')
FeaturePlot(seuset.atac, reduction='umap', features='doubletScore', label=F, cols=doubletScore_colors) + my_theme + ggtitle('')
dev.off()

# pct_reads_in_peaks
pct_in_peaks_colors <- rev(magma(6))
png('Plots/snATAC.pct_in_peaks.png', width=8, height=8, res=1000, units='in')
FeaturePlot(seuset.atac, reduction='umap', features='pct_reads_in_peaks', label=F, cols=pct_in_peaks_colors) + my_theme + ggtitle('')
dev.off()

# feature plot
DefaultAssay(seuset.atac) <- 'MAGIC_RNA'
dir.create('MarkerDistrib')
for(cur_gene in c('AQP4', 'TH', 'SLC6A3','CLDN5','GAD1','GAD2','CD74','MAG','MOBP','PDGFRA','PDGFRB','SYT1')){
  print(cur_gene)
  png(paste0('MarkerDistrib/', cur_gene, '.snATAC.png'), width=5, height=5, res=1000, units='in')
  p <- FeaturePlot(seuset.atac, features=cur_gene, cols=viridis(256), order=T, slot='data') + my_theme + ggtitle('')
  print(p)
  dev.off()
}

##===================================================================================================================================================
##
# NigralN - 15
# Oligo - 0, 1, 2, 3, 4, 5, 10, 11, 14
# OPC - 9
# Ast - 8, 12, 17, 19
# Micro - 6, 7, 13
# Endo - 16
# Peri - 18
# undetermined - 20

## cell-type assignment
seusetAnno <- SetIdent(seuset.atac, WhichCells(seuset.atac, idents=c(15)), value='NigralN')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset.atac, idents=c(0, 1, 2, 3, 4, 5, 10, 11, 14)), value='Oligo')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset.atac, idents=c(9)), value='OPC')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset.atac, idents=c(8, 12, 17, 19)), value='Ast')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset.atac, idents=c(6, 7, 13)), value='Micro')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset.atac, idents=c(16)), value='Endo')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset.atac, idents=c(18)), value='Peri')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset.atac, idents=c(20)), value='undetermined')

seusetAnno$ident <- seusetAnno@active.ident


##===================================================================================================================================================
## subclustering (neuronal cluster)
library(harmony)

neuron_c <- subset(seusetAnno, idents='NigralN')
DefaultAssay(neuron_c) <- 'peaks'

#
neuron_c <- FindTopFeatures(neuron_c, min.cutoff = 10)
neuron_c <- RunTFIDF(neuron_c)
neuron_c <- RunSVD(neuron_c)

neuron_c <- RunHarmony(object = neuron_c, group.by.vars = 'orig.ident', reduction='lsi', assay.use='peaks', project.dim = F)
neuron_c <- RunUMAP(neuron_c, reduction = "harmony", dims = 1:5)
neuron_c <- FindNeighbors(neuron_c, reduction = "harmony", dims = 1:5)
neuron_c <- FindClusters(neuron_c, resolution= 1.0)

#
dir.create('subcluster')

png('subcluster/cluster.png', width=5, height=5, res=1000, units='in')
DimPlot(neuron_c, reduction='umap', label=T) + my_theme + ggtitle('')
dev.off()

png('subcluster/sample.png', width=5, height=5, res=1000, units='in')
DimPlot(neuron_c, reduction='umap', group.by='orig.ident', label=F, cols=sample_colors) + my_theme + ggtitle('')
dev.off()

DefaultAssay(neuron_c) <- 'MAGIC_RNA'
for(cur_gene in c('SYT1', 'TH', 'SLC6A3','GAD1','GAD2')){
  print(cur_gene)
  png(paste0('subcluster/', cur_gene, '.snATAC.png'), width=5, height=5, res=1000, units='in')
  p <- FeaturePlot(neuron_c, features=cur_gene, cols=viridis(256), order=T, slot='data', min.cutoff='q5', max.cutoff='q95') + my_theme + ggtitle('')
  print(p)
  dev.off()
}

## cell-type assignment
# DopaN - 0, 1, 2, 3, 4, 5, 6, 7, 8
# GabaN - 4 
# undetermined - None
##
neuron_c.anno <- SetIdent(neuron_c, WhichCells(neuron_c, idents=c(0, 1, 2, 3, 4, 5, 6, 7, 8)), value='DopaN')
neuron_c.anno <- SetIdent(neuron_c.anno, WhichCells(neuron_c, idents=c(4)), value='GabaN')
#neuron_c.anno <- SetIdent(neuron_c.anno, WhichCells(neuron_c, idents=c()), value='undetermined')
neuron_c.anno$ident <- neuron_c.anno@active.ident
neuron_c.anno.filter <- subset(neuron_c.anno, idents=c('DopaN','GabaN'))
neuron_c.anno.filter@active.ident <- factor(neuron_c.anno.filter@active.ident, levels=c('DopaN', 'GabaN'))
neuron_c.anno.filter$ident <- neuron_c.anno.filter@active.ident

##
png('subcluster/sub_neuron.png', width=5, height=5, res=1000, units='in')
DimPlot(neuron_c.anno.filter, reduction='umap', group.by='ident', cols=c(brewer.pal(8, 'Dark2')[1], brewer.pal(8, 'Dark2')[2]), label=F) + my_theme + ggtitle('')
dev.off()

png('subcluster/sub_neuron.label.png', width=5, height=5, res=1000, units='in')
DimPlot(neuron_c.anno.filter, reduction='umap', group.by='ident', cols=c(brewer.pal(8, 'Dark2')[1], brewer.pal(8, 'Dark2')[2]), label=T) + my_theme + ggtitle('')
dev.off()

DefaultAssay(neuron_c.anno.filter) <- 'MAGIC_RNA'
for(cur_gene in c('TH', 'SLC6A3')){
  print(cur_gene)
  png(paste0('subcluster/', cur_gene, '.snATAC.png'), width=5, height=5, res=1000, units='in')
  p <- FeaturePlot(neuron_c.anno.filter, features=cur_gene, cols=viridis(256), order=T, slot='data', min.cutoff='q20', max.cutoff='q80') + my_theme + ggtitle('')
  print(p)
  dev.off()
}

for(cur_gene in c('GAD1','GAD2')){
  print(cur_gene)
  png(paste0('subcluster/', cur_gene, '.snATAC.png'), width=5, height=5, res=1000, units='in')
  p <- FeaturePlot(neuron_c.anno.filter, features=cur_gene, cols=viridis(256), order=T, slot='data', min.cutoff='q20', max.cutoff='q95') + my_theme + ggtitle('')
  print(p)
  dev.off()
}


##===================================================================================================================================================
## add neuronal info
idx <- match(rownames(neuron_c.anno@meta.data), rownames(seusetAnno@meta.data))
seusetAnno$ident <- as.character(seusetAnno$ident)
seusetAnno$ident[idx] <- as.character(neuron_c.anno$ident)

## Filter
cellList <- c('DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri')
seusetAnno.filter <- subset(seusetAnno, idents=c('undetermined'), invert=TRUE)
seusetAnno.final <- subset(seusetAnno.filter, subset = prediction.score.max < 0.6, invert = TRUE)

seusetAnno.final$ident <- factor(seusetAnno.final$ident, levels=cellList)
seusetAnno.final@active.ident <- seusetAnno.final$ident

seusetAnno.final$disease <- factor(seusetAnno.final$disease, levels=c('NOSN','PDSN'))

# feature plot
DefaultAssay(seusetAnno.final) <- 'MAGIC_RNA'
dir.create('MarkerDistribFinal')
for(cur_gene in c('MAG','MOBP','PDGFRA')){
  print(cur_gene)
  png(paste0('MarkerDistribFinal/', cur_gene, '.snATAC.png'), width=5, height=5, res=1000, units='in')
  p <- FeaturePlot(seusetAnno.final, features=cur_gene, cols=viridis(256), order=T, slot='data', min.cutoff='q20', max.cutoff='q98') + my_theme + ggtitle('')
  print(p)
  dev.off()
}

DefaultAssay(seusetAnno.final) <- 'MAGIC_RNA'
for(cur_gene in c('AQP4','GFAP','CLDN5','CD74','RUNX1','PDGFRB','SYT1')){
  print(cur_gene)
  png(paste0('MarkerDistribFinal/', cur_gene, '.snATAC.png'), width=5, height=5, res=1000, units='in')
  p <- FeaturePlot(seusetAnno.final, features=cur_gene, cols=viridis(256), order=T, slot='data', min.cutoff='q60', max.cutoff='q100') + my_theme + ggtitle('')
  print(p)
  dev.off()
}

DefaultAssay(seusetAnno.final) <- 'MAGIC_RNA'
for(cur_gene in c('TH','SLC6A3')){
  print(cur_gene)
  png(paste0('MarkerDistribFinal/', cur_gene, '.snATAC.png'), width=5, height=5, res=1000, units='in')
  p <- FeaturePlot(seusetAnno.final, features=cur_gene, cols=viridis(256), order=T, slot='data', min.cutoff='q80', max.cutoff='q100') + my_theme + ggtitle('')
  print(p)
  dev.off()
}

DefaultAssay(seusetAnno.final) <- 'RNA'
for(cur_gene in c('GAD1','GAD2')){
  print(cur_gene)
  png(paste0('MarkerDistribFinal/', cur_gene, '.snATAC.png'), width=5, height=5, res=1000, units='in')
  p <- FeaturePlot(seusetAnno.final, features=cur_gene, cols=viridis(256), order=T, slot='data', min.cutoff='q70', max.cutoff='q98') + my_theme + ggtitle('')
  print(p)
  dev.off()
}

# Umap
png('Plots/snATAC_umap_celltype.png', width=6, height=6, res=1000, units='in')
DimPlot(seusetAnno.final, reduction='umap', group.by='ident', cols=celltype_colors, label=F) + my_theme + ggtitle('')
dev.off()

png('Plots/snATAC_umap_celltype.label.png', width=6, height=6, res=1000, units='in')
DimPlot(seusetAnno.final, reduction='umap', group.by='ident', cols=celltype_colors, label=T) + my_theme + ggtitle('')
dev.off()

##
table(seusetAnno.final$ident)
table(seusetAnno.final$ident)/sum(table(seusetAnno.final$ident))*100
table(seusetAnno.final$orig.ident)
table(seusetAnno.final$disease)


##===================================================================================================================================================
## Save the work
saveRDS(seusetAnno.final, 'SeuratObjects/PD.SN.snATAC.postAlign.geneExp.label.anno.rds')


