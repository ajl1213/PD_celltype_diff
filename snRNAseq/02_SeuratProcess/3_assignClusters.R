#!/home/ajl1213/anaconda2/bin/R


library(Seurat)
library(Rmagic)
library(ggplot2)
library(RColorBrewer)
library(viridis)


seuset <- readRDS('SeuratObjects/PD.SN.snRNA.postAlign.rds')

##=======================================================================================================
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
# color:
nCtrl1 <- 6
nCtrl2 <- 7
nDisease <- 6
sample_colors <- c(rev(brewer.pal(9, 'YlGn'))[1:nCtrl1], rev(brewer.pal(9, 'YlGnBu'))[1:nCtrl2], rev(brewer.pal(9, 'RdPu'))[1:nDisease])
celltype_colors <- c(brewer.pal(8, 'Dark2'), brewer.pal(8, 'Accent')[5])
##=======================================================================================================
##
# NigralN - 12, 25
# Oligo - 0, 1, 2, 3, 5, 6, 7, 8, 14, 15, 20, 23, 26
# OPC - 4
# Ast - 10, 11, 16
# Micro - 9, 18, 19
# Endo - 13
# Peri - 17
# undetermined - 21, 22, 24, 27, 28

##=======================================================================================================
## cell-type assignment
seusetAnno <- SetIdent(seuset, WhichCells(seuset, idents=c(12, 25)), value='NigralN')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset, idents=c(0, 1, 2, 3, 5, 6, 7, 8, 14, 15, 20, 23, 26)), value='Oligo')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset, idents=c(4)), value='OPC')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset, idents=c(10, 11, 16)), value='Ast')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset, idents=c(9, 18, 19)), value='Micro')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset, idents=c(13)), value='Endo')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset, idents=c(17)), value='Peri')
seusetAnno <- SetIdent(seusetAnno, WhichCells(seuset, idents=c(21, 22, 24, 27, 28)), value='undetermined')

seusetAnno$ident <- seusetAnno@active.ident

##=======================================================================================================
## subclustering (neuronal cluster)
library(harmony)

neuron_c <- subset(seusetAnno, idents='NigralN')

#neuron_c <- NormalizeData(neuron_c)
neuron_c <- FindVariableFeatures(neuron_c, selection='vst', nfeatures=3000)
neuron_c <- ScaleData(neuron_c, verbose = FALSE)
neuron_c <- RunPCA(neuron_c, npcs = 45, verbose = FALSE)
neuron_c <- RunHarmony(object = neuron_c, group.by.vars = 'orig.ident', project.dim = F)
neuron_c <- RunUMAP(neuron_c, reduction = "harmony", dims = 1:5)
neuron_c <- FindNeighbors(neuron_c, reduction = "harmony", dims = 1:5)
neuron_c <- FindClusters(neuron_c, resolution= 1.0, graph.name='RNA_snn')

#
dir.create('subcluster')

png('subcluster/cluster.png', width=5, height=5, res=1000, units='in')
DimPlot(neuron_c, reduction='umap', label=T) + my_theme + ggtitle('')
dev.off()

png('subcluster/sample.png', width=5, height=5, res=1000, units='in')
DimPlot(neuron_c, reduction='umap', group.by='orig.ident', label=F, cols=sample_colors) + my_theme + ggtitle('')
dev.off()

# imputation
valMat <- t(as.matrix(neuron_c@assays$RNA@data))
magic_data <- magic(valMat)
magic_assay <- CreateAssayObject(counts = t(magic_data$result))
neuron_c[['MAGIC_RNA']] <- magic_assay
DefaultAssay(neuron_c) <- 'MAGIC_RNA'
for(cur_gene in c('SYT1','GAD1','GAD2')){
  print(cur_gene)
  png(paste0('subcluster/', cur_gene, '.snRNA.png'), width=5, height=5, res=1000, units='in')
  p <- FeaturePlot(neuron_c, features=cur_gene, cols=viridis(256), order=T, slot='data', min.cutoff='q0', max.cutoff='q95') + my_theme + ggtitle('')
  print(p)
  dev.off()
}
for(cur_gene in c('TH', 'SLC6A3')){
  print(cur_gene)
  png(paste0('subcluster/', cur_gene, '.snRNA.png'), width=5, height=5, res=1000, units='in')
  p <- FeaturePlot(neuron_c, features=cur_gene, cols=viridis(256), order=T, slot='data', min.cutoff='q0', max.cutoff='q80') + my_theme + ggtitle('')
  print(p)
  dev.off()
}

## cell-type assignment
# DopaN - 0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11
# GabaN - 5
# undetermined - NA

##
neuron_c.anno <- SetIdent(neuron_c, WhichCells(neuron_c, idents=c(0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11)), value='DopaN')
neuron_c.anno <- SetIdent(neuron_c.anno, WhichCells(neuron_c, idents=c(5)), value='GabaN')
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


##=======================================================================================================
## add nueron info
idx <- match(rownames(neuron_c.anno@meta.data), rownames(seusetAnno@meta.data))
seusetAnno$ident <- as.character(seusetAnno$ident)
seusetAnno$ident[idx] <- as.character(neuron_c.anno$ident)

## Filter
cellList <- c('DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri')
seusetAnno.filter <- subset(seusetAnno, subset= ident %in% cellList)
seusetAnno.filter$ident <- factor(seusetAnno.filter$ident, levels=cellList)
seusetAnno.filter@active.ident <- seusetAnno.filter$ident
seusetAnno.filter$disease <- factor(seusetAnno.filter$disease, levels=c('NOSN','PDSN'))

##=======================================================================================================
## Umap
png('Plots/snRNA_umap_celltype.png', width=6, height=6, res=1000, units='in')
DimPlot(seusetAnno.filter, reduction='umap', group.by='ident', cols=celltype_colors, label=F) + my_theme + ggtitle('')
dev.off()

png('Plots/snRNA_umap_celltype.label.png', width=6, height=6, res=1000, units='in')
DimPlot(seusetAnno.filter, reduction='umap', group.by='ident', cols=celltype_colors, label=T) + my_theme + ggtitle('')
dev.off()

## Cell number report
nCellTable <- table(seusetAnno.filter$ident, seusetAnno.filter$orig.ident)
sumVec <- apply(nCellTable, 2, sum)
fracTable <- NULL
for (i in 1:nrow(nCellTable)){
    tmpVec=nCellTable[i,]
    fracVec=tmpVec/sumVec*100
    fracTable=rbind(fracTable, fracVec)
}
sumVec
rownames(fracTable) <- rownames(nCellTable)
pdf('Plots/CellFrac.rna.bar.pdf', height=5, width=10, pointsize=2)
barplot(fracTable, col=celltype_colors, legend=rownames(fracTable), ylab='Proportion (%)', las=2)
dev.off()

table(seusetAnno.filter$disease)
table(seusetAnno.filter$ident)
table(seusetAnno.filter$ident)/sum(table(seusetAnno.filter$ident))*100

##=======================================================================================================
## Save the work
saveRDS(seusetAnno.filter, 'SeuratObjects/PD.SN.snRNA.postAlign.Anno.rds')


