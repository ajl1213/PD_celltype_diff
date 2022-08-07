#!/home/ajl1213/anaconda2/bin/R


library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(viridis)


seuset <- readRDS('SeuratObjects/PD.SN.snRNA.postAlign.rds')

table(seuset@meta.data$orig.ident)
table(seuset@meta.data$disease)

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
##=======================================================================================================
## UMAP
dir.create('Plots')

## by cluster
png('Plots/snRNA_umap_cluster.png', width=10, height=10, res=1000, units='in')
DimPlot(seuset, reduction='umap', label=T) + my_theme + ggtitle('')
dev.off()
## by pathology
png('Plots/snRNA_umap_disease.png', width=10, height=10, res=1000, units='in')
DimPlot(seuset, reduction='umap', group.by='disease', label=F) + scale_color_manual(values=c(brewer.pal(9, 'Set1')[2], brewer.pal(9, 'Set1')[1])) + my_theme + ggtitle('')
dev.off()
## by sample
png('Plots/snRNA_umap_sample.png', width=10, height=10, res=1000, units='in')
DimPlot(seuset, reduction='umap', group.by='orig.ident', label=F, cols=sample_colors) + my_theme + ggtitle('')
dev.off()
## doublet score
png('Plots/snRNA_umap_doubletScore.png', width=10, height=10, res=1000, units='in')
FeaturePlot(seuset, reduction='umap', features='doubletScore', label=F, cols=rev(magma(30))) + my_theme + ggtitle('')
dev.off()

##=======================================================================================================
## CellMarker Distribution
DefaultAssay(seuset) <- 'RNA'
dir.create('CellMarkerDistrib')

## Neuron
png('CellMarkerDistrib/Neuron.SYT1.png',, width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='SYT1', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
## EN
png('CellMarkerDistrib/EN.SLC17A6.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='SLC17A6', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
png('CellMarkerDistrib/EN.SLC17A7.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='SLC17A7', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
## CADPS2Neuron
png('CellMarkerDistrib/CADPS2.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='CADPS2', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
## IN
png('CellMarkerDistrib/IN.GAD2.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='GAD2', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
png('CellMarkerDistrib/IN.GRIK1.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='GRIK1', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
## DopaNeuron
png('CellMarkerDistrib/DopaNeuron.TH.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='TH', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
png('CellMarkerDistrib/DopaNeuron.SLC6A3.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='SLC6A3', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
## Ast
png('CellMarkerDistrib/Ast.GFAP.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='GFAP', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
png('CellMarkerDistrib/Ast.AQP4.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='AQP4', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
## Micro
png('CellMarkerDistrib/Micro.ITGAM.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='ITGAM', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
png('CellMarkerDistrib/Micro.RUNX1.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='RUNX1', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
png('CellMarkerDistrib/Micro.CX3CR1.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='CX3CR1', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
png('CellMarkerDistrib/Micro.CD74.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='CD74', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
## Endo
png('CellMarkerDistrib/Endo.CLDN5.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='CLDN5', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
png('CellMarkerDistrib/Endo.RGS5.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='RGS5', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
png('CellMarkerDistrib/Endo.ICAM2.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='ICAM2', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
## Pericyte
png('CellMarkerDistrib/Peri.PDGFRB.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='PDGFRB', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
## Oligo
png('CellMarkerDistrib/Oligo.OLIG1.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='OLIG1', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
png('CellMarkerDistrib/Oligo.OPALIN.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='OPALIN', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
png('CellMarkerDistrib/Oligo.PLP1.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='PLP1', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
png('CellMarkerDistrib/Oligo.MOBP.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='MOBP', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
png('CellMarkerDistrib/Oligo.MAG.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='MAG', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
## OPC
png('CellMarkerDistrib/OPC.PDGFRA.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='PDGFRA', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()
png('CellMarkerDistrib/OPC.CSPG4.png', width=5, height=5, res=1000, units='in')
FeaturePlot(seuset, features='CSPG4', cols=viridis(256), order=T) + my_theme + ggtitle('')
dev.off()





