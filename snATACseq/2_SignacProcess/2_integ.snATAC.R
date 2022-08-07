
library(Seurat)
library(Signac)
library(Rmagic)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
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
seuset.atac <- readRDS('SeuratObjects/PD.SN.snATAC.postAlign.rds')
maxDim <- 40

## Primary processing
seuset.atac <- FindTopFeatures(seuset.atac, min.cutoff = 50)
seuset.atac <- RunTFIDF(seuset.atac)
seuset.atac <- RunSVD(seuset.atac)

## SeqDepth correlation
pdf('Plots/DepthCor.pdf')
DepthCor(seuset.atac, n=maxDim)
dev.off()

##
seuset.atac <- RunUMAP(seuset.atac, reduction = "lsi", dims = 2:maxDim)

png('Plots/raw.png', height=6, width=6, res=1000, units='in')
p <- DimPlot(seuset.atac, group.by = "orig.ident", reduction='umap')
plot(p)
dev.off()

# Integration
seuset.list <- SplitObject(seuset.atac, split.by='orig.ident') 
integration.anchors <- FindIntegrationAnchors(
  object.list = seuset.list,
  anchor.features = rownames(seuset.atac),
  reduction = "rlsi",
  dims = 2:maxDim
)

seuset.integ <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = seuset.atac[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:maxDim
)

#
seuset.integ <- RunUMAP(object = seuset.integ, reduction = "integrated_lsi", dims = 2:maxDim)
seuset.integ <- FindNeighbors(object = seuset.integ, reduction = "integrated_lsi", dims = 2:maxDim)
seuset.integ <- FindClusters(object=seuset.integ, algorithm=3, resolution=1.5, verbose=FALSE)

png('Plots/integ.png', height=6, width=6, res=1000, units='in')
p <- DimPlot(seuset.integ, group.by = "orig.ident")
plot(p)
dev.off()

## Save the work
saveRDS(seuset.integ, 'SeuratObjects/PD.SN.snATAC.postAlign.integ.rds')

##===================================================================================================================================================
## Gene activity 
gene.activities <- GeneActivity(
  seuset.integ, assay = 'peaks', features = NULL,
  extend.upstream = 5000, extend.downstream = 0,
  biotypes = "protein_coding", max.width = NULL
)

seuset.integ[['RNA']] <- CreateAssayObject(counts = gene.activities)
seuset.integ <- NormalizeData(object=seuset.integ, assay='RNA', normalization.method = 'LogNormalize', scale.factor = median(seuset.integ$nCount_RNA))

## imputation
valMat <- as.matrix(t(seuset.integ@assays$RNA@data))
magic_data <- magic(valMat, )
magic_assay <- CreateAssayObject(counts = t(magic_data$result))
seuset.integ[['MAGIC_RNA']] <- magic_assay

## Save the work
saveRDS(seuset.integ, 'SeuratObjects/PD.SN.snATAC.postAlign.geneExp.rds')
#seuset.integ <- readRDS('SeuratObjects/PD.SN.snATAC.postAlign.geneExp.rds')

##===================================================================================================================================================
## Label transfer 
seuset.atac <- NULL
seuset.atac <- seuset.integ
seuset.integ <- NULL
seuset.rna <- readRDS('/home/ajl1213/Projects/PD/SciAdv_Data/snRNAseq/02_SeuratProcess/SeuratObjects/PD.SN.snRNA.postAlign.Anno.rds')

DefaultAssay(seuset.atac) <- 'RNA'

##
#seuset.rna <- FindVariableFeatures(seuset.rna, nfeatures=5000)
markerDF <- read.table('/home/ajl1213/Projects/PD/SciAdv_Data/snRNAseq/02_SeuratProcess/CellTypeDEGs/snRNA.cellTypeMarkers.txt', header=T)
markerList <- names(table(markerDF$gene))[which(table(markerDF$gene)==1)]
df <- markerDF[markerDF$gene %in% markerList,]
print(table(df$cluster))

transfer.anchors <- FindTransferAnchors(
  reference = seuset.rna, query = seuset.atac,
  features = markerList,
  #features = VariableFeatures(seuset.rna),
  reference.assay='RNA', query.assay='RNA',
  reduction = 'cca', dims = 1:maxDim,
  verbose = T
)

celltype.predictions <- TransferData(
  anchorset=transfer.anchors,
  refdata=seuset.rna$ident,
  weight.reduction=seuset.atac[["integrated_lsi"]],
  dims = 1:maxDim
)

prediction.DF <- data.frame(row.names=rownames(celltype.predictions), predicted.id = celltype.predictions$predicted.id, prediction.score.max = celltype.predictions$prediction.score.max)
seuset.atac <- AddMetaData(seuset.atac, prediction.DF)

seuset.atac$predicted.id <- factor(seuset.atac$predicted.id, levels=levels(seuset.rna$ident))
seuset.atac$orig.ident <- factor(seuset.atac$orig.ident, levels=sampleList)

##===================================================================================================================================================
## Save the work
saveRDS(seuset.atac, 'SeuratObjects/PD.SN.snATAC.postAlign.geneExp.label.rds')



