
library(dplyr)
library(Seurat)
library(patchwork)


args=commandArgs(TRUE)
sampleID=args[1]
print(sampleID)


### Load dataset
inputDIR=paste(sampleID,'/outs/filtered_feature_bc_matrix/', sep='')
data1 = Read10X(data.dir =inputDIR)

### Make Seurat object
seuset = CreateSeuratObject(counts = data1, project = sampleID, min.cells = 3, min.features = 200)
seuset

### mt fraction
seuset[["percent.mt"]] <- PercentageFeatureSet(seuset, pattern = "^MT-")

###
outFile=paste('Plots/',sampleID,'.QC.pdf', sep='')
pdf(outFile)
VlnPlot(seuset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
dev.off()

### Filter
seuset <- subset(seuset, subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)

### Normalization
seuset <- NormalizeData(seuset, normalization.method = "LogNormalize")

### Variable features
seuset <- FindVariableFeatures(seuset, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(1.2, 10))

###
all.genes <- rownames(seuset)
seuset <- ScaleData(seuset, features = all.genes)

### Linear dimensional reduction
seuset <- RunPCA(seuset, features = VariableFeatures(object = seuset))
seuset <- RunTSNE(seuset, reduction = "pca", dims = 1:40)
### Cluster the cells
seuset <- FindNeighbors(seuset, dims = 1:40)
seuset <- FindClusters(seuset, resolution = 0.5)

### Non-linear dimensional reduction (UMAP/tSNE)
outFile=paste('Plots/',sampleID,'.TSNE.pdf', sep='')
pdf(outFile)
TSNEPlot(seuset)
dev.off()

seuset



