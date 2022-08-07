#!/home/ajl1213/anaconda2/bin/R


library(Seurat)
library(dplyr)


dir.create('SeuratObjects')

##===================================================================================================================================================
## Load data
data1 <- Read10X(data.dir='/home/ajl1213/Projects/PD/SciAdv_Data/snRNAseq/01_Process/MergedAll/outs/filtered_feature_bc_matrix') 
seuset <- CreateSeuratObject(counts = data1, project = "PD.SN.snRNA")

## Load sample metadata
splitVec <- strsplit(rownames(seuset@meta.data),'-')
df <- data.frame(matrix(unlist(splitVec),ncol=2, byrow=TRUE))

X5628NOSN    <-which(df$X2==1)
X4689NOSN    <-which(df$X2==2)
X4870NOSN    <-which(df$X2==3)
X4996NOSN    <-which(df$X2==4)
X5006NOSN    <-which(df$X2==5)
X5130NOSN    <-which(df$X2==6)

Public1NOSN  <-which(df$X2==7)
Public2_1NOSN<-which(df$X2==8)
Public2_2NOSN<-which(df$X2==9)
Public3_1NOSN<-which(df$X2==10)
Public3_2NOSN<-which(df$X2==11)
Public4NOSN  <-which(df$X2==12)
Public5NOSN  <-which(df$X2==13)

X5591PDSN    <-which(df$X2==14)
X4831PDSN    <-which(df$X2==15)
X5215PDSN    <-which(df$X2==16)
X5244PDSN    <-which(df$X2==17)
X5742PDSN    <-which(df$X2==18)
X5778PDSN    <-which(df$X2==19)

orig <- c(
rep('X5628NOSN',length(X5628NOSN)),
rep('X4689NOSN',length(X4689NOSN)),
rep('X4870NOSN',length(X4870NOSN)),
rep('X4996NOSN',length(X4996NOSN)),
rep('X5006NOSN',length(X5006NOSN)),
rep('X5130NOSN',length(X5130NOSN)),

rep('Public1NOSN',length(Public1NOSN)),
rep('Public2_1NOSN',length(Public2_1NOSN)),
rep('Public2_2NOSN',length(Public2_2NOSN)),
rep('Public3_1NOSN',length(Public3_1NOSN)),
rep('Public3_2NOSN',length(Public3_2NOSN)),
rep('Public4NOSN',length(Public4NOSN)),
rep('Public5NOSN',length(Public5NOSN)),

rep('X5591PDSN',length(X5591PDSN)),
rep('X4831PDSN',length(X4831PDSN)),
rep('X5215PDSN',length(X5215PDSN)),
rep('X5244PDSN',length(X5244PDSN)),
rep('X5742PDSN',length(X5742PDSN)),
rep('X5778PDSN',length(X5778PDSN))
)

orig <- factor(orig, levels=c(
'X5628NOSN',
'X4689NOSN',
'X4870NOSN',
'X4996NOSN',
'X5006NOSN',
'X5130NOSN',

'Public1NOSN',
'Public2_1NOSN',
'Public2_2NOSN',
'Public3_1NOSN',
'Public3_2NOSN',
'Public4NOSN',
'Public5NOSN',

'X5591PDSN',
'X4831PDSN',
'X5215PDSN',
'X5244PDSN',
'X5742PDSN',
'X5778PDSN'
))

seuset@meta.data$orig.ident <- orig
table(orig)

disease <- c()
disease[grep('NOSN', seuset@meta.data$orig.ident)] <- 'NOSN'
disease[grep('PDSN', seuset@meta.data$orig.ident)] <- 'PDSN'

disease <-factor(disease, levels=c('NOSN','PDSN'))
seuset@meta.data$disease <- disease

##========================================================================================================================================================================================
## Load doublet stat
metaTable <- read.table('/home/ajl1213/Projects/PD/SciAdv_Data/snRNAseq/01_Process/DoubletStat/DoubletStatMerged.txt', header=T)
doubletData <- metaTable
rownames(doubletData) <- doubletData$barcodeID; doubletData$barcodeID <- NULL; doubletData$sampleID <- NULL
doubletData <- doubletData[rownames(doubletData) %in% rownames(seuset@meta.data),]
seuset <- AddMetaData(seuset, metadata=doubletData)
seuset <- subset(seuset, subset= doubletStat!='True')

##========================================================================================================================================================================================
## Quality control
seuset[['percent.mt']] <- PercentageFeatureSet(seuset, pattern = "^MT-")

## Save metadata
outDF <- seuset@meta.data
outDF$barcodeID <- rownames(outDF); rownames(outDF) <- NULL
write.table(outDF, 'SeuratObjects/snRNA.metadata.RAW.txt', row.names=F, col.names=T, quote=F, sep='\t')

seuset <- subset(seuset, subset= nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10 & doubletScore < 0.4)
seuset <- seuset[!grepl("^MT-", rownames(seuset)),]

##========================================================================================================================================================================================
## Split object by sample
seuset.list <- SplitObject(seuset, split.by='orig.ident')
maxDim <- 45

## Normalization and variable feature extraction
seuset.list <- lapply(X = seuset.list, FUN = function(x) {
    x <- NormalizeData(x)
    #x <- FindVariableFeatures(x, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(1.2, 10))
    x <- FindVariableFeatures(x, nfeatures=5000)
})

## Integrate
condition.anchor <- FindIntegrationAnchors(object.list = seuset.list, dims = 1:maxDim)
condition.combined <- IntegrateData(anchorset = condition.anchor, dims = 1:maxDim)
##========================================================================================================================================================================================

## Run the standard workflow for visualization and clustering
condition.combined <- ScaleData(condition.combined, verbose = FALSE)
condition.combined <- RunPCA(condition.combined, npcs = maxDim, verbose = FALSE)

## t-SNE and Clustering
condition.combined <- RunUMAP(condition.combined, reduction = "pca", dims = 1:maxDim)
condition.combined <- FindNeighbors(condition.combined, reduction = "pca", dims = 1:maxDim)
condition.combined <- FindClusters(condition.combined, resolution= 1.5)

## Save the work
saveRDS(condition.combined, 'SeuratObjects/PD.SN.snRNA.postAlign.rds')
##



