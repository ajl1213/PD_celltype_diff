
library(Seurat)
library(Signac)
library(ggplot2)


seuset.atac=readRDS('SeuratObjects/PD.SN.snATAC.postAlign.geneExp.label.anno.rds')
dir.create('CellTypeDEGs')

##===============================================================================================================================================================================
Idents(seuset.atac) <- seuset.atac$ident

## Cell type marker DEGs
snATAC_markers <- FindAllMarkers(object = seuset.atac, test.use = "MAST", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, assay='RNA')
save(snATAC_markers,file="CellTypeDEGs/snATAC.cellTypeMarkers.rda")
write.table(snATAC_markers, file='CellTypeDEGs/snATAC.cellTypeMarkers.txt', col.names=T, row.names=F, sep='\t', quote=F)

##===============================================================================================================================================================================
## Control vs AD in each cell type
celltypes <- unique(seuset.atac$ident)

## Normalized count
combined_df <- data.frame()
for(i in 1:length(celltypes)){
    cur_celltype <- celltypes[[i]]
    print(cur_celltype)

    cur_seurat_obj <- subset(seuset.atac, ident == cur_celltype)
    Idents(cur_seurat_obj) <- cur_seurat_obj$disease
    print(unique(Idents(cur_seurat_obj)))

    #markers <- FindMarkers(object = cur_seurat_obj, test.use = "MAST", ident.1='PDSN', ident.2='NOSN', only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1, assay='RNA')
    markers <- FindMarkers(object = cur_seurat_obj, test.use = "MAST", ident.1='PDSN', ident.2='NOSN', only.pos = FALSE, min.pct = 0, logfc.threshold = 0, assay='RNA')
    markers$diff <- markers$pct.1 - markers$pct.2
    markers$celltype <- cur_celltype

    markers$gene <- rownames(markers); rownames(markers)=NULL
    combined_df <- rbind(combined_df, markers)
}

snATAC_DEGs <- combined_df
save(snATAC_DEGs, file="CellTypeDEGs/snATAC.cellTypeDEGs.rda")
write.table(snATAC_DEGs, file='CellTypeDEGs/snATAC.cellTypeDEGs.txt', col.names=T, row.names=F, sep='\t', quote=F)

##===============================================================================================================================================================================
## Find differentially accessible peaks in each cell type
DefaultAssay(seuset.atac) <- 'peaks'

combined_df <- data.frame()
for (i in 1:length(celltypes)){
    cur_celltype <- celltypes[[i]]
    print(cur_celltype)

    cur_seurat_obj <- subset(seuset.atac, ident == cur_celltype)
    Idents(cur_seurat_obj) <- cur_seurat_obj$disease
    print(unique(Idents(cur_seurat_obj)))

    da_peaks <- FindMarkers(object = cur_seurat_obj, ident.1 = "PDSN", ident.2 = "NOSN", min.pct = 0.05, test.use = 'LR', latent.vars = 'peak_region_fragments')
    da_peaks$diff <- da_peaks$pct.1 - da_peaks$pct.2
    da_peaks$celltype <- cur_celltype

    da_peaks$peakID <- rownames(da_peaks); rownames(da_peaks) <- NULL
    combined_df <- rbind(combined_df, da_peaks)
}

snATAC_DAPeaks <- combined_df
save(snATAC_DAPeaks, file='CellTypeDEGs/snATAC.DAPeaks.rda')
write.table(snATAC_DAPeaks, file='CellTypeDEGs/snATAC.DAPeaks.txt', col.names=T, row.names=F, sep='\t', quote=F)



