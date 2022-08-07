
library(Seurat)
library(ggplot2)


seuset <- readRDS('SeuratObjects/PD.SN.snRNA.postAlign.Anno.rds')
dir.create('CellTypeDEGs')

##===============================================================================================================================================================================
## Cell type marker DEGs
Idents(seuset) <- seuset$ident
## Normalized count
snRNA_markers <- FindAllMarkers(object = seuset, test.use = "MAST", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, assay='RNA')
save(snRNA_markers,file="CellTypeDEGs/snRNA.cellTypeMarkers.rda")
write.table(snRNA_markers, file='CellTypeDEGs/snRNA.cellTypeMarkers.txt', col.names=T, row.names=F, sep='\t', quote=F)

##===============================================================================================================================================================================
## Control vs AD in each cell type
celltypes <- unique(seuset$ident)

## Normalized count
combined_df <- data.frame()
for(i in 1:length(celltypes)){
  cur_celltype <- celltypes[[i]]
  print(cur_celltype)

  cur_seurat_obj <- subset(seuset, ident == cur_celltype)
  Idents(cur_seurat_obj) <- cur_seurat_obj$disease
  print(unique(Idents(cur_seurat_obj)))

  #markers <- FindMarkers(object = cur_seurat_obj, test.use = "MAST", ident.1='PDSN', ident.2='NOSN', only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1, assay='RNA')
  markers <- FindMarkers(object = cur_seurat_obj, test.use = "MAST", ident.1='PDSN', ident.2='NOSN', only.pos = FALSE, min.pct = 0, logfc.threshold = 0, assay='RNA')
  markers$diff <- markers$pct.1 - markers$pct.2
  markers$celltype <- cur_celltype

  markers$gene <- rownames(markers); rownames(markers)=NULL
  combined_df <- rbind(combined_df, markers)
}

snRNA_DEGs <- combined_df
save(snRNA_DEGs, file="CellTypeDEGs/snRNA.cellTypeDEGs.rda")
write.table(snRNA_DEGs, file='CellTypeDEGs/snRNA.cellTypeDEGs.txt', col.names=T, row.names=F, sep='\t', quote=F)





