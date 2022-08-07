library(Seurat)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(viridis)


load('CellTypeDEGs/snATAC.cellTypeMarkers.rda')
load('CellTypeDEGs/snATAC.cellTypePeaks.rda')
load('CellTypeDEGs/snATAC.cellTypeDEGs.rda')

dir.create('CellTypeDEGs/Figs')
cellTypeList=c('DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri')
##===============================================================================================================================================================================
## Dot plot
## V1 (Unique markers)
markerList <- names(table(snATAC_markers$gene))[which(table(snATAC_markers$gene)==1)]
df <- snATAC_markers[snATAC_markers$gene %in% markerList,]
print(table(df$cluster))

deg_list <- c()
for (cellType in cellTypeList){
    tmpDF <- df[which(df$cluster==cellType),]
    tmpList=head(tmpDF, 5)$gene
    deg_list=c(deg_list, tmpList)
}

df <- snATAC_markers[snATAC_markers$gene %in% deg_list,]
df$gene <- factor(df$gene, levels = unique(df$gene))
df$cluster <- factor(df$cluster, levels=rev(cellTypeList))

df$p_val_adj[which(df$p_val_adj == 0)] <- .Machine$double.xmin
df$logAdjPval <- -log10(df$p_val_adj)

p <- ggplot(df) +
  geom_point(aes(x=gene, y=cluster, size=logAdjPval, col=avg_log2FC)) + 
  scale_color_gradientn(colors= viridis_pal()(15)) + 
  labs(x='', y='') + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

pdf('CellTypeDEGs/Figs/snATAC.CellMarker_dotplot.Unique.pdf', width=14, height=3, useDingbats=F)
plot(p)
dev.off()

## V2 (Top markers)
data1 <- read.table('CellTypeDEGs/snATAC.cellTypeMarkers.inputPlot.txt', header=T)
data1$geneID <- factor(data1$geneID, levels = unique(data1$geneID))
data1$celltype <- factor(data1$celltype, levels=rev(cellTypeList))

p <- ggplot(data1) +
  geom_point(aes(x=geneID, y=celltype, size=logAdjPval, col=avgLogFC)) + 
  scale_color_gradientn(colors= viridis_pal()(15)) + 
  labs(x='', y='') + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

pdf('CellTypeDEGs/Figs/snATAC.CellMarker_dotplot.Top.pdf', width=14, height=3, useDingbats=F)
plot(p)
dev.off()

##===============================================================================================================================================================================
## Volcano plot

for (cellType in cellTypeList){
    print(cellType)
    tmpDF=snATAC_DEGs[which(snATAC_DEGs$celltype==cellType),]
    tmpDF$logP=-log10(tmpDF$p_val_adj)
    tmpDF$logP[tmpDF$logP > quantile(tmpDF$logP, 0.9999)]=quantile(tmpDF$logP, 0.9999)

    label=c()
    label[which(tmpDF$avg_log2FC<0)]='neg'
    label[which(tmpDF$avg_log2FC>=0)]='pos'
    tmpDF$label=label

    degList = tmpDF[which(tmpDF$p_val_adj < 0.05 & abs(tmpDF$avg_log2FC) > 0.1),]
    print(table(degList$label))

    negDF=tmpDF[which(tmpDF$avg_log2FC < 0),]
    negDF=negDF[order(negDF$logP, decreasing=T)[1:5],]
    posDF=tmpDF[which(tmpDF$avg_log2FC >= 0),]
    posDF=posDF[order(posDF$logP, decreasing=T)[1:5],]
    topDF=rbind(negDF, posDF)

    p= ggplot()+
        geom_hex(data=tmpDF, aes(x=avg_log2FC, y=logP), bins=200)+
        geom_text_repel(data=topDF, aes(x=avg_log2FC, y=logP, label=gene, color=label))+

        scale_fill_gradientn('Density', colours=c('black', 'navy','red','darkorange','yellow'), limits=c(0,50))+
        geom_vline(xintercept=0, linetype="dashed", color = "black")+
        labs(title=cellType, x='Avg_log2FC', y='-log10_pval_adj')+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(panel.background = element_rect(fill='white', color='black', linetype='solid')) +
        theme(legend.position='none')

    pdf(paste0('CellTypeDEGs/Figs/',cellType,'.scATAC_DEGs.volcano.pdf'), height=5, width=5, pointsize=3)
    plot(p)
    dev.off()
}




