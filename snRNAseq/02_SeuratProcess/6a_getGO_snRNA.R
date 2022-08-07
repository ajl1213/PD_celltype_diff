library(Seurat)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(viridis)



load('CellTypeDEGs/snRNA.cellTypeDEGs.rda')
dir.create('CellTypeDEGs/GO')

minPct <- 0.2
minLog2fc <- 0.2 
qval_thresh <- 0.05

##===============================================================================================================================================================================
## Get GO term
library(enrichR)

#a=listEnrichrDbs() 

filtered <- snRNA_DEGs[which(snRNA_DEGs$p_val_adj < qval_thresh & abs(snRNA_DEGs$avg_log2FC) > minLog2fc), ]
snRNA_DEGs.filtered <- filtered[which(filtered$pct.1 > minPct | filtered$pct.2 > minPct),]

##
#dbs <- c('GO_Biological_Process_2021', 'GO_Cellular_Component_2021', 'GO_Molecular_Function_2021')
dbs <- c('GO_Biological_Process_2021')
#dbs <- c('GO_Biological_Process_2018')
nGene <- 10

cellTypeList <- unique(snRNA_DEGs$celltype)

df.all=NULL
for (cellType in cellTypeList){
    print(cellType)
    tmpDF=snRNA_DEGs.filtered[which(snRNA_DEGs.filtered$celltype==cellType),]
    tmpDF$logP=-log10(tmpDF$p_val_adj)

    negDF=tmpDF[which(tmpDF$avg_log2FC < 0),]
    posDF=tmpDF[which(tmpDF$avg_log2FC >= 0),]

    if ((dim(negDF)[1] > nGene & dim(posDF)[1] > nGene)){

	#negDF=negDF[order(negDF$logP, decreasing=T)[1:nGene],]
	#posDF=posDF[order(posDF$logP, decreasing=T)[1:nGene],]
	negDF=negDF[order(negDF$logP, decreasing=T),]
	posDF=posDF[order(posDF$logP, decreasing=T),]

	negGO=enrichr(negDF$gene, dbs)
	posGO=enrichr(posDF$gene, dbs)

	for (db in dbs){
	    negGO[[db]]$cellType=cellType
	    negGO[[db]]$type='Down'
	    negGO[[db]]$db=db
	    posGO[[db]]$cellType=cellType
	    posGO[[db]]$type='Up'
	    posGO[[db]]$db=db

	    df.all=rbind(df.all, negGO[[db]])
	    df.all=rbind(df.all, posGO[[db]])
	}
    } else {
	print(paste0('insufficient # gene: ', cellType))
	}
}

#df.filtered <- df.all[which(df.all$Adjusted.P.value<0.05),]
df.filtered <- df.all[which(df.all$P.value<0.05),]

downDF <- df.filtered[which(df.filtered$type=='Down'),]
upDF <- df.filtered[which(df.filtered$type=='Up'),]

write.table(downDF, 'CellTypeDEGs/GO/DownDEGs.GO.txt', row.names=F, col.names=T, sep='\t', quote=F)
write.table(upDF, 'CellTypeDEGs/GO/UpDEGs.GO.txt', row.names=F, col.names=T, sep='\t', quote=F)



##===============================================================================================================================================================================
## Organize GO with python script (make_GO_data.py)
## make plot
cellList <- c('DopaN','Oligo','OPC','Ast','Micro','Endo','Peri')
colfunc <- colorRampPalette((brewer.pal(9, 'YlOrRd')[2:9]))

#
downPlotDF=read.table('CellTypeDEGs/GO/DownDEGs.plotData.txt', header=T)
downPlotDF$cellType=factor(downPlotDF$cellType, levels=cellList)
downPlotDF$termID=factor(downPlotDF$termID, levels=rev(unique(downPlotDF$termID)))
downPlotDF$logPval[downPlotDF$logPval > 5] <- 5
downPlotDF$oddRatio[downPlotDF$oddRatio > 300] <- 300

p <- ggplot() +
    geom_point(data=downPlotDF, aes(x=cellType, y=termID, size=oddRatio, col=logPval)) +

    scale_color_gradientn(colors=colfunc(30)) +

    labs(x='', y='') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf('CellTypeDEGs/GO/DownDEGs.GO.pdf', width=10, height=8, pointsize=3)
plot(p)
dev.off()

#
upPlotDF=read.table('CellTypeDEGs/GO/UpDEGs.plotData.txt', header=T)
upPlotDF$cellType=factor(upPlotDF$cellType, levels=cellList)
upPlotDF$termID=factor(upPlotDF$termID, levels=rev(unique(upPlotDF$termID)))
upPlotDF$logPval[upPlotDF$logPval > 5] <- 5
upPlotDF$oddRatio[upPlotDF$oddRatio > 300] <- 300

p <- ggplot() +
    geom_point(data=upPlotDF, aes(x=cellType, y=termID, size=oddRatio, col=logPval)) +

    scale_color_gradientn(colors=colfunc(30)) +

    labs(x='', y='') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf('CellTypeDEGs/GO/UpDEGs.GO.pdf', width=10, height=8, pointsize=3)
plot(p)
dev.off()




