
library(Seurat)
library(Signac)
library(ggplot2)
library(RColorBrewer)
library(viridis)
 

seuset.atac <- readRDS('SeuratObjects/PD.SN.snATAC.postAlign.geneExp.label.anno.rds')

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
## Save metadata
outDF <- seuset.atac@meta.data
outDF$barcodeID <- rownames(outDF); rownames(outDF) <- NULL
write.table(outDF, 'SeuratObjects/snATAC.metadata.FINAL.txt', row.names=F, col.names=T, quote=F, sep='\t')

##===================================================================================================================================================
## umap by sample
png('Plots/snATAC_umap_sample.final.png', width=8, height=8, res=1000, units='in')
DimPlot(seuset.atac, reduction='umap', group.by='orig.ident', label=F, cols=sample_colors) + my_theme + ggtitle('')
dev.off()

## sample barplot
seuset_meta <- seuset.atac@meta.data
variable <- 'orig.ident'
clusters <- unique(seuset_meta$ident)
df <- data.frame()
for(i in 1:length(clusters)){
  idx=which(seuset_meta$ident==clusters[i])
  cur_df <- table(seuset_meta$orig.ident[idx])
  cur_df <- as.data.frame(cur_df / table(seuset_meta[[variable]])[names(cur_df)])
  cur_df$Freq <- cur_df$Freq * 1/(sum(cur_df$Freq))
  cur_df$cluster <- clusters[i]
  df <- rbind(df, cur_df)
}

df$Var1 <- factor(df$Var1, levels=sampleList)

pdf('Plots/snATAC_barplot_sample.pdf', height=4, width=9)
p <- ggplot(df, aes(y=Freq, x=cluster, fill=Var1)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values=sample_colors) +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank()
  )
print(p)
dev.off()

##===================================================================================================================================================
## snATAC-seq umap
png('Plots/snATAC_umap_disease.final.png', width=8, height=8, res=1000, units='in')
DimPlot(seuset.atac, reduction='umap', group.by='disease', label=F) + scale_color_manual(values=c(brewer.pal(9, 'Set1')[2], brewer.pal(9, 'Set1')[1])) + my_theme + ggtitle('')
dev.off()

##===================================================================================================================================================
## doublet score
doublet_colors <- rev(magma(6))
png('Plots/snATAC_umap_doubletScore.final.png', width=8, height=8, res=1000, units='in')
FeaturePlot(seuset.atac, reduction='umap', features='doubletScore', label=F, cols=doublet_colors) + my_theme + ggtitle('')
dev.off()

## doublet score barplot
seuset_meta <- seuset.atac@meta.data
seuset_meta$doubletInterval <- cut(seuset_meta$doubletScore, breaks=6)
variable <- 'doubletInterval'
cellTypeList <- unique(seuset_meta$ident)
df <- data.frame()
for(i in 1:length(cellTypeList)){
  idx=which(seuset_meta$ident==cellTypeList[i])
  cur_df <- table(seuset_meta[[variable]][idx])
  #cur_df <- as.data.frame(cur_df / table(seuset_meta[[variable]])[names(cur_df)])
  cur_df <- as.data.frame(cur_df)
  cur_df$Freq <- cur_df$Freq * 1/(sum(cur_df$Freq))
  cur_df$cluster <- cellTypeList[i]
  df <- rbind(df, cur_df)
}
df$Var1 <- factor(df$Var1, levels=levels(seuset_meta[[variable]]))

pdf('Plots/snATAC_barplot_doublet.pdf', height=4, width=9)
p <- ggplot(df, aes(y=Freq, x=cluster, fill=Var1)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values=doublet_colors) +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank()
  )
plot(p)
dev.off()

##===================================================================================================================================================
## sex
sex_colors <- c('gold','dodgerblue')
seuset.atac$sex=factor(seuset.atac$sex, levels=c('Male','Female'))
png('Plots/snATAC_umap_sex.final.png', width=8, height=8, res=1000, units='in')
DimPlot(seuset.atac, reduction='umap', group.by='sex', cols=sex_colors, label=F) + my_theme + ggtitle('')
dev.off()

## sex barplot
seuset_meta <- seuset.atac@meta.data
variable <- 'sex'
cellTypeList <- unique(seuset_meta$ident)
df <- data.frame()
for(i in 1:length(cellTypeList)){
  idx=which(seuset_meta$ident==cellTypeList[i])
  cur_df <- table(seuset_meta[[variable]][idx])
  cur_df <- as.data.frame(cur_df / table(seuset_meta[[variable]])[names(cur_df)])
  cur_df$Freq <- cur_df$Freq * 1/(sum(cur_df$Freq))
  cur_df$cluster <- cellTypeList[i]
  df <- rbind(df, cur_df)
}
df$Var1 <- factor(df$Var1, levels=unique(seuset_meta[[variable]]))

pdf('Plots/snATAC_barplot_sex.pdf', height=4, width=9)
p <- ggplot(df, aes(y=Freq, x=cluster, fill=Var1)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values=sex_colors) +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank()
  )
plot(p)
dev.off()

##===================================================================================================================================================
## age
age_colors <- rev(magma(12)[3:7])
png('Plots/snATAC_umap_age.final.png', width=8, height=8, res=1000, units='in')
FeaturePlot(seuset.atac, reduction='umap', features='age', label=F, cols= age_colors) + my_theme + ggtitle('')
dev.off()

## age barplot
seuset_meta <- seuset.atac@meta.data
seuset_meta$ageInterval <- cut(seuset_meta$age, breaks=5)
seuset_meta$ageInterval=factor(seuset_meta$ageInterval, levels=rev(levels(seuset_meta$ageInterval)))

variable <- 'ageInterval'
cellTypeList <- unique(seuset_meta$ident)
df <- data.frame()
for(i in 1:length(cellTypeList)){
  idx=which(seuset_meta$ident==cellTypeList[i])
  cur_df <- table(seuset_meta[[variable]][idx])
  #cur_df <- as.data.frame(cur_df / table(seuset_meta[[variable]])[names(cur_df)])
  cur_df <- as.data.frame(cur_df)
  cur_df$Freq <- cur_df$Freq * 1/(sum(cur_df$Freq))
  cur_df$cluster <- cellTypeList[i]
  df <- rbind(df, cur_df)
}
df$Var1 <- factor(df$Var1, levels=levels(seuset_meta[[variable]]))

pdf('Plots/snATAC_barplot_age.pdf', height=4, width=9)
p <- ggplot(df, aes(y=Freq, x=cluster, fill=Var1)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values=age_colors) +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank()
  )
plot(p)
dev.off()


##=======================================================================================================
## post-mortem hour
pmi_colors <- rev(brewer.pal(9, 'RdPu'))
png('Plots/snATAC_umap_PMI.final.png', width=8, height=8, res=1000, units='in')
FeaturePlot(seuset.atac, reduction='umap', features='PMI', label=F, cols= rev(pmi_colors)) + my_theme + ggtitle('')
dev.off()

## post-mortem interval barplot
seuset_meta <- seuset.atac@meta.data
seuset_meta$pmiInterval <- cut(seuset_meta$PMI, breaks=9)
seuset_meta$pmiInterval=factor(seuset_meta$pmiInterval, levels=rev(levels(seuset_meta$pmiInterval)))

variable <- 'pmiInterval'
cellTypeList <- unique(seuset_meta$ident)
df <- data.frame()
for(i in 1:length(cellTypeList)){
  idx=which(seuset_meta$ident==cellTypeList[i])
  cur_df <- table(seuset_meta[[variable]][idx])
  #cur_df <- as.data.frame(cur_df / table(seuset_meta[[variable]])[names(cur_df)])
  cur_df <- as.data.frame(cur_df)
  cur_df$Freq <- cur_df$Freq * 1/(sum(cur_df$Freq))
  cur_df$cluster <- cellTypeList[i]
  df <- rbind(df, cur_df)
}
df$Var1 <- factor(df$Var1, levels=levels(seuset_meta[[variable]]))

pdf('Plots/snATAC_barplot_PMI.pdf', height=4, width=9)
p <- ggplot(df, aes(y=Freq, x=cluster, fill=Var1)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values=pmi_colors) +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank()
  )
plot(p)
dev.off()


##===================================================================================================================================================
## QC violin plots
# peak region fragments
pdf(paste0('Plots/snATAC_peakFrags.sample.violin.pdf'), width=8, height=4)
VlnPlot(seuset.atac, features='peak_region_fragments', pt.size=0, group.by='orig.ident', cols=sample_colors, ncol=1) +
  geom_boxplot(fill='white', outlier.shape=NA) + RotatedAxis() + NoLegend()
dev.off()

pdf(paste0('Plots/snATAC_peakFrags.celltype.violin.pdf'), width=8, height=4)
VlnPlot(seuset.atac, features='peak_region_fragments', pt.size=0, group.by='ident', cols=celltype_colors, ncol=1) +
  geom_boxplot(fill='white', outlier.shape=NA) + RotatedAxis() + NoLegend()
dev.off()

# nucleosome position
DefaultAssay(seuset.atac) <- 'peaks'
pdf('QualityControl/nucleosome.pdf')
FragmentHistogram(object = seuset.atac, group.by = 'orig.ident') + scale_fill_manual(values=sample_colors)
dev.off()


##===================================================================================================================================================
## Cell type proportions
library(ggpubr)
library(dplyr)

cellTypeList <- c('DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri')

seuset.atac$ident <- factor(seuset.atac$ident, levels=cellTypeList)
meta_list <- seuset.atac@meta.data %>% dplyr::group_split(orig.ident)

temp <- lapply(meta_list, function(meta){
  print(table(meta$disease))
  df <- as.data.frame(meta$ident %>% table / nrow(meta))
  colnames(df) <- c('cluster', 'proportion')
  df$sampleID <- unique(meta$orig.ident)
  df$disease <- unique(meta$disease)
  df
})
proportion_df <- Reduce(rbind, temp)
proportion_df$cluster_num <- as.numeric(proportion_df$cluster)

# relevel
proportion_df$cluster<- factor(
  as.character(proportion_df$cluster), levels=cellTypeList)

proportion_df$cluster_num <- as.numeric(proportion_df$cluster)
proportion_df$disease<- factor(proportion_df$disease, levels=rev(c('NOSN','PDSN')))

# box plot
p <- ggplot(proportion_df, aes(y=proportion, x=reorder(cluster, -cluster_num), fill=disease)) +
  geom_boxplot(outlier.shape=NA, color='black') +
  coord_flip() +
  stat_compare_means(method='wilcox.test', label='p.signif', label.y=0.9) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position="bottom",
  ) + xlab('') + theme_bw()

pdf('Plots/Celltype_composition_boxplot.pdf', width=5, height=4)
plot(p)
dev.off()

# pval
for (curr_cell in cellTypeList){
  print(curr_cell)
  tmpDF <- proportion_df[which(proportion_df$cluster == curr_cell),]
  head(tmpDF)
  ctrlVals <- tmpDF$proportion[which(tmpDF$disease=='NOSN')]
  expVals <- tmpDF$proportion[which(tmpDF$disease=='PDSN')]

  #
  print('t-test')
  test_res1 <- t.test(ctrlVals, expVals, alternative='less')
  test_res2 <- t.test(ctrlVals, expVals, alternative='greater')
  print(test_res1$p.value)
  print(test_res2$p.value)

  #
  print('wilcox-test')
  test_res1 <- wilcox.test(ctrlVals, expVals, alternative='less')
  test_res2 <- wilcox.test(ctrlVals, expVals, alternative='greater')
  print(test_res1$p.value)
  print(test_res2$p.value)
}


