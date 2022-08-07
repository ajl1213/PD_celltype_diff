
library(plotrix)
library(RColorBrewer)



df <- read.table('GWAS_eQTL.asso.txt', header=T)


GWAS_overlap <- df$overlapGWAS[1]
GWAS_spec <- df$totalGWAS[1] - df$overlapGWAS[1]

dataTypeList <- c('GWAS_common','GWAS_spec')
countVal <- c(GWAS_overlap, GWAS_spec)
fracVal <- round(countVal / df$totalGWAS[1], 4)

fracVal <- round(df$frac, 4)
label <- paste(dataTypeList,'\n',fracVal*100,'%', sep='')

#
pdf('Plots/GWAS_eQTL.asso.3dPie.pdf', height=3, width=4, pointsize=3)
pie3D(rev(countVal), labels=rev(label), radius=1,  explode=0.08, col=c(brewer.pal(8, 'Set2')[8], brewer.pal(8, 'Set2')[3]), main='GWAS_eQTL_asso')
dev.off()




