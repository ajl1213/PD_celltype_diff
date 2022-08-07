
library(plotrix)
library(RColorBrewer)



data1 <- read.table('GWAS_eQTL_HiC.comparison.count.txt', header=T)

#
eqtl_spec_target <- data1$eqtl_spec[1]
common_target <- data1$common[1]
total <- eqtl_spec_target+common_target

#
dataType <- c('common_target', 'eqtl_spec_target')
count <- c(common_target, eqtl_spec_target)
frac <- c(common_target/total, eqtl_spec_target/total)

#
df <- data.frame(dataType, count, frac)

fracVal <- round(df$frac, 4)
label <- paste(df$dataType,'\n',fracVal*100,'%', sep='')

#
pdf('Plots/GWAS_eQTL.target.3dPie.pdf', height=3, width=4, pointsize=3)
pie3D(rev(df$count), labels=rev(label), radius=1,  explode=0.08, col=c(brewer.pal(8, 'Set2')[8], brewer.pal(8, 'Set2')[1]), main='GWAS_eQTL_target')
dev.off()


hyperP <- sum(dhyper(data1$common[1]:data1$hic_spec[1], data1$eqtl_spec[1], data1$allAsso[1]-data1$eqtl_spec[1], data1$hic_spec[1]))

print(hyperP)


