

library(edgeR)


dir.create('ValTable')

data1 <- as.data.frame(read.table('/home/ajl1213/Projects/PD/data/ChIP/2_CountMatrix/AdjustedBulk/H3K27ac.cellAdjusted.qq.combat.txt', header=T, row.names=1))

# 
diagnosis=c()
diagnosis[grep('NOSN', colnames(data1))] <- 'NOSN'
diagnosis[grep('PDSN', colnames(data1))] <- 'PDSN'

# Create a DGEList object
y = DGEList(counts=data1, group=diagnosis)
y

# TMM normalization 
#y_n <- calcNormFactors(y)
#y_n$samples

# Design matrix
#design <- model.matrix(~0+Gender+diagnosis)
design <- model.matrix(~diagnosis)
rownames(design) <- colnames(y)
design

# Estimate dispersion
y_nd <- estimateDisp(y, design, robust=TRUE)
y_nd$common.dispersion

# GLM fitting
fit <- glmQLFit(y_nd, design, robust=TRUE)

# identify differential peaks
qlf <- glmQLFTest(fit)

adjPval <- p.adjust(qlf$table$PValue, method="BH")
print(sum(qlf$table$PValue < 0.05))
print(sum(adjPval < 0.05))

##
df_final=data.frame(peakID=row.names(qlf$table), qlf$table, adjPval)
write.table(df_final, file=paste0('ValTable/H3K27ac.valTable.txt'), row.names=F, col.names=T, sep="\t", quote=F)




