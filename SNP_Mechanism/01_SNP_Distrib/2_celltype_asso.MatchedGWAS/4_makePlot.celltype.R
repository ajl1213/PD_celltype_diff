



data1 <- read.table('matchedGWAS.celltype_asso.txt', header=T)


orderedDF <- data1[order(data1$log2Enrich, decreasing=T), ]
print(orderedDF)

pdf('Plots/matchedGWAS.celltype_asso.bar.pdf')
barplot(orderedDF$log2Enrich, names.arg=orderedDF$CellType, las=2, ylab='Log2 enrichment of matched PD GWAS-SNPs')
dev.off()



