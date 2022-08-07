#!/env/maftools



library(maftools)

dir.create('Plots')

genomeVer <- 'hg19'
input_files <- c(
'annovar/X4870NOSN.hg19_multianno.exonic.txt',
'annovar/X4689NOSN.hg19_multianno.exonic.txt',
'annovar/X5628NOSN.hg19_multianno.exonic.txt',
'annovar/X5130_1NOSN.hg19_multianno.exonic.txt',
'annovar/X4996NOSN.hg19_multianno.exonic.txt',
'annovar/X5006_1NOSN.hg19_multianno.exonic.txt',
'annovar/X5244PDSN.hg19_multianno.exonic.txt',
'annovar/X4831PDSN.hg19_multianno.exonic.txt',
'annovar/X5215PDSN.hg19_multianno.exonic.txt',
'annovar/X5627PDSN.hg19_multianno.exonic.txt',
'annovar/X5778PDSN.hg19_multianno.exonic.txt',
'annovar/X5742PDSN.hg19_multianno.exonic.txt',
'annovar/X5649PDSN.hg19_multianno.exonic.txt',
'annovar/X5591PDSN.hg19_multianno.exonic.txt'
)

annovar.maf <- annovarToMaf(input_files, refBuild=genomeVer, table='refGene')
maf.data <- read.maf(maf=annovar.maf)

pdf('Plots/maf_summary.pdf')
plotmafSummary(maf=maf.data, rmOutlier=TRUE, addStat='median', dashboard=TRUE, titvRaw=FALSE, showBarcodes=TRUE)
dev.off()


#oncoplot(maf = maf.data, top = 20)





