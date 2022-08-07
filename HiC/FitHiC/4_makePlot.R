
dir.create('Plots')


#=================================================================================================================================================
## Number of valid read-pairs
sampleList <- c(
'X4870NOSN','X4689NOSN','X5628NOSN','X5130NOSN','X5006NOSN','X4996NOSN',
'X5244PDSN','X4831PDSN','X5215PDSN','X5649PDSN','X5591PDSN'
)

nPairVec <- c() 
for (sampleID in sampleList){
    data1 <- read.table(paste0('/home/ajl1213/Projects/PD/SciAdv_Data/HiC/Mapping/TransCis/', sampleID, '.TransCis.txt'), header=T)
    nPair <- data1$cis[1]
    nPairVec <- c(nPairVec, nPair)
}

df=data.frame(sampleList=sampleList, nPair=nPairVec)
df.order=df[order(df$nPair, decreasing=T), ]

print(df.order)

pdf('Plots/nReadPair.bar.pdf')
barplot(rev(df.order$nPair), names.arg=rev(df.order$sampleList), horiz=T, las=2, xlab='nReadPair', col='darkgrey')
dev.off()


#=================================================================================================================================================
## SigInter Distance Histogram
# Normal SN
data1 <- read.table('SigInter/MergedNOSN.5kb.1e-2.FitHiC.all.txt', header=T)
pdf('Plots/NormalSN.SigInterDist.hist.pdf', width=6, height=3, pointsize=3)
hist(data1$dist, breaks=150, col='darkgrey', xlab='Genomic distance (bp)', ylab='Frequency', main='Normal SN')
abline(v=mean(data1$dist), lty=2, lwd=1.5)
abline(v=median(data1$dist), lty=2, lwd=1.5)
print(length(data1$dist))
print(summary(data1$dist))
# PD SN
data1 <- read.table('SigInter/MergedPDSN.5kb.1e-2.FitHiC.all.txt', header=T)
pdf('Plots/PDSN.SigInterDist.hist.pdf', width=6, height=3, pointsize=3)
hist(data1$dist, breaks=150, col='darkgrey', xlab='Genomic distance (bp)', ylab='Frequency', main='PD SN')
abline(v=mean(data1$dist), lty=2, lwd=1.5)
abline(v=median(data1$dist), lty=2, lwd=1.5)
print(length(data1$dist))
print(summary(data1$dist))


#=================================================================================================================================================
## SigInter Categorization PieChart
data1=read.table('SigInter/MergedTotal.5kb.1e-2.FitHiC.ProCre.txt', header=T)

df=data.frame(table(data1$interLabel))
df$FracVal=round(df$Freq/nrow(data1)*100, digits=2)
df

ccFrac=df$FracVal[which(df$Var1=='CC')]
coFrac=df$FracVal[which(df$Var1=='CO')]
ooFrac=df$FracVal[which(df$Var1=='OO')]
poFrac=df$FracVal[which(df$Var1=='PO')]
ppFrac=df$FracVal[which(df$Var1=='PP')]
pcFrac=df$FracVal[which(df$Var1=='PC')]

vals=c(ccFrac, coFrac, ooFrac, poFrac, ppFrac, pcFrac)
label=c('CC','CO','OO','PO','PP','PC')
label=paste(label, '(', vals,'%)', sep='')

pdf('Plots/SigInter.pie.pdf', height=4, width=4, pointsize=3)
pie(vals, label, border='white')
dev.off()





