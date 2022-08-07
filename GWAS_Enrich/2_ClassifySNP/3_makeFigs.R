

dir.create('Plots')



##
data1=read.table('SnpTable/MergedPD.PeakAsso.txt', header=T)
nTotal=nrow(data1)
nTotal

nonPeakSNP=data1[which(data1$peakLabel=='none'),]
peakSNP=data1[which(data1$peakLabel!='none'),]


##
nNonPeakSNP=nrow(nonPeakSNP)
nPeakSNP=nrow(peakSNP)

fracNonPeakSNP=round(nNonPeakSNP/nTotal*100, digits=2)
fracPeakSNP=round(nPeakSNP/nTotal*100, digits=2)

vals=c(fracNonPeakSNP, fracPeakSNP)
label=c('nonPeakSNP','PeakSNP')
label=paste(label, '(', vals,'%)', sep='')

pdf('Plots/PeakSNP.pie.pdf')
pie(vals, label, border='white')
dev.off()


##
nNonPeakSNP=nrow(nonPeakSNP)
nPlainPeakSNP=nrow(peakSNP[which(peakSNP$peakType=='nonDysPeak'),])
nSigPeakSNP=nrow(peakSNP[which(peakSNP$peakType=='dysPeak'),])

fracNonPeakSNP=round(nNonPeakSNP/nTotal*100, digits=2)
fracPlainPeakSNP=round(nPlainPeakSNP/nTotal*100, digits=2)
fracSigPeakSNP=round(nSigPeakSNP/nTotal*100, digits=2)

vals=c(fracNonPeakSNP, fracPlainPeakSNP, fracSigPeakSNP)
label=c('nonPeakSNP','PlainPeakSNP','SigPeakSNP')
label=paste(label, '(', vals,'%)', sep='')

print(vals)

pdf('Plots/PeakSNP.detail.pie.pdf')
pie(vals, label, border='white')
dev.off()
##






