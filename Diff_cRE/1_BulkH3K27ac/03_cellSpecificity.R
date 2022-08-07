


CellLabelFile <- read.table('/home/ajl1213/Projects/PD/data/snATAC/04_GetAtacCount/CellLabel/AtacPeak.CellLabel.txt', header=T)
BulkDiffPeakFile <- read.table('/home/ajl1213/Projects/PD/data/ChIP/2_CountMatrix/ValTable/H3K27ac.valTable.filter.txt', header=T)

downPeaks <- BulkDiffPeakFile$peakID[which(BulkDiffPeakFile$label =='Down')]
upPeaks <- BulkDiffPeakFile$peakID[which(BulkDiffPeakFile$label =='Up')]

downDF <- CellLabelFile[CellLabelFile$peakID %in% downPeaks,]
upDF <-  CellLabelFile[CellLabelFile$peakID %in% upPeaks,]


table(downDF$nCellType)
table(upDF$nCellType)



table(downDF$nCellType) / sum(table(downDF$nCellType))
table(upDF$nCellType) / sum(table(upDF$nCellType))

