
library(Vennerable)

inputDIR <- '/home/ajl1213/Projects/PD/SciAdv_Anal/CombinedTargetID/01_FindTarget/RESULT'

dir.create('Plots')
minAbcScore <- 10


#======================================== Chow-Ruskey Plot
DopaN <- read.table(paste0(inputDIR,'/DopaN.selectList.txt'), header=T)
DopaN <- DopaN[which(DopaN$abcScore > minAbcScore),]

Oligo <- read.table(paste0(inputDIR, '/Oligo.selectList.txt'), header=T)
Oligo <- Oligo[which(Oligo$abcScore > minAbcScore),]

Ast <- read.table(paste0(inputDIR, '/Ast.selectList.txt'), header=T)
Ast <- Ast[which(Ast$abcScore > minAbcScore),]

Micro <- read.table(paste0(inputDIR, '/Micro.selectList.txt'), header=T)
Micro <- Micro[which(Micro$abcScore > minAbcScore),]

Endo <- read.table(paste0(inputDIR, '/Endo.selectList.txt'), header=T)
Endo <- Endo[which(Endo$abcScore > minAbcScore),]

vennList <- list(
Endo = Endo$geneID,
DopaN = DopaN$geneID, 
Oligo = Oligo$geneID, 
Ast = Ast$geneID, 
Micro = Micro$geneID
)
vennData <- Venn(vennList)

pdf('Plots/CellTypeTargetGene.ChowRusky.pdf')
plot(vennData, doWeights = TRUE, type = "ChowRuskey")
dev.off()




