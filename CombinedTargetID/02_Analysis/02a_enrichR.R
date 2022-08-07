#!/conda2/scEnv


library(RColorBrewer)
library(enrichR)


inputDIR <- '/home/ajl1213/Projects/PD/SciAdv_Anal/CombinedTargetID/01_FindTarget/RESULT'

#
#a <- listEnrichrDbs()
#print(unique(a$libraryName))

#
cellTypeList <- c('DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri')
minAbcScore <- 10

dir.create('GOTABLE')
geneList <- c()

## BP
dbs <- c('GO_Biological_Process_2021')
#dbs <- c('GO_Biological_Process_2018')
#dbs <- c('GO_Biological_Process_2015')

for (cellType in cellTypeList){
    print(cellType)

    data1 <- read.table(paste0(inputDIR, '/', cellType, '.selectList.txt'), header=T)
    data1 <- data1[which(data1$abcScore > minAbcScore),]
    geneList <- c(geneList, data1$geneID)

    downDF <- data1[which(data1$creType == 'DownPeak'),]
    upDF <- data1[which(data1$creType == 'UpPeak'),]

    downGO <- enrichr(downDF$geneID, dbs)
    upGO <- enrichr(upDF$geneID, dbs)

    df <- NULL
    for (db in dbs){
	downGO[[db]]$cellType=cellType
	downGO[[db]]$type='Down'
	downGO[[db]]$db=db
	upGO[[db]]$cellType=cellType
	upGO[[db]]$type='Up'
	upGO[[db]]$db=db

	df=rbind(df, downGO[[db]])
	df=rbind(df, upGO[[db]])
	}

    df <- df[order(df$P.value),]
    write.table(df, file=paste0('GOTABLE/', cellType,'.BP.txt'), row.names=F, col.names=T, sep='\t', quote=F)

    # metascape input
    print(length(downDF$geneID))
    print(length(upDF$geneID))
    write.table(downDF$geneID, file=paste0('GOTABLE/', cellType,'.Down.geneList.txt'), row.names=F, col.names=F, sep='\t',quote=F)
    write.table(upDF$geneID, file=paste0('GOTABLE/', cellType,'.Up.geneList.txt'), row.names=F, col.names=F, sep='\t',quote=F)
}


## Phenotype
dbs <- c('MGI_Mammalian_Phenotype_Level_4_2021')

geneList <- unique(geneList)
print(length(geneList))

MP <- enrichr(geneList, dbs)
df <- NULL
for (db in dbs){
    df <- rbind(df, MP[[db]])
}

write.table(df, file=paste0('GOTABLE/Merged.MP.txt'), row.names=F, col.names=T, sep='\t', quote=F)

#
for (cellType in cellTypeList){
    print(cellType)

    data1 <- read.table(paste0(inputDIR, '/', cellType, '.selectList.txt'), header=T)
    data1 <- data1[which(data1$abcScore > minAbcScore),]

    MP <- enrichr(data1$geneID, dbs)

    df <- NULL
    for (db in dbs){
	df <- rbind(df, MP[[db]])
	}
    write.table(df, file=paste0('GOTABLE/', cellType,'.MP.txt'), row.names=F, col.names=T, sep='\t', quote=F)
}

## OMIM_Disease
dbs <- c('OMIM_Disease')

geneList <- unique(geneList)
print(length(geneList))

OD <- enrichr(geneList, dbs)
df <- NULL
for (db in dbs){
    df <- rbind(df, OD[[db]])
}

write.table(df, file=paste0('GOTABLE/Merged.OMIM_Disease.txt'), row.names=F, col.names=T, sep='\t', quote=F)




