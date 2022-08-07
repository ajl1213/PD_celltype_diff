

library(rGREAT)



dataTypeList <- c('Down','Up')
cellTypeList <- c('DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo','Peri')


dir.create('GOBP')
dir.create('GOMF')
dir.create('MousePheno')


for (dataType in dataTypeList){
    for (cellType in cellTypeList){
	print(paste0(dataType, ', ', cellType))

	peakFile <- read.table(paste0('/home/ajl1213/Projects/PD/SciAdv_Anal/Diff_cRE/1_BulkH3K27ac/BedFiles', dataType, '.', cellType,'.bed'))
	peakList <- peakFile[, 1:3]
	greatJob <- submitGreatJob(peakList, species='hg19') 

	#availableOntologies(greatJob)
	valTable <- getEnrichmentTables(greatJob, download_by='tsv', ontology=c("GO Molecular Function", "GO Biological Process", "GO Cellular Component", "Mouse Phenotype" ))

	# GOBP
	allList <- valTable[['GO Biological Process']]
	selectedList <- allList[, c('Ontology','ID','Desc', 'BinomP', 'RegionFoldEnrich','GeneFoldEnrich','ObsGenes','TotalGenes','Regions','Genes')]
	write.table(selectedList, paste0('GOBP/', dataType, '.', cellType, '.txt'), row.names=F, col.names=T, sep='\t', quote=F) 

	# GOMF
	allList <- valTable[['GO Molecular Function']]
	selectedList <- allList[, c('Ontology','ID','Desc', 'BinomP', 'RegionFoldEnrich','GeneFoldEnrich','ObsGenes','TotalGenes','Regions','Genes')]
	write.table(selectedList, paste0('GOMF/', dataType, '.', cellType, '.txt'), row.names=F, col.names=T, sep='\t', quote=F) 

	# MousePheno
	allList <- valTable[['Mouse Phenotype']]
	selectedList <- allList[, c('Ontology','ID','Desc', 'BinomP', 'RegionFoldEnrich','GeneFoldEnrich','ObsGenes','TotalGenes','Regions','Genes')]
	write.table(selectedList, paste0('MousePheno/', dataType, '.', cellType, '.txt'), row.names=F, col.names=T, sep='\t', quote=F) 

    }
}





