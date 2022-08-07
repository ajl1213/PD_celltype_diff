#!/home/ajl1213/anaconda2/bin/python


import os
import numpy


DIR=os.getcwd()
GTF='/home/ajl1213/genome.info/gencode/hg19/gencode.v19.ptn_coding.lincRNA.level_1_2'
TargetGeneDIR=DIR+'/MergedTargetList'
#
GeneLabelDIR='/home/ajl1213/Projects/PD/SciAdv_Anal/Diff_RNA/4_DEG_overlap/'
FuncGeneFile='/home/ajl1213/Projects/PD/data/snRNA/03_GetRnaCount/CellLabel/GeneExp.CellLabel.txt'
GeneCountFile='/home/ajl1213/Projects/PD/data/RNAseq/2_Cellularity/AdjustedBulk/RNA.cellAdjusted.combat.qq.txt'
#
FuncPeakFile='/home/ajl1213/Projects/PD/data/snATAC/04_GetAtacCount/CellLabel/AtacPeak.CellLabel.txt'
H3K27acCountFile='/home/ajl1213/Projects/PD/SciAdv_Data/ChIP/02_Processing/AdjustedBulk/H3K27ac.cellAdjusted.combat.qq.txt'
#
AbcScoreDIR='/home/ajl1213/Projects/PD/SciAdv_Data/HiC/ActivityByContact/ActivityByContactScore'
#
dataTypeList=['Down','Up']
cellTypeList=['DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']
minAbcScore=10
minCor=0.3

os.system('mkdir '+DIR+'/RESULT')



def loadGeneInfo():
    geneDict={}
    input1=open(GTF,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        ensembleID=each[0]
        geneID=each[1]
        geneDict[ensembleID]=geneID
    input1.close()

    return geneDict


def loadGeneLabel():
    #   
    geneLabelDict={}
    for dataType in dataTypeList:
        for cellType in cellTypeList:
            keyID=dataType+';'+cellType

            input1=open(GeneLabelDIR+'/'+dataType+'.'+cellType+'.txt','r')
            all_input1=input1.readlines()
            for line in all_input1[1:]:
                each=line.strip().split('\t')
                geneID=each[0]
                geneLabel=each[1]

                if geneLabelDict.has_key(keyID):
                    tmp=geneLabelDict[keyID]
                    tmp.append(geneID)
                    geneLabelDict[keyID]=tmp
                else:
                    geneLabelDict[keyID]=[geneID]
            input1.close()

    #
    funcGeneDict={}
    input1=open(FuncGeneFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        geneID=each[0]
        cellList=each[1].split(',')
        funcGeneDict[geneID]=cellList
    input1.close()

    #
    PdGeneList=['SNCA','UCHL1','PINK1','PARK7','LRRK2',
    'ATP13A2','FBXO7','VPS35','DNAJC6','SYNJ1','CHCHD2',
    'VPS13C','GBA','MAPT','GAK','SMPD1','SCARB2','SLC17A5','ATP6V0A1','CTSB'
    ]

    return geneLabelDict, funcGeneDict, PdGeneList


def loadPeakLabel():
    #
    funcPeakDict={}
    input1=open(FuncPeakFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        peakID=each[0]
        cellList=each[1].split(',')
        funcPeakDict[peakID]=cellList
    input1.close()

    return funcPeakDict


def getCor(elementList):
    rnaCountDict={}
    input1=open(GeneCountFile,'r')
    all_input1=input1.readlines()
    rnaSamples=all_input1[0].strip().split('\t')[1:]
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	geneID=each[0]
	vals=numpy.array(each[1:], dtype='float')
	tmpDict={}
	idx=0
	for val in vals:
	    tmpDict[rnaSamples[idx]]=val
	    idx+=1
	rnaCountDict[geneID]=tmpDict
    input1.close()

    h3k27acCountDict={}
    input1=open(H3K27acCountFile,'r')
    all_input1=input1.readlines()
    h3k27acSamples=all_input1[0].strip().split('\t')[1:]
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	peakID=each[0]
	vals=numpy.array(each[1:], dtype='float')
	tmpDict={}
	idx=0
	for val in vals:
	    tmpDict[h3k27acSamples[idx]]=val
	    idx+=1
	h3k27acCountDict[peakID]=tmpDict
    input1.close()
    h3k27acCountDict['chr7:23144778-23146888']=h3k27acCountDict['chr7:23155188-23156134']

    #
    print 'getting correlation values..'
    commonSamples= set(rnaSamples) & set(h3k27acSamples)
    corDict={}
    idx=1
    for i in set(elementList):
	print str(idx)+' / '+str(len(set(elementList)))
	idx+=1
	
	geneID=i.split(';')[0]
	peakID=i.split(';')[1]

	if rnaCountDict.has_key(geneID):
	    rnaValList=[]
	    h3k27acValList=[]

	    for sampleID in commonSamples:
		rnaVal= rnaCountDict[geneID][sampleID]
		h3k27acVal= h3k27acCountDict[peakID][sampleID]

		rnaValList.append(rnaVal)
		h3k27acValList.append(h3k27acVal)

	    rnaZscoreList=[]
	    meanRnaVal=numpy.mean(rnaValList)
	    stdRnaVal=numpy.std(rnaValList)
	    for j in rnaValList:
		zVal=(j-meanRnaVal)/stdRnaVal
		rnaZscoreList.append(zVal)

	    h3k27acZscoreList=[]
	    meanRnaVal=numpy.mean(h3k27acValList)
	    stdRnaVal=numpy.std(h3k27acValList)
	    for j in h3k27acValList:
		zVal=(j-meanRnaVal)/stdRnaVal
		h3k27acZscoreList.append(zVal)

	    corVal=str(numpy.corrcoef(rnaZscoreList, h3k27acZscoreList)[0,1])
	    keyID=geneID+';'+peakID
	    corDict[keyID]=corVal

    print 'done'

    return corDict


def loadContactScore():
    print 'loading contact scores..'
    contactValDict={}

    for cellType in cellTypeList:
        tmpDict={}

        input1=open(AbcScoreDIR+'/'+cellType+'.ActivityByContactScore.txt','r')
        all_input1=input1.readlines()
        for line in all_input1[1:]:
            each=line.strip().split('\t')
            ensembleID=each[0]
            creID=each[3]
            contactVal=float(each[8])
            abcScore=float(each[9])

            keyID=ensembleID+';'+creID
            tmpDict[keyID]=[contactVal, abcScore]

        contactValDict[cellType]=tmpDict

        input1.close()

    print 'done'

    return contactValDict


def attach_info():
    geneLabelDict, funcGeneDict, PdGeneList=loadGeneLabel()
    funcPeakDict=loadPeakLabel()
    contactValDict=loadContactScore()

    for cellType in cellTypeList:
	print cellType
	#
	elementList=[]
	input1=open(DIR+'/MergedTargetList/'+cellType+'.MergedTargetList.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    ensembleID=each[0]
	    geneID=each[1]
	    creID=each[2]

	    valID=geneID+';'+creID
	    elementList.append(valID)
	input1.close()

	corDict=getCor(elementList)

	#
	input1=open(DIR+'/MergedTargetList/'+cellType+'.MergedTargetList.txt','r')
	all_input1=input1.readlines()
	output1=open(DIR+'/RESULT/'+cellType+'.fullList.txt','w')
	output1.write('ensembleID\tgeneID\tgeneLabel\tPdGeneLabel\tfuncGeneLabel\tcreID\tfuncPeakLabel\tcreType\tdist\tcorVal\tcontactVal\tabcScore\n')
	for line in all_input1[1:]:
	    each=line.strip().split('\t')

	    ensembleID=each[0]
	    geneID=each[1]
	    creID=each[2]
	    creType=each[3]
	    dist=each[4]
	
	    #
	    if not (geneID in geneLabelDict['Down;'+cellType] or geneID in geneLabelDict['Up;'+cellType]):
		geneLabel='nonDEG'
	    else:
		if geneID in geneLabelDict['Down;'+cellType]:
		    geneLabel='downDEG'

		elif geneID in geneLabelDict['Up;'+cellType]:
		    geneLabel='upDEG'
	    #
	    if funcGeneDict.has_key(geneID):
		funcGeneLabel=funcGeneDict[geneID]
	    else:
		funcGeneLabel='none'

	    #
	    if geneID in PdGeneList:
		PdGeneLabel='knownPdGene'
	    else:
		PdGeneLabel='nonPdGene'

	    #
	    funcPeakLabel=funcPeakDict[creID]
	    #
	    if corDict.has_key(geneID+';'+creID):
		corVal=corDict[geneID+';'+creID]
	    else:
		corVal=0
	    #
	    contactVal=contactValDict[cellType][ensembleID+';'+creID][0] if contactValDict[cellType].has_key(ensembleID+';'+creID) else 0
	    abcScore=contactValDict[cellType][ensembleID+';'+creID][1] if contactValDict[cellType].has_key(ensembleID+';'+creID) else 0

	    # 
	    if not funcGeneLabel == 'none':
		newLine=[ensembleID, geneID, geneLabel, PdGeneLabel, ','.join(funcGeneLabel), creID, ','.join(funcPeakLabel), creType, dist, str(corVal), str(contactVal), str(abcScore)]
		output1.write('\t'.join(newLine)+'\n')
	    else:
		newLine=[ensembleID, geneID, geneLabel, PdGeneLabel, 'none', creID, ','.join(funcPeakLabel), creType, dist, str(corVal), str(contactVal), str(abcScore)]
		output1.write('\t'.join(newLine)+'\n')


	output1.close()



def select():
    for cellType in cellTypeList:
	print cellType

	input1=open(DIR+'/RESULT/'+cellType+'.fullList.txt','r')
	all_input1=input1.readlines()
	output1=open(DIR+'/RESULT/'+cellType+'.selectList.txt','w')
	output1.write('ensembleID\tgeneID\tgeneLabel\tPdGeneLabel\tfuncGeneLabel\tcreID\tfuncPeakLabel\tcreType\tdist\tcorVal\tcontactVal\tabcScore\n')
	for line in all_input1[1:]:
	    [ensembleID, geneID, geneLabel, PdGeneLabel, funcGeneLabel, creID, funcPeakLabel, creType, dist, corVal, contactVal, abcScore]=line.strip().split('\t')
	    if (cellType in funcGeneLabel.split(',') and cellType in funcPeakLabel.split(',')):
		if float(corVal) > minCor:
		    if float(abcScore) > minAbcScore:
			newLine=[ensembleID, geneID, geneLabel, PdGeneLabel, funcGeneLabel, creID, funcPeakLabel, creType, dist, corVal, contactVal, abcScore]
			output1.write('\t'.join(newLine)+'\n')

			if float(abcScore) > minAbcScore:
			    if PdGeneLabel == 'knownPdGene':
				print '\t'.join(newLine)				

	input1.close()
	output1.close()



attach_info()
select()



