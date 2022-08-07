#!/home/ajl1213/anaconda2/bin/python


import os
import numpy



DIR=os.getcwd()
resolution='5kb'
gtf='/home/ajl1213/genome.info/gencode/hg19/gencode.v19.ptn_coding.lincRNA.level_1_2'
RnaCountFile='/home/ajl1213/Projects/PD/data/snRNA/03_GetRnaCount/CountByCellType/RawCount.byCellType.qq.txt'
AbcScoreDIR=DIR+'/ActivityByContactScore'
cellTypeList=['DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']

os.system('mkdir '+DIR+'/RnaEffect')


def loadGeneInfo():
    geneDict={}
    input1=open(gtf,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        ensembleID=each[0]
        geneID=each[1]
        geneDict[ensembleID]=geneID
    input1.close()

    return geneDict


def loadRnaCount():
    rnaCountDict={}
    input1=open(RnaCountFile,'r')
    all_input1=input1.readlines()
    cellTypeList=all_input1[0].strip().split('\t')[1:]
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        geneID=each[0]
        vals=numpy.array(each[1:], dtype='float')
        tmpDict={}
        index=0
        for val in vals:
            tmpDict[cellTypeList[index]]=val
            index+=1
	rnaCountDict[geneID]=tmpDict
    input1.close()


    return rnaCountDict


def loadInterVal(cellType):
    creSumDict={}
    interSumDict={}
    contactSumDict={}
    input1=open(AbcScoreDIR+'/'+cellType+'.ActivityByContactScore.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')

	ensembleID=each[0]
	geneID=each[1]
	peakCount=float(each[6])
	interCount=float(each[7])
	contactVal=float(each[8])

	if creSumDict.has_key(geneID):
	    creSumDict[geneID]+=peakCount
	else:
	    creSumDict[geneID]=peakCount
	
	if interSumDict.has_key(geneID):
	    interSumDict[geneID]+=interCount
	else:
	    interSumDict[geneID]=interCount

	if contactSumDict.has_key(geneID):
	    contactSumDict[geneID]+=contactVal
	else:
	    contactSumDict[geneID]=contactVal	
    input1.close()

    return creSumDict, interSumDict, contactSumDict


def makeData():
    geneDict=loadGeneInfo()
    rnaCountDict=loadRnaCount()

    for cellType in cellTypeList:
	print cellType

	creSumDict, interSumDict, contactSumDict=loadInterVal(cellType)

	output1=open(DIR+'/RnaEffect/'+cellType+'.ContactValuePerGene.txt','w')
	output1.write('ensembleID\tgeneID\trnaCount\tinterSum\tcreSum\ttotalContactVal\n')
	for ensembleID in geneDict:
	    geneID=geneDict[ensembleID]

	    if rnaCountDict.has_key(geneID):
		rnaCount=str(rnaCountDict[geneID][cellType])
		if interSumDict.has_key(geneID):
		    interSum=str(interSumDict[geneID])
		    creSum=str(creSumDict[geneID])
		    contactSum=str(contactSumDict[geneID])

		    newLine=[ensembleID, geneID, rnaCount, interSum, creSum, contactSum]
		    output1.write('\t'.join(newLine)+'\n')
	output1.close()






makeData()


