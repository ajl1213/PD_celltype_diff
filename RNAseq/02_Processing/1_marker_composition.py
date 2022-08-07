#!/python2


import os
import numpy
import operator


DIR=os.getcwd()
CellTypeRnaFile='/home/ajl1213/Projects/PD/data/snRNA/03_GetRnaCount/CountByCellType/RawCount.byCellType.qq.txt'
BulkRnaFile='/home/ajl1213/Projects/PD/SciAdv_Data/RNA/01_Mapping/CountMatrix/RNA.rawCount.qq.txt'
MarkerFile='/home/ajl1213/Projects/PD/data/snRNA/02a_Analysis.Original/CellTypeDEGs/snRNA.cellTypeMarkers.unique.txt'
cellTypeList=['DopaN','GabaN','Oligo','OPC','Ast','Micro','Endo','Peri']

os.system('mkdir '+DIR+'/MarkerGeneRatio')


def getGeneList():
    # 
    markerDict={}
    input1=open(MarkerFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	geneID=each[6]
	cellType=each[5]
	markerDict[geneID]=cellType
    input1.close()

    #
    bulkGeneList=[]
    input1=open(BulkRnaFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	ensembleID=each[0]
	geneID=each[1]
	bulkGeneList.append(geneID)
    input1.close()

    snRnaGeneList=[]
    input1=open(CellTypeRnaFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	geneID=each[0]
	snRnaGeneList.append(geneID)
    input1.close()

    commonGeneSet=set(bulkGeneList)&set(snRnaGeneList)

    #
    output1=open(DIR+'/MarkerGeneRatio/geneList.txt','w')
    output1.write('geneID\tgeneLabel\n')
    for geneID in commonGeneSet:
	if markerDict.has_key(geneID):
	    geneLabel=markerDict[geneID]
	else:
	    geneLabel='none'
	output1.write(geneID+'\t'+geneLabel+'\n')
    output1.close()


def loadBulkSample():
    #
    geneList=[]
    geneLabelDict={}
    input1=open(DIR+'/MarkerGeneRatio/geneList.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	geneID=each[0]
	geneLabel=each[1]
	geneList.append(geneID)	
	geneLabelDict[geneID]=geneLabel
    input1.close()

    #
    bulkDict={}
    input1=open(BulkRnaFile,'r')
    all_input1=input1.readlines()
    sampleList=all_input1[0].strip().split('\t')[2:]
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	ensembleID=each[0]
	geneID=each[1]
	if geneID in geneList:
	    vals=numpy.array(each[2:], dtype='float')
	    tmpDict={}
	    index=0
	    for val in vals:
		tmpDict[sampleList[index]]=val
		index+=1
	    bulkDict[geneID]=tmpDict
	input1.close()

    #
    allDict={}
    for sampleID in sampleList:
	sumDict={}
	totalSum=0
	for geneID in geneList:
	    val = float(bulkDict[geneID][sampleID])
	    geneLabel = geneLabelDict[geneID]
	    if not geneLabel =='none':
		totalSum+=val
		if sumDict.has_key(geneLabel):
		    sumDict[geneLabel]+=val
		else:
		    sumDict[geneLabel]=val
	ratioDict={}
	for geneLabel in sumDict:
	    ratioDict[geneLabel]=sumDict[geneLabel]/totalSum
	allDict[sampleID]=ratioDict

    #
    output1=open(DIR+'/MarkerGeneRatio/marker_composition.txt','w')
    output1.write('sampleID\tcellType\tratio\n')
    for sampleID in sampleList:
	for cellType in cellTypeList:
	    ratioVal=allDict[sampleID][cellType]
	    newLine=[sampleID, cellType, str(ratioVal)]
	    output1.write('\t'.join(newLine)+'\n')
    output1.close()

    output1=open(DIR+'/MarkerGeneRatio/avg_marker_composition.txt','w')
    output1.write('cellType\tavgRatio\n')
    for cellType in cellTypeList:
	valList=[]
	for sampleID in sampleList:
	    ratioVal=allDict[sampleID][cellType]
	    valList.append(ratioVal)
	meanVal=numpy.mean(valList)
	output1.write(cellType+'\t'+str(meanVal)+'\n')
    output1.close()




getGeneList()
loadBulkSample()



