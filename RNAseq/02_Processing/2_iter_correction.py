#!/python2


import os
import numpy
import operator


DIR=os.getcwd()
BulkRnaFile='/home/ajl1213/Projects/PD/SciAdv_Data/RNA/01_Mapping/CountMatrix/RNA.rawCount.qq.txt'
RefRatioFile='/home/ajl1213/Projects/PD/data/snRNA/03_GetRnaCount/CountByCellType/ExpRatio.byCellType.txt'
GeneListFile=DIR+'/MarkerGeneRatio/geneList.txt'
AvgMarkerCompositionFile=DIR+'/MarkerGeneRatio/avg_marker_composition.txt'
MarkerCompositionFile=DIR+'/MarkerGeneRatio/marker_composition.txt'
maxIter=3

os.system('mkdir '+DIR+'/AdjustedBulk/')


def loadBulkRNA():
    valDict={}
    input1=open(BulkRnaFile,'r')
    all_input1=input1.readlines()
    sampleList=all_input1[0].strip().split('\t')[2:]
    for line in all_input1[1:]:
        each=line.strip().split('\t')
	ensembleID=each[0]
        geneID=each[1]
        vals=numpy.array(each[2:], dtype='float')
        tmpDict={}
        index=0
        for val in vals:
            tmpDict[sampleList[index]]=val
            index+=1
        valDict[geneID]=tmpDict
    input1.close()

    return sampleList, valDict


def loadRefRatio():
    refRatioDict={}
    input1=open(RefRatioFile,'r')
    all_input1=input1.readlines()
    cellTypeList=all_input1[0].strip().split('\t')[1:]
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        geneID=each[0]
        fracVals=numpy.array(each[1:], dtype='float')
        tmpDict={}
        index=0
        for fracVal in fracVals:
            tmpDict[cellTypeList[index]]=fracVal
            index+=1
        refRatioDict[geneID]=tmpDict
    input1.close()

    return cellTypeList, refRatioDict


def loadMarkerComposition():
    avgMarkerCompDict={}
    input1=open(AvgMarkerCompositionFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	cellType=each[0]
	avgRatio=float(each[1])
	avgMarkerCompDict[cellType]=avgRatio
    input1.close()

    markerCompDict={}
    input1=open(MarkerCompositionFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	sampleID=each[0]
	cellType=each[1]
	ratioVal=float(each[2])
	markerCompDict[sampleID+';'+cellType]=ratioVal
    input1.close()

    return avgMarkerCompDict, markerCompDict

    
def adjustCellularity():
    sampleList, valDict = loadBulkRNA()
    cellTypeList, refRatioDict = loadRefRatio()
    avgMarkerCompDict, markerCompDict = loadMarkerComposition()

    # gene list
    input1=open(GeneListFile,'r')
    all_input1=input1.readlines()
    geneList=[]
    geneLabelDict={}
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	geneID=each[0]
	geneLabel=each[1]
	geneList.append(geneID)
	geneLabelDict[geneID]=geneLabel
    input1.close()

    for i in range(maxIter):
	nIter = i+1
	print str(nIter)

	if nIter == 1:
	    # correction
	    output1=open(DIR+'/AdjustedBulk/RNA.cellAdjusted.iter_'+str(nIter)+'.txt','w')
	    output1.write('geneID\t'+'\t'.join(sampleList)+'\n')
	    for geneID in geneList:
		valList=[]
		for sampleID in sampleList:
		    finalVal=0
		    for cellType in cellTypeList:
			val=valDict[geneID][sampleID]
			adjVal=val * (refRatioDict[geneID][cellType]/100) * (avgMarkerCompDict[cellType]/markerCompDict[sampleID+';'+cellType])
			finalVal+=adjVal
		    valList.append(str(finalVal))

		output1.write(geneID+'\t'+'\t'.join(valList)+'\n')
	    output1.close()

	    # compute composition
	    adjustedDict={}
	    input1=open(DIR+'/AdjustedBulk/RNA.cellAdjusted.iter_'+str(nIter)+'.txt','r')
	    all_input1=input1.readlines()
	    sampleList=all_input1[0].strip().split('\t')[1:]
	    for line in all_input1[1:]:
		each=line.strip().split('\t')
		geneID=each[0]
		if geneID in geneList:
		    vals=numpy.array(each[1:], dtype='float')
		    tmpDict={}
		    index=0
		    for val in vals:
			tmpDict[sampleList[index]]=val
			index+=1
		    adjustedDict[geneID]=tmpDict
		input1.close()

	    allDict={}
	    for sampleID in sampleList:
		sumDict={}
		totalSum=0
		for geneID in geneList:
		    val = float(adjustedDict[geneID][sampleID])
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

	    output1=open(DIR+'/AdjustedBulk/marker_composition.iter_'+str(nIter)+'.txt','w')
	    output1.write('sampleID\tcellType\tratio\n')
	    for sampleID in sampleList:
		for cellType in cellTypeList:
		    ratioVal=allDict[sampleID][cellType]
		    newLine=[sampleID, cellType, str(ratioVal)]
		    output1.write('\t'.join(newLine)+'\n')
	    output1.close()

	if not nIter == 1:
	    valDict={}
	    markerCompDict={}
	    
	    # load from previous correction
	    input1=open(DIR+'/AdjustedBulk/RNA.cellAdjusted.iter_'+str(nIter-1)+'.txt','r')
	    all_input1=input1.readlines()
	    sampleList=all_input1[0].strip().split('\t')[1:]
	    for line in all_input1[1:]:
		each=line.strip().split('\t')
		geneID=each[0]
		vals=numpy.array(each[1:], dtype='float')
		tmpDict={}
		index=0
		for val in vals:
		    tmpDict[sampleList[index]]=val
		    index+=1
		valDict[geneID]=tmpDict
	    input1.close()

	    input1=open(DIR+'/AdjustedBulk/marker_composition.iter_'+str(nIter-1)+'.txt','r')
	    all_input1=input1.readlines()
	    for line in all_input1[1:]:
		each=line.strip().split('\t')
		sampleID=each[0]
		cellType=each[1]
		ratioVal=float(each[2])
		markerCompDict[sampleID+';'+cellType]=ratioVal
	    input1.close()

	    # correction
	    output1=open(DIR+'/AdjustedBulk/RNA.cellAdjusted.iter_'+str(nIter)+'.txt','w')
	    output1.write('geneID\t'+'\t'.join(sampleList)+'\n')
	    for geneID in geneList:
		valList=[]
		for sampleID in sampleList:
		    finalVal=0
		    for cellType in cellTypeList:
			val=valDict[geneID][sampleID]
			adjVal=val * (refRatioDict[geneID][cellType]/100) * (avgMarkerCompDict[cellType]/markerCompDict[sampleID+';'+cellType])
			finalVal+=adjVal
		    valList.append(str(finalVal))

		output1.write(geneID+'\t'+'\t'.join(valList)+'\n')
	    output1.close()

	    # compute composition
	    adjustedDict={}
	    input1=open(DIR+'/AdjustedBulk/RNA.cellAdjusted.iter_'+str(nIter)+'.txt','r')
	    all_input1=input1.readlines()
	    sampleList=all_input1[0].strip().split('\t')[1:]
	    for line in all_input1[1:]:
		each=line.strip().split('\t')
		geneID=each[0]
		if geneID in geneList:
		    vals=numpy.array(each[1:], dtype='float')
		    tmpDict={}
		    index=0
		    for val in vals:
			tmpDict[sampleList[index]]=val
			index+=1
		    adjustedDict[geneID]=tmpDict
		input1.close()

	    allDict={}
	    for sampleID in sampleList:
		sumDict={}
		totalSum=0
		for geneID in geneList:
		    val = float(adjustedDict[geneID][sampleID])
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

	    output1=open(DIR+'/AdjustedBulk/marker_composition.iter_'+str(nIter)+'.txt','w')
	    output1.write('sampleID\tcellType\tratio\n')
	    for sampleID in sampleList:
		for cellType in cellTypeList:
		    ratioVal=allDict[sampleID][cellType]
		    newLine=[sampleID, cellType, str(ratioVal)]
		    output1.write('\t'.join(newLine)+'\n')
	    output1.close()





adjustCellularity()




