#!/python2


import os
import numpy
import operator


DIR=os.getcwd()
BulkH3K27acFile=DIR+'/CountMatrix/H3K27ac.readCount.denoised.qq.txt'
RefRatioFile='/home/ajl1213/Projects/PD/data/snATAC/04_GetAtacCount/PeakCount/atacRatio.byCellType.txt'
PeakListFile=DIR+'/MarkerPeakRatio/peakList.txt'
AvgMarkerCompositionFile=DIR+'/MarkerPeakRatio/avgMarkerComposition.H3K27ac.txt'
MarkerCompositionFile=DIR+'/MarkerPeakRatio/markerComposition.H3K27ac.txt'
maxIter=3

os.system('mkdir '+DIR+'/AdjustedBulk/')


def loadBulkH3K27ac():
    valDict={}
    input1=open(BulkH3K27acFile,'r')
    all_input1=input1.readlines()
    sampleList=all_input1[0].strip().split('\t')[1:]
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        peakID=each[0]
        vals=numpy.array(each[1:], dtype='float')
        tmpDict={}
        index=0
        for val in vals:
            tmpDict[sampleList[index]]=val
            index+=1
        valDict[peakID]=tmpDict
    input1.close()

    return sampleList, valDict


def loadRefRatio():
    refRatioDict={}
    input1=open(RefRatioFile,'r')
    all_input1=input1.readlines()
    cellTypeList=all_input1[0].strip().split('\t')[1:]
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        peakID=each[0]
        fracVals=numpy.array(each[1:], dtype='float')
        tmpDict={}
        index=0
        for fracVal in fracVals:
            tmpDict[cellTypeList[index]]=fracVal
            index+=1
        refRatioDict[peakID]=tmpDict
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
    sampleList, valDict = loadBulkH3K27ac()
    cellTypeList, refRatioDict = loadRefRatio()
    avgMarkerCompDict, markerCompDict = loadMarkerComposition()

    # peak list
    input1=open(PeakListFile,'r')
    all_input1=input1.readlines()
    peakList=[]
    peakLabelDict={}
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	peakID=each[0]
	peakLabel=each[1]
	peakList.append(peakID)
	peakLabelDict[peakID]=peakLabel
    input1.close()

    for i in range(maxIter):
	nIter = i+1
	print str(nIter)

	if nIter == 1:
	    # correction
	    output1=open(DIR+'/AdjustedBulk/H3K27ac.cellAdjusted.iter_'+str(nIter)+'.txt','w')
	    output1.write('peakID\t'+'\t'.join(sampleList)+'\n')
	    for peakID in peakList:
		valList=[]
		for sampleID in sampleList:
		    finalVal=0
		    for cellType in cellTypeList:
			val=valDict[peakID][sampleID]
			adjVal=val * (refRatioDict[peakID][cellType]/100) * (avgMarkerCompDict[cellType]/markerCompDict[sampleID+';'+cellType])
			finalVal+=adjVal
		    valList.append(str(finalVal))

		output1.write(peakID+'\t'+'\t'.join(valList)+'\n')
	    output1.close()

	    # compute composition
	    adjustedDict={}
	    input1=open(DIR+'/AdjustedBulk/H3K27ac.cellAdjusted.iter_'+str(nIter)+'.txt','r')
	    all_input1=input1.readlines()
	    sampleList=all_input1[0].strip().split('\t')[1:]
	    for line in all_input1[1:]:
		each=line.strip().split('\t')
		peakID=each[0]
		if peakID in peakList:
		    vals=numpy.array(each[1:], dtype='float')
		    tmpDict={}
		    index=0
		    for val in vals:
			tmpDict[sampleList[index]]=val
			index+=1
		    adjustedDict[peakID]=tmpDict
		input1.close()

	    allDict={}
	    for sampleID in sampleList:
		sumDict={}
		totalSum=0
		for peakID in peakList:
		    val = float(adjustedDict[peakID][sampleID])
		    peakLabel = peakLabelDict[peakID]
		    if not peakLabel =='none':
			totalSum+=val
			if sumDict.has_key(peakLabel):
			    sumDict[peakLabel]+=val
			else:
			    sumDict[peakLabel]=val
		ratioDict={}
		for peakLabel in sumDict:
		    ratioDict[peakLabel]=sumDict[peakLabel]/totalSum
		allDict[sampleID]=ratioDict

	    output1=open(DIR+'/AdjustedBulk/markerComposition.H3K27ac.iter_'+str(nIter)+'.txt','w')
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
	    input1=open(DIR+'/AdjustedBulk/H3K27ac.cellAdjusted.iter_'+str(nIter-1)+'.txt','r')
	    all_input1=input1.readlines()
	    sampleList=all_input1[0].strip().split('\t')[1:]
	    for line in all_input1[1:]:
		each=line.strip().split('\t')
		peakID=each[0]
		vals=numpy.array(each[1:], dtype='float')
		tmpDict={}
		index=0
		for val in vals:
		    tmpDict[sampleList[index]]=val
		    index+=1
		valDict[peakID]=tmpDict
	    input1.close()

	    input1=open(DIR+'/AdjustedBulk/markerComposition.H3K27ac.iter_'+str(nIter-1)+'.txt','r')
	    all_input1=input1.readlines()
	    for line in all_input1[1:]:
		each=line.strip().split('\t')
		sampleID=each[0]
		cellType=each[1]
		ratioVal=float(each[2])
		markerCompDict[sampleID+';'+cellType]=ratioVal
	    input1.close()

	    # correction
	    output1=open(DIR+'/AdjustedBulk/H3K27ac.cellAdjusted.iter_'+str(nIter)+'.txt','w')
	    output1.write('peakID\t'+'\t'.join(sampleList)+'\n')
	    for peakID in peakList:
		valList=[]
		for sampleID in sampleList:
		    finalVal=0
		    for cellType in cellTypeList:
			val=valDict[peakID][sampleID]
			adjVal=val * (refRatioDict[peakID][cellType]/100) * (avgMarkerCompDict[cellType]/markerCompDict[sampleID+';'+cellType])
			finalVal+=adjVal
		    valList.append(str(finalVal))

		output1.write(peakID+'\t'+'\t'.join(valList)+'\n')
	    output1.close()

	    # compute composition
	    adjustedDict={}
	    input1=open(DIR+'/AdjustedBulk/H3K27ac.cellAdjusted.iter_'+str(nIter)+'.txt','r')
	    all_input1=input1.readlines()
	    sampleList=all_input1[0].strip().split('\t')[1:]
	    for line in all_input1[1:]:
		each=line.strip().split('\t')
		peakID=each[0]
		if peakID in peakList:
		    vals=numpy.array(each[1:], dtype='float')
		    tmpDict={}
		    index=0
		    for val in vals:
			tmpDict[sampleList[index]]=val
			index+=1
		    adjustedDict[peakID]=tmpDict
		input1.close()

	    allDict={}
	    for sampleID in sampleList:
		sumDict={}
		totalSum=0
		for peakID in peakList:
		    val = float(adjustedDict[peakID][sampleID])
		    peakLabel = peakLabelDict[peakID]
		    if not peakLabel =='none':
			totalSum+=val
			if sumDict.has_key(peakLabel):
			    sumDict[peakLabel]+=val
			else:
			    sumDict[peakLabel]=val
		ratioDict={}
		for peakLabel in sumDict:
		    ratioDict[peakLabel]=sumDict[peakLabel]/totalSum
		allDict[sampleID]=ratioDict

	    output1=open(DIR+'/AdjustedBulk/markerComposition.H3K27ac.iter_'+str(nIter)+'.txt','w')
	    output1.write('sampleID\tcellType\tratio\n')
	    for sampleID in sampleList:
		for cellType in cellTypeList:
		    ratioVal=allDict[sampleID][cellType]
		    newLine=[sampleID, cellType, str(ratioVal)]
		    output1.write('\t'.join(newLine)+'\n')
	    output1.close()



adjustCellularity()



