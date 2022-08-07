#!/python2


import os
import numpy
import operator


DIR=os.getcwd()
CellTypeAtacFile='/home/ajl1213/Projects/PD/data/snATAC/04_GetAtacCount/PeakCount/atacCount.byCellType.qq.txt'
BulkH3K27acFile=DIR+'/CountMatrix/H3K27ac.readCount.denoised.qq.txt'
TotalPeakBED='/home/ajl1213/Projects/PD/data/snATAC/04_GetAtacCount/CellLabel/AtacPeak.TOTAL.bed'
MarkerFile='/home/ajl1213/Projects/PD/data/snATAC/04_GetAtacCount/CellLabel/MarkerPeaks.ratio.txt'
cellTypeList=['DopaN','GabaN','Oligo','OPC','Ast','Micro','Endo','Peri']

os.system('mkdir '+DIR+'/MarkerPeakRatio')


def getPeakList():
    #
    markerDict={}
    input1=open(MarkerFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        peakID=each[0]
        cellType=each[9]
        markerDict[peakID]=cellType
    input1.close()

    #
    peakList=[]
    input1=open(TotalPeakBED,'r')
    all_input1=input1.readlines()
    for line in all_input1:
	each=line.strip().split('\t')
	chrID=each[0]
	pt1=each[1]
	pt2=each[2]
	peakID=chrID+':'+pt1+'-'+pt2
	peakList.append(peakID)
    input1.close()

    #
    output1=open(DIR+'/MarkerPeakRatio/peakList.txt','w')
    output1.write('peakID\tpeakLabel\n')
    for peakID in peakList:
	if markerDict.has_key(peakID):
	    peakLabel=markerDict[peakID]
	else:
	    peakLabel='none'
	output1.write(peakID+'\t'+peakLabel+'\n')
    output1.close()


def loadBulkSample():
    #   
    peakList=[]
    peakLabelDict={}
    input1=open(DIR+'/MarkerPeakRatio/peakList.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        peakID=each[0]
        peakLabel=each[1]
        peakList.append(peakID) 
        peakLabelDict[peakID]=peakLabel
    input1.close()

    # H3K27ac
    bulkDict={}
    input1=open(BulkH3K27acFile,'r')
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
	    bulkDict[peakID]=tmpDict
	input1.close()

    allDict={}
    for sampleID in sampleList:
        sumDict={}
        totalSum=0
        for peakID in peakList:
            val = float(bulkDict[peakID][sampleID])
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

    output1=open(DIR+'/MarkerPeakRatio/markerComposition.H3K27ac.txt','w')
    output1.write('sampleID\tcellType\tratio\n')
    for sampleID in sampleList:
        for cellType in cellTypeList:
            ratioVal=allDict[sampleID][cellType]
	    newLine=[sampleID, cellType, str(ratioVal)]
	    output1.write('\t'.join(newLine)+'\n')
    output1.close()

    output1=open(DIR+'/MarkerPeakRatio/avgMarkerComposition.H3K27ac.txt','w')
    output1.write('cellType\tavgRatio\n')
    for cellType in cellTypeList:
	valList=[]
	for sampleID in sampleList:
	    ratioVal=allDict[sampleID][cellType]
	    valList.append(ratioVal)
	meanVal=numpy.mean(valList)
	output1.write(cellType+'\t'+str(meanVal)+'\n')
    output1.close()




getPeakList()
loadBulkSample()



