#!/python2


import os
import numpy
import operator



DIR=os.getcwd()
BulkDiffPeakFile='/home/ajl1213/Projects/PD/SciAdv_Data/ChIP/02_Processing/ValTable/H3K27ac.valTable.filter.txt'
CellLabelFile='/home/ajl1213/Projects/PD/data/snATAC/04_GetAtacCount/CellLabel/AtacPeak.CellLabel.txt'
DataTypeList=['Down','Up']
CellTypeList=['DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']
nMaxCellType=8

os.system('mkdir '+DIR+'/InputData')
os.system('mkdir '+DIR+'/BedFiles')



def loadBulkPeaks():
    dysPeakDict={}
    input1=open(BulkDiffPeakFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        peakID=each[0]
        peakLabel=each[7]

	if dysPeakDict.has_key(peakLabel):
	    tmp=dysPeakDict[peakLabel]
	    tmp.append(peakID)
	    dysPeakDict[peakLabel]=tmp
	else:
	    dysPeakDict[peakLabel]=[peakID]
    input1.close()

    return dysPeakDict


def loadCellLabel():
    labelDict={}
    input1=open(CellLabelFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	peakID=each[0]
	funcCell=each[1]
	labelDict[peakID]=funcCell.split(',')
    input1.close()

    return labelDict


def makeData():
    dysPeakDict=loadBulkPeaks()
    labelDict=loadCellLabel()

    for dataType in dysPeakDict:
	output1=open(DIR+'/InputData/'+dataType+'.funcCell.txt','w')
	output1.write('peakID\t'+'\t'.join(CellTypeList)+'\n')
	for peakID in dysPeakDict[dataType]:
	    if not len(labelDict[peakID]) > nMaxCellType:
		valList=[]
		for cellType in CellTypeList:
		    if cellType in labelDict[peakID]:
			binary='1'
		    else:
			binary='0'
		    valList.append(binary)
		output1.write(peakID+'\t'+'\t'.join(valList)+'\n')
	output1.close()


def makeBED():
    dysPeakDict=loadBulkPeaks()
    labelDict=loadCellLabel()

    ##
    for dataType in dysPeakDict:
	for cellType in CellTypeList:
	    output1=open(DIR+'/BedFiles/'+dataType+'.'+cellType+'.bed','w')
	    for peakID in dysPeakDict[dataType]:
		if not len(labelDict[peakID]) > nMaxCellType:
		    if cellType in labelDict[peakID]:

			chrID=peakID.split(':')[0]
			pt1=peakID.split(':')[1].split('-')[0]
			pt2=peakID.split(':')[1].split('-')[1]


			newLine=[chrID, pt1, pt2, peakID]

			output1.write('\t'.join(newLine)+'\n')
	    output1.close()

    ##
    for dataType in dysPeakDict:
	output1=open(DIR+'/BedFiles/'+dataType+'.Merged.bed','w')
	for peakID in dysPeakDict[dataType]:
	    if not len(labelDict[peakID]) > nMaxCellType:

		chrID=peakID.split(':')[0]
		pt1=peakID.split(':')[1].split('-')[0]
		pt2=peakID.split(':')[1].split('-')[1]


		newLine=[chrID, pt1, pt2, peakID]

		output1.write('\t'.join(newLine)+'\n')
	output1.close()




makeData()
makeBED()




