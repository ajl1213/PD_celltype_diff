#!/home/ajl1213/anaconda2/bin/python


import os
import numpy


DIR=os.getcwd()
inputDIR='/home/ajl1213/Projects/PD/SciAdv_Anal/CombinedTargetID/01_FindTarget/RESULT'
cellTypeList=['DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']
minAbcScore=10


##
uniqGeneList=[]
allDict={}
for cellType in cellTypeList:
    print cellType

    tmpList=[]
    input1=open(inputDIR+'/'+cellType+'.selectList.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	[ensembleID, geneID, geneLabel, PdGeneLabel, funcGeneLabel, creID, funcPeakLabel, creType, dist, corVal, contactVal, abcScore]=line.strip().split('\t')
	if float(abcScore) > minAbcScore:
	    tmpList.append(geneID)
	    if not geneID in uniqGeneList:
		uniqGeneList.append(geneID)
    input1.close()

    allDict[cellType]=tmpList


##
output1=open(DIR+'/CellSpecificity.txt','w')
output1.write('geneID\tfuncCells\tnCellType\n')
for geneID in uniqGeneList:

    cellList=[]

    for cellType in cellTypeList:
	if geneID in allDict[cellType]:
	    cellList.append(cellType)

    funcCells=';'.join(cellList)
    nCellType=str(len(cellList))

    newLine=[geneID, funcCells, nCellType]
    output1.write('\t'.join(newLine)+'\n')

output1.close()










