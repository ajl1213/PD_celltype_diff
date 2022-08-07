#!python2


import os
import numpy
import operator


DIR=os.getcwd()
snRNA_file='/home/ajl1213/Projects/PD/SciAdv_Data/snRNAseq/02_SeuratProcess/CellTypeDEGs/snRNA.cellTypeDEGs.txt'
cellTypeList=['DopaN','GabaN','Oligo', 'OPC', 'Ast', 'Micro','Endo','Peri']

minPct = 0.2
minLog2fc = 0.2
qval_thresh = 0.05

os.system('mkdir '+DIR+'/DegMatrix')


def load_rna():
    valDict={}
    degDict={}
    input1=open(snRNA_file,'r')
    all_input1=input1.readlines()
    for cellType in cellTypeList:
	tmp1Dict={}
	tmp2Dict={}
	for line in all_input1[1:]:
	    each=line.strip().split()

	    pct1=float(each[2])
	    pct2=float(each[3])
	    log2fc=float(each[1])
	    adjPval=float(each[4])
	    cluster=each[6]
	    geneID=each[7]

	    if cluster==cellType:
		if max(pct1, pct2) > minPct:
		    tmp1Dict[geneID]=[log2fc, adjPval]

		    if (adjPval < qval_thresh and abs(log2fc) > minLog2fc):
			if log2fc < 0:
			    degLabel='Down'
			else:
			    degLabel='Up'

			if tmp2Dict.has_key(degLabel):
			    tmp=tmp2Dict[degLabel]
			    tmp.append(geneID)
			    tmp2Dict[degLabel]=tmp
			else:
			    tmp2Dict[degLabel]=[geneID]
		else: 
		    tmp1Dict[geneID]=[0, 0]

	valDict[cellType]=tmp1Dict
	degDict[cellType]=tmp2Dict

    input1.close()

    return degDict, valDict


def count_DEG():
    degDict, valDict = load_rna()

    # Count table
    output1=open(DIR+'/CountTable_snRNA.txt','w')
    output1.write('cellType\tdownDEG\tupDEG\n')
    downList=[]
    upList=[]
    for cellType in cellTypeList:
	tmpList=[]
	for dataType in ['Down','Up']:
	    if degDict[cellType].has_key(dataType):
		nDEG = len(degDict[cellType][dataType])
	    else:
		nDEG = 0
	    tmpList.append(str(nDEG))

	    if dataType=='Down':
		if degDict[cellType].has_key(dataType):
		    for geneID in degDict[cellType][dataType]:
			downList.append(geneID)
	    if dataType=='Up':
		if degDict[cellType].has_key(dataType):
		    for geneID in degDict[cellType][dataType]:
			upList.append(geneID)
	
	output1.write(cellType+'\t'+'\t'.join(tmpList)+'\n')

    downList=set(downList)
    upList=set(upList)

    output1.write('Unique\t'+str(len(downList))+'\t'+str(len(upList))+'\n')
    output1.close()

    # DEG list
    output1=open(DIR+'/snRNA.DEG_list.txt','w')
    output1.write('cellType\tDEG_type\tgeneID\tlog2fc\tadjPval\n')
    for cellType in cellTypeList:
	for dataType in ['Down','Up']:
	    if degDict[cellType].has_key(dataType):
		for geneID in degDict[cellType][dataType]:
		    log2fc=str(valDict[cellType][geneID][0])
		    adjPval=str(valDict[cellType][geneID][1])

		    newLine=[cellType, dataType, geneID, log2fc, adjPval]
		    output1.write('\t'.join(newLine)+'\n')
    output1.close()

    return downList, upList


def makeDegMatrix():
    degDict, valDict = load_rna()
    downList, upList = count_DEG()
    
    output1=open(DIR+'/DegMatrix/DegMatrix.Down.txt','w')
    output1.write('geneID\t'+'\t'.join(cellTypeList)+'\n')
    for geneID in downList:
	valList=[]
	for cellType in cellTypeList:
	    idx = cellTypeList.index(cellType)+1
	    if degDict[cellType].has_key('Down'):
		if geneID in degDict[cellType]['Down']:
		    valList.append(str(idx))
		else:
		    valList.append(str(0.0))
	    else:
		valList.append(str(0.0))
	output1.write(geneID+'\t'+'\t'.join(valList)+'\n')
    output1.close()

    output1=open(DIR+'/DegMatrix/DegMatrix.Up.txt','w')
    output1.write('geneID\t'+'\t'.join(cellTypeList)+'\n')
    for geneID in upList:
	valList=[]
	for cellType in cellTypeList:
	    idx = cellTypeList.index(cellType)+1
	    if degDict[cellType].has_key('Up'):
		if geneID in degDict[cellType]['Up']:
		    valList.append(str(idx))
		else:
		    valList.append(str(0.0))
	    else:
		valList.append(str(0.0))
	output1.write(geneID+'\t'+'\t'.join(valList)+'\n')
    output1.close()


def makeLog2Matrix():
    degDict, valDict = load_rna()
    downList, upList = count_DEG()

    ##
    output1=open(DIR+'/DegMatrix/DegMatrix.log2fc.Down.txt','w')
    output1.write('geneID\t'+'\t'.join(cellTypeList)+'\n')
    for geneID in downList:
	valList=[]
	for cellType in cellTypeList:
	    log2fc = str(valDict[cellType][geneID][0])
	    valList.append(log2fc)

	output1.write(geneID+'\t'+'\t'.join(valList)+'\n')
    output1.close()

    ##
    output1=open(DIR+'/DegMatrix/DegMatrix.log2fc.Up.txt','w')
    output1.write('geneID\t'+'\t'.join(cellTypeList)+'\n')
    for geneID in upList:
	valList=[]
	for cellType in cellTypeList:
	    log2fc = str(valDict[cellType][geneID][0])
	    valList.append(log2fc)

	output1.write(geneID+'\t'+'\t'.join(valList)+'\n')
    output1.close()


def getOrder():
    degDict, valDict = load_rna()
    downList, upList = count_DEG()

    ## Down
    allDict={}
    for geneID in downList:
	tmpList=[]
	for cellType in cellTypeList:
	    if degDict[cellType].has_key('Down'):
		if geneID in degDict[cellType]['Down']:
		    tmpList.append(cellType)
	allDict[geneID]=tmpList

    assignDict={}
    for geneID in allDict:
	if len(allDict[geneID]) == 1:
	    assignDict[geneID]=allDict[geneID][0]
	else:
	    finalCellType = 'none'
	    cur_val = 0
	    for cellType in allDict[geneID]:
		fcVal=abs(valDict[cellType][geneID][0])
		if fcVal > cur_val:
		    finalCellType = cellType
		    cur_val=fcVal
	    assignDict[geneID]=finalCellType

    output1=open(DIR+'/DegOrdered.Down.txt','w')
    output1.write('geneID\torder\n')
    idx=1
    for cellType in cellTypeList:
	tmpDict={}
	for geneID in assignDict:
	    if cellType == assignDict[geneID]:
		tmpDict[geneID]=valDict[cellType][geneID][0]

	for key in sorted(tmpDict.items(), key=operator.itemgetter(1)):
	    geneID=key[0]
	    output1.write(geneID+'\t'+str(idx)+'\n')
	    idx+=1
    output1.close()
	    
    ## Up
    allDict={}
    for geneID in upList:
	tmpList=[]
	for cellType in cellTypeList:
	    if degDict[cellType].has_key('Up'):
		if geneID in degDict[cellType]['Up']:
		    tmpList.append(cellType)
	allDict[geneID]=tmpList

    assignDict={}
    for geneID in allDict:
	if len(allDict[geneID]) == 1:
	    assignDict[geneID]=allDict[geneID][0]
	else:
	    finalCellType = 'none'
	    cur_val = 0
	    for cellType in allDict[geneID]:
		fcVal=abs(valDict[cellType][geneID][0])
		if fcVal > cur_val:
		    finalCellType = cellType
		    cur_val=fcVal
	    assignDict[geneID]=finalCellType

    output1=open(DIR+'/DegOrdered.Up.txt','w')
    output1.write('geneID\torder\n')
    idx=1
    for cellType in cellTypeList:
	tmpDict={}
	for geneID in assignDict:
	    if cellType == assignDict[geneID]:
		tmpDict[geneID]=valDict[cellType][geneID][0]

	for key in sorted(tmpDict.items(), key=operator.itemgetter(1), reverse=True):
	    geneID=key[0]
	    output1.write(geneID+'\t'+str(idx)+'\n')
	    idx+=1
    output1.close()




count_DEG()
makeDegMatrix()
makeLog2Matrix()
getOrder()




