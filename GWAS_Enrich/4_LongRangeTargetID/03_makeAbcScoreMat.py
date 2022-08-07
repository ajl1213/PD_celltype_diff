
import os
import numpy
import operator


DIR=os.getcwd()
inputDIR=DIR+'/ContactScoreList'
NearbyTargetFile=DIR+'/NearbyTargetList/MergedPD.NearbyTargetList.5kb.txt'
LongRangeTargetFile=DIR+'/LongRangeTargetList/MergedPD.LongRangeTargetList.txt'

cellTypeList=['DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']
minScore=20


def loadTargetType():

    nearbyTargetList=[]
    input1=open(NearbyTargetFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	geneID=each[1] 
	nearbyTargetList.append(geneID)
    input1.close()   

    longRangeTargetList=[]
    input1=open(LongRangeTargetFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	geneID=each[1] 
	longRangeTargetList.append(geneID)
    input1.close()   

    return nearbyTargetList, longRangeTargetList


def makeValMat():
    nearbyTargetList, longRangeTargetList = loadTargetType()

    valDict={}
    geneList=[]
    for cellType in cellTypeList:

	tmpDict={}
	input1=open(inputDIR+'/MergedPD.'+cellType+'.ContactScore.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    geneID=each[1]
	    if not geneID in geneList:
		geneList.append(geneID)

	    abcVal=float(each[4])
	    tmpDict[geneID]=abcVal
	valDict[cellType]=tmpDict
	input1.close()

    ## all matrix
    output1=open(DIR+'/AbcScoreMat.all.txt','w')
    output1.write('geneID\t'+'\t'.join(cellTypeList)+'\n')
    for geneID in geneList:
	valList=[]
	for cellType in cellTypeList:
	    val=str(valDict[cellType][geneID])
	    valList.append(val)
	output1.write(geneID+'\t'+'\t'.join(valList)+'\n')
    output1.close()


    ## filter
    filteredList=[]
    for geneID in geneList:
	tmpList=[]
	for cellType in cellTypeList:
	    val=valDict[cellType][geneID]
	    tmpList.append(val)

	if max(tmpList) > minScore:
	    filteredList.append(geneID)

    ## order
    assignDict={}
    for geneID in filteredList:
	finalCellType='none'
	cur_val=0
	for cellType in cellTypeList:
	    val=float(valDict[cellType][geneID])
	    if val > cur_val:
		finalCellType=cellType
		cur_val=val
	assignDict[geneID]=finalCellType

    output1=open(DIR+'/AbcScoreMat.filtered.ordered.txt','w')
    output1.write('geneID\tmaxCell\tNearbyTarget\tLongRangeTarget\t'+'\t'.join(cellTypeList)+'\n')
    for cellType in cellTypeList:
	tmpDict={}
	for geneID in filteredList:
	    if cellType == assignDict[geneID]:
		tmpDict[geneID]=float(valDict[cellType][geneID])

	for key in sorted(tmpDict.items(), key=operator.itemgetter(1), reverse=True):
	    geneID=key[0]

	    if geneID in nearbyTargetList:
		nearbyLabel='Target'
	    else:
		nearbyLabel='nonTarget'

	    if geneID in longRangeTargetList:
		longRangeLabel='Target'
	    else:
		longRangeLabel='nonTarget'

	    valList=[]
	    for cellType in cellTypeList:
		val=str(valDict[cellType][geneID])
		valList.append(val)
	    output1.write(geneID+'\t'+assignDict[geneID]+'\t'+nearbyLabel+'\t'+longRangeLabel+'\t'+'\t'.join(valList)+'\n')
    output1.close()
	    

	
	


makeValMat()



