#!/home/ajl1213/anaconda2/bin/python


import os
import numpy


DIR=os.getcwd()
resolution='5kb'
NearbyWindow=[5000,'5kb']
GTF='/home/ajl1213/genome.info/gencode/hg19/gencode.v19.ptn_coding.lincRNA.level_1_2'
gwasID='MergedPD'

NearbyTargetDIR=DIR+'/NearbyTargetList'
LongRangeTargetDIR=DIR+'/LongRangeTargetList'
AbcScoreDIR='/home/ajl1213/Projects/PD/data/HiC/ActivityByContact/ActivityByContactScore'
cellTypeList=['DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']

os.system('mkdir '+DIR+'/ContactScoreList')



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


def loadContactScore(cellType):
    contactValDict={}
    input1=open(AbcScoreDIR+'/'+cellType+'.ActivityByContactScore.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	ensembleID=each[0]
	creID=each[3]
	contactVal=float(each[8])
	abcScore=float(each[9])
	
	keyID=ensembleID+';'+creID
	valID=[contactVal, abcScore]

	contactValDict[keyID]=valID
    input1.close()

    return contactValDict


def arrangeTargetList():
    geneDict=loadGeneInfo()

    for cellType in cellTypeList:
	contactValDict=loadContactScore(cellType)

	allDict={}
	## load nearby cREs
	input1=open(NearbyTargetDIR+'/'+gwasID+'.NearbyTargetList.'+NearbyWindow[1]+'.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    ensembleID=each[0]
	    geneID=each[1]
	    creList=each[2].split(',')
	    for creID in creList:
		if allDict.has_key(ensembleID):
		    tmp=allDict[ensembleID]
		    tmp.append(creID)
		    allDict[ensembleID]=tmp
		else:
		    allDict[ensembleID]=[creID]
	input1.close()

	## load long-range cREs
	input1=open(LongRangeTargetDIR+'/'+gwasID+'.LongRangeTargetList.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    ensembleID=each[0]
	    geneID=each[1]
	    creList=each[2].split(',')
	    for creID in creList:
		if allDict.has_key(ensembleID):
		    tmp=allDict[ensembleID]
		    tmp.append(creID)
		    allDict[ensembleID]=tmp
		else:
		    allDict[ensembleID]=[creID]
	input1.close()

	##
	output1=open(DIR+'/ContactScoreList/'+gwasID+'.'+cellType+'.ContactScore.txt','w')
	output1.write('ensembleID\tgeneID\tcreList\tsumContactVal\tsumAbcScore\n')
	for ensembleID in allDict:
	    geneID=geneDict[ensembleID]

	    sumContactVal=0
	    sumAbcScore=0
	    for creID in set(allDict[ensembleID]):
		if contactValDict.has_key(ensembleID+';'+creID):
		    contactVal=contactValDict[ensembleID+';'+creID][0]
		    abcScore=contactValDict[ensembleID+';'+creID][1]

		    sumContactVal+=contactVal
		    sumAbcScore+=abcScore

	    creList=','.join(set(allDict[ensembleID]))
	    newLine=[ensembleID, geneID, creList, str(sumContactVal), str(sumAbcScore)]
	    output1.write('\t'.join(newLine)+'\n')
	output1.close()



arrangeTargetList()



