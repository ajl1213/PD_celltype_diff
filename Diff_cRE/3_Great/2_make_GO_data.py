#!/home/ajl1213/anaconda2/bin/python


import os
import numpy
import operator


DIR=os.getcwd()
inputDIR=DIR+'/GOBP'
#inputDIR=DIR+'/GOMF'
#inputDIR=DIR+'/MousePheno'

dataTypeList=['Down','Up']
cellTypeList=['DopaN','GabaN','Oligo','OPC','Ast','Micro','Endo','Peri']
nTerm=5


for dataType in dataTypeList:
    ## Load data
    allDict={}
    pvalDict={}
    for cellType in cellTypeList:
	tmp1Dict={}
	tmp2Dict={}
	input1=open(inputDIR+'/'+dataType+'.'+cellType+'.txt', 'r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    termID=each[2].replace('\'','')
	    termID='_'.join(termID.split(' '))
	    pval=each[3]
	    regionFC=each[4]
	    geneFC=each[5]
	    regions=each[8]
	    genes=each[9]

	    tmp1Dict[termID]=[pval, regionFC, geneFC, regions, genes]
	    tmp2Dict[termID]=-numpy.log10(float(pval))

	input1.close()

	allDict[cellType]=tmp1Dict
	pvalDict[cellType]=tmp2Dict

    ## Select GO terms
    termList=[]
    for cellType in cellTypeList:
	idx=0
	for key in sorted(pvalDict[cellType].items(), key=operator.itemgetter(1), reverse=True):
	    termID=key[0]
	    logPval=key[1]

	    if idx < nTerm:
		if not termID in termList:
		    print cellType, termID, str(logPval)
		    termList.append(termID)
		    idx+=1
		
    ## Make data
    output1=open(DIR+'/'+dataType+'.plotData.txt','w')
    output1.write('termID\tcellType\tlogPval\tregionFC\tgeneFC\tregions\tgenes\n')
    for termID in termList:

	for cellType in cellTypeList:

	    if allDict[cellType].has_key(termID):
		[pval, regionFC, geneFC, regions, genes]=allDict[cellType][termID]
		logPval=-numpy.log10(float(pval))
		newLine=[termID, cellType, str(logPval), regionFC, geneFC, regions, genes]
		output1.write('\t'.join(newLine)+'\n')
	    else:
		newLine=[termID, cellType, '0.0', '0.0', '0.0', 'none', 'none']
		output1.write('\t'.join(newLine)+'\n')
    output1.close()

    ## plot matrix
    output1=open(DIR+'/'+dataType+'.plotMat.txt','w')
    output1.write('termID\t'+'\t'.join(cellTypeList)+'\n')
    for termID in termList:

	valList=[]
	for cellType in cellTypeList:

	    if allDict[cellType].has_key(termID):
		[pval, regionFC, geneFC, regions, genes]=allDict[cellType][termID]
		logPval=-numpy.log10(float(pval))
		valList.append(str(logPval))

	    else:
		valList.append('0.0')
	output1.write(termID+'\t'+'\t'.join(valList)+'\n')

    output1.close()






