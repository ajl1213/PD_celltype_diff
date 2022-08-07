#!/home/ajl1213/anaconda2/bin/python


import os
import numpy
import operator


DIR=os.getcwd()
inputDIR=DIR+'/CellTypeDEGs/GO'
dataTypeList=['Down','Up']
cellTypeList=['DopaN','Oligo','OPC','Ast','Micro','Endo','Peri']
nTerm=5



for dataType in dataTypeList:
    print dataType

    input1=open(inputDIR+'/'+dataType+'DEGs.GO.txt', 'r')
    all_input1=input1.readlines()

    ## Load data
    allDict={}
    pvalDict={}
    for cellType in cellTypeList:
	tmp1Dict={}
	tmp2Dict={}
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    termID=each[0].replace('\'','')
	    termID='_'.join(termID.split(' ')[:-1])
	    pval=each[2]
	    adjPval=each[3]
	    oddRatio=each[6]
	    combinedScore=each[7]
	    genes=each[8]
	    cellLabel=each[9]
	    db=each[11]

	    if cellType==cellLabel:
		tmp1Dict[termID]=[pval, oddRatio, combinedScore, genes]
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
    output1=open(inputDIR+'/'+dataType+'DEGs.plotData.txt','w')
    output1.write('termID\tcellType\tlogPval\toddRatio\tcombinedScore\tgenes\n')
    for termID in termList:

	for cellType in cellTypeList:

	    if allDict[cellType].has_key(termID):
		[pval, oddRatio, combinedScore, genes]=allDict[cellType][termID]
		logPval=-numpy.log10(float(pval))
		newLine=[termID, cellType, str(logPval), oddRatio, combinedScore, genes]
		output1.write('\t'.join(newLine)+'\n')
	    else:
		newLine=[termID, cellType, '0.0', '0.0', '0.0', 'none']
		output1.write('\t'.join(newLine)+'\n')
    output1.close()




