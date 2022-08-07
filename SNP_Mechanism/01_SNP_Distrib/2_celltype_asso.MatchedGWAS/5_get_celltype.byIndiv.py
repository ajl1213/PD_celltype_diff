#!/python2


import os
import numpy
from scipy import stats


DIR=os.getcwd()
inputFile='/home/ajl1213/Projects/PD/SciAdv_Anal/SNP_Mechanism/01_SNP_Distrib/2_celltype_asso.MatchedGWAS/BED/matchedGWAS.peakIntersect.bed'
totalPeakBED='/home/ajl1213/Projects/PD/data/snATAC/04_GetAtacCount/CellLabel/TOTAL.atacPeak.bed'
cellTypeList=['DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']
minCount=100


def peakFrac():
    #
    totalSize=0
    peakSizeDict={}
    input1=open(totalPeakBED,'r')
    all_input1=input1.readlines()
    for line in all_input1:
	each=line.strip().split('\t')
	peakSize=int(each[2])-int(each[1])
	cellTypeList=each[3]
	for cellType in cellTypeList.split(','):
	    if peakSizeDict.has_key(cellType):
		peakSizeDict[cellType]+=peakSize
	    else:
		peakSizeDict[cellType]=peakSize
	    totalSize+=peakSize
    input1.close()

    #
    peakFracDict={}
    for cellType in peakSizeDict:
	frac= float(peakSizeDict[cellType]) / float(totalSize)
	peakFracDict[cellType]=frac

    return peakFracDict


def indivSNP():
    peakFracDict = peakFrac()

    #
    sampleDict={}
    celltypeDict={}
    input1=open(inputFile,'r')
    all_input1=input1.readlines()
    for line in all_input1:
	each=line.strip().split('\t')
	snpID=each[0]+':'+each[1]
	sampleID=each[3]
	celltypes=each[7]

	if sampleDict.has_key(sampleID):
	    tmp=sampleDict[sampleID]
	    if not snpID in tmp:
		tmp.append(snpID)
	    sampleDict[sampleID]=tmp
	else:
	    sampleDict[sampleID]=[snpID]

	celltypeDict[snpID]=celltypes

    input1.close()

    #
    countDict={}
    nTotalDict={}
    for sampleID in sampleDict:
	tmpDict={}
	idx=0
	for snpID in sampleDict[sampleID]:
	    for celltype in celltypeDict[snpID].split(','):

		if tmpDict.has_key(celltype):
		    tmpDict[celltype]+=1
		else:
		    tmpDict[celltype]=1

		idx+=1

	countDict[sampleID]=tmpDict
	nTotalDict[sampleID]=idx

    for sampleID in nTotalDict:
	print sampleID, nTotalDict[sampleID]

    #
    output1=open(DIR+'/matchedGWAS.celltype_asso.sample.txt','w')
    output1.write('CellType\tSampleID\ttotalCount\tnVar\texpected\tfoldEnrich\tBinomPval\n')
    for sampleID in countDict:
	if nTotalDict[sampleID] > minCount:

	    for cellType in cellTypeList:

		if countDict[sampleID].has_key(cellType):
		    nVar=countDict[sampleID][cellType]
		else:
		    nVar=0

		totalCount=nTotalDict[sampleID]

		nExp=float(totalCount) * peakFracDict[cellType]

		fc=nVar/nExp
		#log2Enrich=numpy.log2(nVar/nExp)

		binomP = stats.binom_test(x=nVar, n=float(totalCount), p=peakFracDict[cellType], alternative='greater')
		newLine=[cellType, sampleID, str(totalCount), str(nVar), str(nExp), str(fc), str(binomP)]
		output1.write('\t'.join(newLine)+'\n')

    output1.close()



indivSNP()


