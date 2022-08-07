#!/python2


import os
import numpy
from scipy import stats


DIR=os.getcwd()
inputFile='/home/ajl1213/Projects/PD/SciAdv_Anal/SNP_Mechanism/01_SNP_Distrib/2a_cRE_asso.MatchedGWAS/BED/matchedGWAS.peakIntersect.bed'
totalPeakBED='/home/ajl1213/Projects/PD/data/snATAC/04_GetAtacCount/CellLabel/TOTAL.atacPeak.bed'


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

#
totalCount=0
countDict={}
input1=open(inputFile,'r')
all_input1=input1.readlines()
for line in all_input1:
    each=line.strip().split('\t')
    cellTypeList=each[-1]

    for cellType in cellTypeList.split(','):
	if countDict.has_key(cellType):
	    countDict[cellType]+=1
	else:
	    countDict[cellType]=1
	totalCount+=1
input1.close()

output1=open(DIR+'/matchedGWAS.celltype_asso.txt','w')
output1.write('CellType\tnVar\texpected\tlog2Enrich\tBinomPval\n')
for cellType in countDict:

    nVar=countDict[cellType]
    nExp=float(totalCount) * peakFracDict[cellType]

    log2Enrich=numpy.log2(nVar/nExp)
    binomP = stats.binom_test(x=nVar, n=float(totalCount), p=peakFracDict[cellType], alternative='greater')
    newLine=[cellType, str(nVar), str(nExp), str(log2Enrich), str(binomP)]
    #print sampleID, peakType, nVar, nExp, log2Enrich, binomP
    output1.write('\t'.join(newLine)+'\n')

output1.close()



