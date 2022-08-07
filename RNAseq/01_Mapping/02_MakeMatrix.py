#!/home/ajl1213/anaconda2/bin/python


import os
import numpy


DIR=os.getcwd()
gtf='/home/ajl1213/genome.info/gencode/hg19.release38/GTF/gencode.v38.ptncoding.lncRNA.level_1_2.txt'
sampleList=[
'X4870NOSN',
'X4689NOSN', 
'X5628NOSN', 
'X5130_1NOSN',
'X5130_2NOSN',
'X4996NOSN', 
'X5006_1NOSN',
'X5006_2NOSN',
'X5244PDSN',
'X4831PDSN',
'X5215PDSN',
'X5627PDSN',
'X5778PDSN',
'X5742PDSN',
'X5649PDSN',
'X5591PDSN'
]

print len(sampleList)

os.system('mkdir '+DIR+'/CountMatrix')


def getGeneSet():
    GeneDict={}
    input1=open(gtf,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	ensembleID=each[0]
	geneID=each[1]

	GeneDict[ensembleID]=geneID
    input1.close()

    return GeneDict


def makeSubjectInfo():
    infoDict={}
    input1=open(DIR+'/Metadata/metadata.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        sampleID=each[0]
        info=each[1:]
        infoDict[sampleID]=info
    input1.close()

    output1=open(DIR+'/Metadata/metadata.final.txt','w')
    output1.write('SampleID\tIndivID\tDiagnosis\tAGE\tGender\tBraakCombo\tMMSE\tDRS\n')
    for i in sampleList:
	if i.count('_')>0:
	    if i.count('NOSN')>0:
		tail='NOSN'
	    if i.count('PDSN')>0:
		tail='PDSN'

	    sampleID=i.split('_')[0]+tail

	else:
	    sampleID=i

        sampleID=sampleID.split(';')[0]
        output1.write(i+'\t'+sampleID+'\t')
        output1.write('\t'.join(infoDict[sampleID])+'\n')
    output1.close()


def getStat():
    output1=open(DIR+'/AlignedBam/MapStat.txt','w')
    output1.write('sampleID\tReadTotal\tUniqMap(%)\tUniqMapRead\n')
    for sampleID in sampleList:
	input1=open(DIR+'/AlignedBam/'+sampleID+'/Log.final.out','r')
	all_input1=input1.readlines()
	for line in all_input1:
	    each=line.strip().split('\t')
	    if line.count('Number of input reads')>0:
		inputReadNum=each[1]
	    if line.count('Uniquely mapped reads %')>0:
		mapPercent=each[1]
	    if line.count('Uniquely mapped reads number')>0:
		uniqReadNum=each[1]
	output1.write(sampleID+'\t'+inputReadNum+'\t'+mapPercent+'\t'+uniqReadNum+'\n')
	input1.close()
    output1.close()



def getCountInfo():
    GeneDict=getGeneSet()

    totalCountDict={}
    ValDict={}
    for sampleID in sampleList:
	tmpDict={}
	tmpCount=0
	input1=open(DIR+'/CountInfo/'+sampleID+'.genes.results','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    ensembleID=each[0].split('_')[0]
	    transcriptID=each[1]
	    count_val=each[4]
	    tpm_val=each[5]
	    fpkm_val=each[6]

	    tmpCount+=float(count_val)
	    tmpDict[ensembleID]=count_val
	input1.close()

	totalCountDict[sampleID]=tmpCount
	ValDict[sampleID]=tmpDict

    return totalCountDict, ValDict


def makeMatrix():
    GeneDict=getGeneSet()
    totalCountDict, ValDict=getCountInfo()

    output1=open(DIR+'/CountMatrix/RNA.rawCount.txt','w')
    output1.write('ensembleID\tgeneID\t'+'\t'.join(sampleList)+'\n')
    for ensembleID in GeneDict:
	geneID=GeneDict[ensembleID]
	if ValDict[sampleList[0]].has_key(ensembleID):
	    output1.write(ensembleID+'\t'+geneID)
	    for sampleID in sampleList:
		count_val=ValDict[sampleID][ensembleID]
		output1.write('\t'+count_val)
	    output1.write('\n')
    output1.close()





makeSubjectInfo()
getStat()
makeMatrix()


