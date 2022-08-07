#!/home/ajl1213/anaconda2/bin/python


import os
import sys
import numpy
import random



DIR=os.getcwd()
dataType=sys.argv[1]
cellType=sys.argv[2]
winDist=sys.argv[3]
peakDIR=sys.argv[4]
maxIter=int(sys.argv[5])
fai=sys.argv[6]


def getDegNumber():
    # Merged DEGs
    nTotalDEG=0
    input1=open(DIR+'/DegList/'+dataType+'.'+cellType+'.'+winDist+'.MergedDEG.bed','r')
    all_input1=input1.readlines()
    for line in all_input1:
	nTotalDEG+=1
    input1.close()

    return nTotalDEG


def loadSigPeaks():
    PeakDict={}
    input1=open(peakDIR+'/'+dataType+'.'+cellType+'.bed','r')
    all_input1=input1.readlines()
    for line in all_input1:
	each=line.strip().split('\t')
	chrID=each[0]
	pt1=each[1]
	pt2=each[2]
	peakSize=int(pt2)-int(pt1)

	if PeakDict.has_key(chrID):
	    tmp=PeakDict[chrID]
	    tmp.append(peakSize)
	    PeakDict[chrID]=tmp
	else:
	    PeakDict[chrID]=[peakSize]
    input1.close()

    return PeakDict


## Number of peaks matched in each chr
def genRandPeak():
    PeakDict = loadSigPeaks()

    genomeDict={}
    input1=open(fai,'r')
    all_input1=input1.readlines()
    for line in all_input1:
        each=line.strip().split('\t')
        chrID=each[0]
        chrSize=each[1]
	genomeDict[chrID]=chrSize
    input1.close()

    output1=open(DIR+'/RandPeakPermute/'+dataType+'.'+cellType+'.'+winDist+'.RandPeak.bed','w')
    for chrID in PeakDict:
        for peakSize in PeakDict[chrID]:
            pt1=random.randrange(0, int(genomeDict[chrID])-peakSize)
            pt2=pt1+peakSize
            output1.write(chrID+'\t'+str(pt1)+'\t'+str(pt2)+'\n')
    output1.close()


def permuteRandPeaks():
    nTotalDEG = getDegNumber()

    output1=open(DIR+'/RandPeakPermute/'+dataType+'.'+cellType+'.'+winDist+'.RandPeakPermuteRecord.txt','w')
    output1.write('nPermute\tDEG_type\tnDEG\tfraction\n')
    for nIter in range(0, maxIter):
	print str(nIter+1)+'/'+str(maxIter)

	genRandPeak()

	## Merged DEGs
	runScript='intersectBed -a '+DIR+'/RandPeakPermute/'+dataType+'.'+cellType+'.'+winDist+'.RandPeak.bed -b '+DIR+'/DegList/'+dataType+'.'+cellType+'.'+winDist+'.MergedDEG.bed -wa -wb | cut -f7 | sort -u | wc -l > '+DIR+'/RandPeakPermute/'+dataType+'.'+cellType+'.'+winDist+'.MergedDEG.RandPeakIntersect.bed'
	print runScript
	os.system(runScript)
	#
	input1=open(DIR+'/RandPeakPermute/'+dataType+'.'+cellType+'.'+winDist+'.MergedDEG.RandPeakIntersect.bed', 'r')
	all_input1=input1.readlines()
	nDEG=all_input1[0].strip().split('\t')[0]
	DegFrac=float(nDEG)/float(nTotalDEG)*100 if not float(nTotalDEG)==0 else 0
	input1.close()

	newLine=[str(nIter+1), 'MergedDEG', nDEG, str(DegFrac)]
	output1.write('\t'.join(newLine)+'\n')

    output1.close()




permuteRandPeaks()


