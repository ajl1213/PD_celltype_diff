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
gtf=sys.argv[6]

windowDict={
'100kb':100000
}


def getDegNumber():
    # bulkDEG
    nTotalBulkDEG=0
    input1=open(DIR+'/DegList/'+dataType+'.'+cellType+'.'+winDist+'.bulkDEG.bed','r')
    all_input1=input1.readlines()
    for line in all_input1:
	nTotalBulkDEG+=1
    input1.close()

    # snRNA DEG
    nTotalSnDEG=0
    input1=open(DIR+'/DegList/'+dataType+'.'+cellType+'.'+winDist+'.snDEG.bed','r')
    all_input1=input1.readlines()
    for line in all_input1:
	nTotalSnDEG+=1
    input1.close()

    return nTotalBulkDEG, nTotalSnDEG


def genRandDEG():
    nTotalBulkDEG, nTotalSnDEG = getDegNumber()    

    # Load TSS
    tssDict={}
    totalGeneList=[]
    input1=open(gtf,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        geneID=each[1]
        chrID=each[2]
        tss=int(each[3])
        geneLoc=chrID+':'+str(tss)
        tssDict[geneID]=geneLoc
	totalGeneList.append(geneID)
    input1.close()

    #
    output1=open(DIR+'/RandDegPermute/'+dataType+'.'+cellType+'.'+winDist+'.RandBulkDEG.bed','w')
    for geneID in random.sample(totalGeneList, nTotalBulkDEG):
	tss=tssDict[geneID]
	chrID=tss.split(':')[0]
	pt=int(tss.split(':')[1])

	pt1=pt-windowDict[winDist]
	if pt1 < 0:
	    pt1=0
	pt2=pt+windowDict[winDist]

	newLine=[chrID, str(pt1), str(pt2), geneID]
	output1.write('\t'.join(newLine)+'\n')
    output1.close()

    #
    output1=open(DIR+'/RandDegPermute/'+dataType+'.'+cellType+'.'+winDist+'.RandSnDEG.bed','w')
    for geneID in random.sample(totalGeneList, nTotalSnDEG):
	tss=tssDict[geneID]
	chrID=tss.split(':')[0]
	pt=int(tss.split(':')[1])

	pt1=pt-windowDict[winDist]
	if pt1 < 0:
	    pt1=0
	pt2=pt+windowDict[winDist]

	newLine=[chrID, str(pt1), str(pt2), geneID]
	output1.write('\t'.join(newLine)+'\n')
    output1.close()


def permuteRandDEGs():
    nTotalBulkDEG, nTotalSnDEG = getDegNumber()    

    output1=open(DIR+'/RandDegPermute/'+dataType+'.'+cellType+'.'+winDist+'.RandDegPermuteRecord.txt','w')
    output1.write('nPermute\tDEG_type\tnDEG\tfraction\n')
    for nIter in range(0, maxIter):
	print str(nIter+1)+'/'+str(maxIter)

	genRandDEG()

	## Bulk DEG
	runScript='intersectBed -a '+DIR+'/RandDegPermute/'+dataType+'.'+cellType+'.'+winDist+'.RandBulkDEG.bed -b '+peakDIR+'/'+dataType+'.'+cellType+'.bed -wa -wb | cut -f4 | sort -u | wc -l > '+DIR+'/RandDegPermute/'+dataType+'.'+cellType+'.'+winDist+'.bulkDEG.RandDegIntersect.bed'
	print runScript
	os.system(runScript)
	#
	input1=open(DIR+'/RandDegPermute/'+dataType+'.'+cellType+'.'+winDist+'.bulkDEG.RandDegIntersect.bed', 'r')
	all_input1=input1.readlines()
	nBulkDEG=all_input1[0].strip().split('\t')[0]
	nBulkDegFrac=float(nBulkDEG)/float(nTotalBulkDEG)*100 if not float(nTotalBulkDEG)==0 else 0
	input1.close()

	newLine=[str(nIter+1), 'bulkDEG', nBulkDEG, str(nBulkDegFrac)]
	output1.write('\t'.join(newLine)+'\n')

	## snRNA DEG
	runScript='intersectBed -a '+DIR+'/RandDegPermute/'+dataType+'.'+cellType+'.'+winDist+'.RandSnDEG.bed -b '+peakDIR+'/'+dataType+'.'+cellType+'.bed -wa -wb | cut -f4 | sort -u | wc -l > '+DIR+'/RandDegPermute/'+dataType+'.'+cellType+'.'+winDist+'.snDEG.RandDegIntersect.bed'
	print runScript
	os.system(runScript)
	#
	input1=open(DIR+'/RandDegPermute/'+dataType+'.'+cellType+'.'+winDist+'.snDEG.RandDegIntersect.bed', 'r')
	all_input1=input1.readlines()
	nSnDEG=all_input1[0].strip().split('\t')[0]
	nSnDegFrac=float(nSnDEG)/float(nTotalSnDEG)*100 if not float(nTotalSnDEG)==0 else 0
	input1.close()

	newLine=[str(nIter+1), 'snDEG', nSnDEG, str(nSnDegFrac)]
	output1.write('\t'.join(newLine)+'\n')

    output1.close()




permuteRandDEGs()



