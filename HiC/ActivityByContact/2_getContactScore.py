#!/home/ajl1213/anaconda2/bin/python


import os
import numpy
import gzip



DIR=os.getcwd()
gtf='/home/ajl1213/genome.info/gencode/hg19/gencode.v19.ptn_coding.lincRNA.level_1_2'
InterCountDIR='/home/ajl1213/Projects/PD/SciAdv_Data/HiC/CovNorm/FeatureVec'
PeakCountFile='/home/ajl1213/Projects/PD/data/snATAC/04_GetAtacCount/PeakCount/atacCount.byCellType.qq.txt'
chrList=['chr1', 'chr2', 'chr3','chr4', 'chr5', 'chr6','chr7', 'chr8', 'chr9','chr10', 'chr11', 'chr12','chr13', 'chr14', 'chr15','chr16', 'chr17', 'chr18','chr19', 'chr20', 'chr21','chr22', 'chrX']
inputDIR=DIR+'/CreListPerGene'
resolution=['5kb',5000]
sampleList=['MergedNOSN','MergedPDSN']
cellTypeList=['DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']

os.system('mkdir '+DIR+'/ActivityByContactScore')



def loadGeneInfo():
    geneDict={}
    input1=open(gtf,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        ensembleID=each[0]
        geneID=each[1]
        geneDict[ensembleID]=geneID
    input1.close()

    return geneDict


def loadInterCount():
    interCountDict={}
    print 'loading HiC count..'
    for chrID in chrList:
	print chrID
	tmpDict={}
	for sampleID in sampleList:
	    input1=gzip.open(InterCountDIR+'/'+sampleID+'.'+chrID+'.'+resolution[0]+'.count.gz','r')
	    all_input1=input1.readlines()
	    for line in all_input1[1:]:
		each=line.strip().split('\t')
		frag1=each[0]
		frag2=each[1]
		freqVal=float(each[4])

		interID1=frag1+'-'+frag2
		interID2=frag2+'-'+frag1

		##
		if tmpDict.has_key(interID1):
		    tmp=tmpDict[interID1]
		    tmp.append(freqVal)
		    tmpDict[interID1]=tmp
		else:
		    tmpDict[interID1]=[freqVal]
		##
		if tmpDict.has_key(interID2):
		    tmp=tmpDict[interID2]
		    tmp.append(freqVal)
		    tmpDict[interID2]=tmp
		else:
		    tmpDict[interID2]=[freqVal]
		input1.close()

	for interID in tmpDict:
	    sumInter=sum(tmpDict[interID])
	    interCountDict[interID]=sumInter

    print 'done'

    return interCountDict


def loadPeakCount():
    print 'loading Peak count..'
    peakCountDict={}
    input1=open(PeakCountFile,'r')
    all_input1=input1.readlines()
    cellTypeList=all_input1[0].strip().split('\t')[1:]
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        peakID=each[0]
        vals=numpy.array(each[1:], dtype='float')
        tmpDict={}
        index=0
        for val in vals:
            tmpDict[cellTypeList[index]]=val
            index+=1
        peakCountDict[peakID]=tmpDict
    input1.close()

    print 'done'

    return peakCountDict


def computeContactScore():
    geneDict=loadGeneInfo()
    interCountDict=loadInterCount()
    peakCountDict=loadPeakCount()

    print 'computing ABC score..'
    for cellType in cellTypeList:
	print cellType

	totalScoreDict={}
	input1=open(inputDIR+'/'+cellType+'.CreListPerGene.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    ensembleID=each[0]
	    tssBin=each[2]
	    creID=each[3]
	    creBin=each[4]

	    peakCount=peakCountDict[creID][cellType]

	    interID=tssBin+'-'+creBin
	    dist=int(each[5])

	    if interCountDict.has_key(interID):
		if dist > 0:
		    hicCount=float(interCountDict[interID])/2
		else:
		    frag1=tssBin.split('.')[0]+'.'+str(int(tssBin.split('.')[1])-resolution[1])+'.'+str(int(tssBin.split('.')[2])-resolution[1])
		    frag2=tssBin.split('.')[0]+'.'+str(int(tssBin.split('.')[1])+resolution[1])+'.'+str(int(tssBin.split('.')[2])+resolution[1])
		    interID1=frag1+'-'+tssBin
		    interID2=frag2+'-'+tssBin

		    count1=float(interCountDict[interID1])/2
		    count2=float(interCountDict[interID2])/2

		    hicCount=float(count1+count2)/2  ###

		contactVal=hicCount*peakCount
		if totalScoreDict.has_key(ensembleID):
		    totalScoreDict[ensembleID]+=contactVal
		else:
		    totalScoreDict[ensembleID]=contactVal

	    #else:
	#	print 'no Hi-C count : '+ensembleID+' '+tssBin+' '+creID+' '+creBin+' '+str(peakCount)+' '+str(dist)

	input1.close()

	##
	output1=open(DIR+'/ActivityByContactScore/'+cellType+'.SumContactPerGene.txt','w')
	output1.write('ensembleID\tgeneID\tsumContactVal\n')
	for ensembleID in totalScoreDict:
	    geneID=geneDict[ensembleID]
	    sumContactVal=str(totalScoreDict[ensembleID])
	    newLine=[ensembleID, geneID, sumContactVal]
	    output1.write('\t'.join(newLine)+'\n')
	output1.close()

	##
	output1=open(DIR+'/ActivityByContactScore/'+cellType+'.ActivityByContactScore.txt','w')
	output1.write('ensembleID\tgeneID\ttssBin\tcreID\tcreBin\tdist\tpeakCount\tinterCount\tcontactVal\tActivityByContactScore\n')
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    ensembleID=each[0]
	    geneID=each[1]
	    tssBin=each[2]
	    creID=each[3]
	    creBin=each[4]

	    peakCount=peakCountDict[creID][cellType]

	    interID=tssBin+'-'+creBin
	    dist=int(each[5])

	    if interCountDict.has_key(interID):
		if dist > 0:
		    hicCount=float(interCountDict[interID])/2
		else:
		    frag1=tssBin.split('.')[0]+'.'+str(int(tssBin.split('.')[1])-resolution[1])+'.'+str(int(tssBin.split('.')[2])-resolution[1])
		    frag2=tssBin.split('.')[0]+'.'+str(int(tssBin.split('.')[1])+resolution[1])+'.'+str(int(tssBin.split('.')[2])+resolution[1])
		    interID1=frag1+'-'+tssBin
		    interID2=frag2+'-'+tssBin

		    count1=float(interCountDict[interID1])/2
		    count2=float(interCountDict[interID2])/2

		    hicCount=float(count1+count2)/2  ##

		contactVal=hicCount*peakCount
		sumScore=totalScoreDict[ensembleID]
		activityByContactScore= contactVal/sumScore*100 if not sumScore==0 else 0

		newLine=[ensembleID, geneID, tssBin, creID, creBin, str(dist), str(peakCount), str(hicCount), str(contactVal), str(activityByContactScore)]
		output1.write('\t'.join(newLine)+'\n')

	output1.close()



computeContactScore()



