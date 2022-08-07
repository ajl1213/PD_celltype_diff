#!/home/ajl1213/anaconda2/bin/python


import os
import gzip
import numpy
import operator



DIR=os.getcwd()
gtf='/home/ajl1213/genome.info/gencode/hg19/gencode.v19.ptn_coding.lincRNA.level_1_2'
inputDIR='/home/ajl1213/Projects/PD/SciAdv_Data/HiC/FitHiC/OUTPUT'
PeakFile='/home/ajl1213/Projects/PD/data/snATAC/04_GetAtacCount/CellLabel/AtacPeak.TOTAL.bed'
sampleList=['MergedNOSN','MergedPDSN']
newID='MergedTotal'
resolution=[5000, '5kb']
qvalCut=[0.01,'1e-2']
#qvalCut=[0.001,'1e-3']
binFile='/home/ajl1213/Projects/PD/data/HiC/CovNorm/Bin/allchr.'+resolution[1]+'.bin'

os.system('mkdir '+DIR+'/SigInter/')
os.system('mkdir '+DIR+'/BedFiles/')



def sortFitHiC():
    for sampleID in sampleList:
	input1=open(inputDIR+'/'+sampleID+'.'+resolution[1]+'.'+qvalCut[1]+'.SigInter','r') 
	output1=open(DIR+'/SigInter/'+sampleID+'.'+resolution[1]+'.'+qvalCut[1]+'.FitHiC.all.txt','w')
	output1.write('frag1\tfrag2\tdist\trawCount\tqval\n')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')

            chrID=each[0]
            pt1=int(each[1])-(resolution[0]/2)
            pt2=int(each[3])-(resolution[0]/2)
	    dist=pt2-pt1
            rawCount=each[4]
            qval=each[6]

            frag1=chrID+'.'+str(pt1)+'.'+str(pt1+resolution[0])
            frag2=chrID+'.'+str(pt2)+'.'+str(pt2+resolution[0])


	    newLine=[frag1, frag2, str(dist), rawCount, qval]
	    output1.write('\t'.join(newLine)+'\n')

        input1.close()
        output1.close()



def mergeFitHiC():
    interDict={}
    countDict={}
    qvalDict={}

    for sampleID in sampleList:
	input1=open(DIR+'/SigInter/'+sampleID+'.'+resolution[1]+'.'+qvalCut[1]+'.FitHiC.all.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    [frag1, frag2, dist, rawCount, qval]=line.strip().split('\t')
	    interID=frag1+'-'+frag2
	    logQval=-numpy.log10(float(qval))

	    if interDict.has_key(interID):
		tmp=interDict[interID]
		tmp.append(sampleID)
		interDict[interID]=tmp
	    else:
		interDict[interID]=[sampleID]

	    if countDict.has_key(interID):
		tmp=countDict[interID]
		tmp.append(float(rawCount))
		countDict[interID]=tmp
	    else:
		countDict[interID]=[float(rawCount)]

	    if qvalDict.has_key(interID):
		tmp=qvalDict[interID]
		tmp.append(logQval)
		qvalDict[interID]=tmp
	    else:
		qvalDict[interID]=[logQval]

	input1.close()

    output1=open(DIR+'/SigInter/'+newID+'.'+resolution[1]+'.'+qvalCut[1]+'.FitHiC.all.txt','w')
    output1.write('frag1\tfrag2\tsampleList\tlogQval\tsumCount\n')
    for interID in interDict:
	frag1=interID.split('-')[0]
	frag2=interID.split('-')[1]
	samples=';'.join(interDict[interID])
	logQval=numpy.mean(qvalDict[interID])
	countSum=sum(countDict[interID])
	newLine=[frag1, frag2, samples, str(logQval),str(countSum)]
	output1.write('\t'.join(newLine)+'\n')
    output1.close()


def getTssBin():     
    input1=open(gtf,'r')
    output1=open(DIR+'/BedFiles/TSS.bed','w')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	ensembleID=each[0]
	geneID=each[1]
	chrID=each[2]
	tss=each[3]
	output1.write(chrID+'\t'+tss+'\t'+str(int(tss)+1)+'\t'+ensembleID+'\n')
    input1.close()
    output1.close()

    runScript='intersectBed -a '+binFile+' -b '+DIR+'/BedFiles/TSS.bed -wa -wb > '+DIR+'/BedFiles/TSS.binIntersect.bed'
    print runScript
    os.system(runScript)

    tssFragDict={}
    input1=open(DIR+'/BedFiles/TSS.binIntersect.bed','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	tssBin=each[0]+'.'+each[1]+'.'+each[2]
	ensembleID=each[6]
	if tssFragDict.has_key(tssBin):
	    tmp=tssFragDict[tssBin]
	    tmp.append(ensembleID)
	    tssFragDict[tssBin]=tmp
	else:
	    tssFragDict[tssBin]=[ensembleID]
    input1.close()

    return tssFragDict


def getPeakBin():     ## get peak bins
    input1=open(PeakFile,'r')
    output1=open(DIR+'/BedFiles/cRE.bed','w')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	chrID=each[0]
	pt1=each[1]
	pt2=each[2]

	creID=chrID+':'+pt1+'-'+pt2
	output1.write(chrID+'\t'+pt1+'\t'+pt2+'\n')
    input1.close()
    output1.close()

    runScript='cat '+binFile+' | cut -f1,2,3 | intersectBed -a stdin -b '+DIR+'/BedFiles/cRE.bed -wa -wb > '+DIR+'/BedFiles/cRE.binIntersect.bed'
    print runScript
    os.system(runScript)

    creFragDict={}
    input1=open(DIR+'/BedFiles/cRE.binIntersect.bed','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	frag1=each[0]+'.'+each[1]+'.'+each[2]
	creID=each[3]+':'+each[4]+'-'+each[5]
	if creFragDict.has_key(frag1):
	    tmp=creFragDict[frag1]
	    tmp.append(creID)
	    creFragDict[frag1]=tmp
	else:
	    creFragDict[frag1]=[creID]
    input1.close()

    return creFragDict



def getProCreInt():
    tssFragDict=getTssBin()
    creFragDict=getPeakBin()

    input1=open(DIR+'/SigInter/'+newID+'.'+resolution[1]+'.'+qvalCut[1]+'.FitHiC.all.txt','r')
    all_input1=input1.readlines()
    output1=open(DIR+'/SigInter/'+newID+'.'+resolution[1]+'.'+qvalCut[1]+'.FitHiC.ProCre.txt','w')
    output1.write('interLabel\tfrag1\tfrag1Info\tfrag2\tfrag2Info\tdist\tsumCount\n')
    for line in all_input1[1:]:
	[frag1, frag2, samples, logQval, sumCount]=line.strip().split('\t')

	pt1=int(frag1.split('.')[1])
	pt2=int(frag2.split('.')[1])
	dist=pt2-pt1

	if (tssFragDict.has_key(frag1) or tssFragDict.has_key(frag2)):
	    if (tssFragDict.has_key(frag1) and tssFragDict.has_key(frag2)):
		interLabel='PP'
		frag1Info=';'.join(tssFragDict[frag1])
		frag2Info=';'.join(tssFragDict[frag2])
		newLine=[interLabel, frag1, frag1Info, frag2, frag2Info, str(dist), sumCount]
		output1.write('\t'.join(newLine)+'\n')

	    else:
		if tssFragDict.has_key(frag1):
		    frag1Info=';'.join(tssFragDict[frag1])
		    if creFragDict.has_key(frag2):
			interLabel='PC'
			frag2Info=';'.join(creFragDict[frag2])
			newLine=[interLabel, frag1, frag1Info, frag2, frag2Info, str(dist), sumCount]
			output1.write('\t'.join(newLine)+'\n')

		    else:
			interLabel='PO'
			frag2Info='none'
			newLine=[interLabel, frag1, frag1Info, frag2, frag2Info, str(dist), sumCount]
			output1.write('\t'.join(newLine)+'\n')

		if tssFragDict.has_key(frag2):
		    frag2Info=';'.join(tssFragDict[frag2])
		    if creFragDict.has_key(frag1):
			interLabel='PC'
			frag1Info=';'.join(creFragDict[frag1])
			newLine=[interLabel, frag2, frag2Info, frag1, frag1Info, str(dist), sumCount]
			output1.write('\t'.join(newLine)+'\n')

		    else:
			interLabel='PO'
			frag1Info='none'
			newLine=[interLabel, frag2, frag2Info, frag1, frag1Info, str(dist), sumCount]
			output1.write('\t'.join(newLine)+'\n')

	elif (creFragDict.has_key(frag1) or creFragDict.has_key(frag2)):
	    if (creFragDict.has_key(frag1) and creFragDict.has_key(frag2)):
		interLabel='CC'
		frag1Info=';'.join(creFragDict[frag1])
		frag2Info=';'.join(creFragDict[frag2])
		newLine=[interLabel, frag1, frag1Info, frag2, frag2Info, str(dist), sumCount]
		output1.write('\t'.join(newLine)+'\n')

	    else:
		if creFragDict.has_key(frag1):
		    interLabel='CO'
		    frag1Info=';'.join(creFragDict[frag1])
		    frag2Info='none'
		    newLine=[interLabel, frag1, frag1Info, frag2, frag2Info, str(dist), sumCount]
		    output1.write('\t'.join(newLine)+'\n')

		if creFragDict.has_key(frag2):
		    interLabel='CO'
		    frag2Info=';'.join(creFragDict[frag2])
		    frag1Info='none'
		    newLine=[interLabel, frag2, frag2Info, frag1, frag1Info, str(dist), sumCount]
		    output1.write('\t'.join(newLine)+'\n')

	else:
	    interLabel='OO'
	    frag1Info='none'
	    frag2Info='none'
	    newLine=[interLabel, frag1, frag1Info, frag2, frag2Info, str(dist), sumCount]
	    output1.write('\t'.join(newLine)+'\n')

    input1.close()
    output1.close()


def countProInter():

    input1=open(DIR+'/SigInter/'+newID+'.'+resolution[1]+'.'+qvalCut[1]+'.FitHiC.ProCre.txt','r')
    all_input1=input1.readlines()
    targetGeneList=[]
    interCount=0
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	interLabel=each[0]
	if (interLabel=='PC' or interLabel=='PO'):
	    geneList=each[2]
	    interCount+=1
	    for ensembleID in geneList.split(';'):
		targetGeneList.append(ensembleID)
    input1.close()

    print len(set(targetGeneList))
    print str(interCount)




#sortFitHiC()
#mergeFitHiC()
#getProCreInt()
countProInter()



