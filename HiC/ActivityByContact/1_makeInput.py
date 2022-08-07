#!/home/ajl1213/anaconda2/bin/python


import os
import numpy



DIR=os.getcwd()
resolution=['5kb',5000]
winSize=['1mb', 1000000]
gtf='/home/ajl1213/genome.info/gencode/hg19/gencode.v19.ptn_coding.lincRNA.level_1_2'
peakDIR='/home/ajl1213/Projects/PD/data/snATAC/04_GetAtacCount/CellLabel'
cellTypeList=['DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']


os.system('mkdir '+DIR+'/BedFiles')
os.system('mkdir '+DIR+'/CreListPerGene')


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


def getTssBin():    
    input1=open(gtf,'r')
    output1=open(DIR+'/BedFiles/TSS.'+winSize[0]+'.bed','w')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        ensembleID=each[0]
        chrID=each[2]
        tss=int(each[3])
	tssBin=chrID+'.'+str(int(tss/resolution[1])*resolution[1])+'.'+str(int(tss/resolution[1]+1)*resolution[1])
	pt1=tss-winSize[1] if not tss-winSize[1] < 0 else 0
	pt2=tss+winSize[1]
        output1.write(chrID+'\t'+str(pt1)+'\t'+str(pt2)+'\t'+ensembleID+'\t'+tssBin+'\n')
    input1.close()
    output1.close()


def getCreBin():     ### get cRE bins
    for cellType in cellTypeList:
	input1=open(peakDIR+'/'+cellType+'.atacPeak.bed','r')
	output1=open(DIR+'/BedFiles/'+cellType+'.cRE.bed','w')
	all_input1=input1.readlines()
	for line in all_input1:
	    each=line.strip().split('\t')
	    chrID=each[0]
	    pt1=int(each[1])
	    pt2=int(each[2])
	    creID=chrID+':'+str(pt1)+'-'+str(pt2)
	    crePos=(pt1+pt2)/2
	    creBin=chrID+'.'+str(int(crePos/resolution[1])*resolution[1])+'.'+str(int(crePos/resolution[1]+1)*resolution[1])
	    newLine=[chrID, str(crePos), str(crePos+1), creID, creBin]
	    output1.write('\t'.join(newLine)+'\n')
	input1.close()
	output1.close()


def getCreListPerGene():
    geneDict=loadGeneInfo()
    getTssBin()
    getCreBin()

    for cellType in cellTypeList:
	runScript='intersectBed -a '+DIR+'/BedFiles/TSS.'+winSize[0]+'.bed -b '+DIR+'/BedFiles/'+cellType+'.cRE.bed -wa -wb > '+DIR+'/BedFiles/TSS.cRE.'+cellType+'.'+winSize[0]+'.Intersect.bed'
	print runScript
	os.system(runScript)

	allDict={}
	input1=open(DIR+'/BedFiles/TSS.cRE.'+cellType+'.'+winSize[0]+'.Intersect.bed','r')
	all_input1=input1.readlines()
	for line in all_input1:
	    each=line.strip().split('\t')
	    ensembleID=each[3]
	    tssBin=each[4]
	    creID=each[8]
	    creBin=each[9]

	    keyID=ensembleID+';'+tssBin
	    valID=creID+';'+creBin

	    if allDict.has_key(keyID):
		tmp=allDict[keyID]
		tmp.append(valID)
		allDict[keyID]=tmp
	    else:
		allDict[keyID]=[valID]
	input1.close()

	output1=open(DIR+'/CreListPerGene/'+cellType+'.CreListPerGene.txt','w')
	output1.write('ensembleID\tgeneID\ttssBin\tcreID\tcreBin\tdist\n')
	for keyID in allDict:
	    [ensembleID, tssBin]=keyID.split(';')
	    for valID in allDict[keyID]:
		[creID, creBin]=valID.split(';')
		tssBin1=int(tssBin.split('.')[1])
		creBin1=int(creBin.split('.')[1])
		dist=str(abs(tssBin1-creBin1))

		geneID=geneDict[ensembleID]
		newLine=[ensembleID, geneID, tssBin, creID, creBin, dist]
		output1.write('\t'.join(newLine)+'\n')
	output1.close()



getCreListPerGene()



