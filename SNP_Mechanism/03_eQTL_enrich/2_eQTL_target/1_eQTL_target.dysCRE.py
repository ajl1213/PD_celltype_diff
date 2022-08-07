#!/home/ajl1213/anaconda2/bin/python


import os
import numpy



DIR=os.getcwd()
GTF='/home/ajl1213/genome.info/gencode/hg19/gencode.v19.ptn_coding.lincRNA.level_1_2'
eQTL_DIR='/home/ajl1213/Projects/PD/SciAdv_Anal/SNP_Mechanism/03_eQTL_enrich/1_eQTL_enrich/BED'
InterTarget_DIR='/home/ajl1213/Projects/PD/SciAdv_Anal/CombinedTargetID/01_FindTarget/RESULT'
DysPeak_DIR='/home/ajl1213/Projects/PD/NeuronJ/DiffAnal/02_H3K27ac/1_BulkH3K27ac/BedFiles'

dataTypeList=['Down','Up','MergedDysPeak']
cellTypeList=['DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']
winSize=[1000000,'1mb']
minAbcScore=0

os.system('mkdir '+DIR+'/BED/')


def loadGeneInfo():
    geneDict={}
    input1=open(GTF,'r')
    all_input1=input1.readlines()
    output1=open(DIR+'/BED/GeneTSS.'+winSize[1]+'.bed','w')
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        ensembleID=each[0]
        geneID=each[1]
        chrID=each[2]
        tssPos=each[3]
        if each[6]=='protein_coding':
            geneDict[ensembleID]=geneID

	    if (int(tssPos)-winSize[0]) < 0:
		pt1=0
	    else:
		pt1=int(tssPos)-winSize[0]

	    pt2=int(tssPos)+winSize[0]

            output1.write(chrID+'\t'+str(pt1)+'\t'+str(pt2)+'\t'+ensembleID+'\t'+geneID+'\t'+tssPos+'\n')
    input1.close()
    output1.close()

    return geneDict


def load_eQTL():
    geneDict=loadGeneInfo()

    eqtlDict={}
    for dataType in dataTypeList:
	for cellType in cellTypeList:
	    tmpDict={}
	    input1=open(eQTL_DIR+'/'+dataType+'.'+cellType+'.peakIntersect.bed','r')
	    all_input1=input1.readlines()
	    for line in all_input1:
		each=line.strip().split('\t')
		eqtlID=each[0]+';'+each[1]+';'+each[3]+';'+each[4]
		ensembleID=each[5]
		peakID=each[10]
		assoID=ensembleID+';'+peakID
		if geneDict.has_key(ensembleID):
		    if tmpDict.has_key(assoID):
			tmp=tmpDict[assoID]
			if not eqtlID in tmp:
			    tmp.append(eqtlID)
			tmpDict[assoID]=tmp
		    else:
			tmpDict[assoID]=[eqtlID]
	    input1.close()

	    eqtlDict[dataType+';'+cellType]=tmpDict

    return eqtlDict


def load_HiC():
    interDict={}
    for i in dataTypeList:
	for cellType in cellTypeList:
	    tmpDict={}
	    input1=open(InterTarget_DIR+'/'+cellType+'.fullList.txt','r')
	    all_input1=input1.readlines()
	    for line in all_input1[1:]:
		each=line.strip().split('\t')

		ensembleID=each[0]
		geneID=each[1]
		peakID=each[5]
		cellTypes=each[6]
		dataType=each[7].split('Peak')[0]
		abcScore=each[11]
		if float(abcScore) > minAbcScore:
		    assoID=ensembleID+';'+peakID
		    if dataType==i:
			tmpDict[assoID]=abcScore
	    input1.close()

	    interDict[i+';'+cellType]=tmpDict

	if i =='MergedDysPeak':
	    for cellType in cellTypeList:
		tmpDict={}
		input1=open(InterTarget_DIR+'/'+cellType+'.fullList.txt','r')
		all_input1=input1.readlines()
		for line in all_input1[1:]:
		    each=line.strip().split('\t')

		    ensembleID=each[0]
		    geneID=each[1]
		    peakID=each[5]
		    cellTypes=each[6]
		    dataType=each[7]
		    if dataType.count('Peak') > 0:

			abcScore=each[11]
			if float(abcScore) > minAbcScore:
			    assoID=ensembleID+';'+peakID
			    tmpDict[assoID]=abcScore
		input1.close()

		interDict[i+';'+cellType]=tmpDict

    return interDict


def getAllAsso():
    distDict={}
    for dataType in dataTypeList:
	for cellType in cellTypeList:
	    runLine='intersectBed -a '+DIR+'/BED/GeneTSS.'+winSize[1]+'.bed -b '+DysPeak_DIR+'/'+dataType+'.'+cellType+'.bed -wa -wb > '+DIR+'/BED/'+dataType+'.'+cellType+'.GeneTSS.Peak.intersect.bed'
	    print runLine
	    os.system(runLine)

	    tmpDict={}
	    input1=open(DIR+'/BED/'+dataType+'.'+cellType+'.GeneTSS.Peak.intersect.bed', 'r')
	    all_input1=input1.readlines()
	    for line in all_input1:
		each=line.strip().split('\t')
		ensembleID=each[3]
		geneID=each[4]
		tssPos=each[5]
		peakID=each[6]+':'+each[7]+'-'+each[8]

		assoID=ensembleID+';'+peakID
		dist=min(abs(int(tssPos)-int(each[7])), abs(int(tssPos)-int(each[8])))

		tmpDict[assoID]=dist
	    input1.close()

	    output1=open(DIR+'/BED/'+dataType+'.'+cellType+'.AllAsso.txt','w')
	    for assoID in tmpDict:
		output1.write(assoID+'\n')
	    output1.close()

	    distDict[dataType+';'+cellType]=tmpDict

    return distDict


def getMatch():
    geneDict=loadGeneInfo()
    eqtlDict = load_eQTL()
    interDict = load_HiC()
    distDict=getAllAsso()

    output1=open(DIR+'/eQTL_HiC.comparison.list.txt','w')
    output1.write('celltype\tdatatype\tpeakID\teQTL\tensembleID\tgeneID\tassoLabel\tabcScore\tdist\n')
    for dataType in dataTypeList:
	for cellType in cellTypeList:

	    unionSet=set(eqtlDict[dataType+';'+cellType]) | set(interDict[dataType+';'+cellType])
	    for assoID in unionSet:
		ensembleID=assoID.split(';')[0]
		peakID=assoID.split(';')[1]
		if distDict[dataType+';'+cellType].has_key(assoID):
		    dist=str(distDict[dataType+';'+cellType][assoID])

		    if (eqtlDict[dataType+';'+cellType].has_key(assoID) and interDict[dataType+';'+cellType].has_key(assoID)):
			assoLabel='common'

			for eqtlID in eqtlDict[dataType+';'+cellType][assoID]: 
			    newLine=[cellType, dataType, peakID, eqtlID, ensembleID, geneDict[ensembleID], assoLabel, interDict[dataType+';'+cellType][assoID], dist]
			    output1.write('\t'.join(newLine)+'\n')

		    else:
			if eqtlDict[dataType+';'+cellType].has_key(assoID):
			    assoLabel='eqtl_spec'

			    for eqtlID in eqtlDict[dataType+';'+cellType][assoID]: 
				newLine=[cellType, dataType, peakID, eqtlID, ensembleID, geneDict[ensembleID], assoLabel, '.', dist]
				output1.write('\t'.join(newLine)+'\n')

			else:
			    assoLabel='hic_spec'

			    newLine=[cellType, dataType, peakID, '.', ensembleID, geneDict[ensembleID], assoLabel, interDict[dataType+';'+cellType][assoID], dist]
			    output1.write('\t'.join(newLine)+'\n')

    output1.close()



def makeData():

    output1=open(DIR+'/eQTL_HiC.comparison.count.txt','w')
    output1.write('Datatype\tCelltype\teqtl_spec\thic_spec\tcommon\tallAsso\n')
    for dataType in dataTypeList:
	for cellType in cellTypeList:
	    #
	    allcount=0
	    input1=open(DIR+'/BED/'+dataType+'.'+cellType+'.AllAsso.txt','r')
	    all_input1=input1.readlines()
	    for line in all_input1:
		allcount+=1
	    input1.close() 

	    #
	    countDict={}
	    input1=open(DIR+'/eQTL_HiC.comparison.list.txt','r')
	    all_input1=input1.readlines()
	    for line in all_input1[1:]:
		each=line.strip().split('\t')

		assoLabel=each[6]
		if dataType+';'+cellType==each[1]+';'+each[0]:

		    if countDict.has_key(assoLabel):
			countDict[assoLabel]+=1
		    else:
			countDict[assoLabel]=1
	    input1.close()

	    # 
	    newLine=[dataType, cellType, str(countDict['eqtl_spec']), str(countDict['hic_spec']), str(countDict['common']), str(allcount)]
	    output1.write('\t'.join(newLine)+'\n')
    output1.close()



getMatch()
makeData()


