#!/home/ajl1213/anaconda2/bin/python


import os
import numpy


DIR=os.getcwd()
GTF='/home/ajl1213/genome.info/gencode/hg19/gencode.v19.ptn_coding.lincRNA.level_1_2'
GWAS_file='/home/ajl1213/Projects/PD/SciAdv_Anal/GWAS_Enrich/1_CollectSNP/LdFinal/MergedPD.allchr.txt'
sigEQTL_file='/home/ajl1213/Projects/PD/SciAdv_Anal/SNP_Mechanism/03_eQTL_enrich/1_eQTL_enrich/BED/BrainSN.sigEQTL.bed'
totalPeak_file='/home/ajl1213/Projects/PD/data/snATAC/04_GetAtacCount/CellLabel/AtacPeak.TOTAL.bed'
InterTarget_DIR='/home/ajl1213/Projects/PD/SciAdv_Anal/CombinedTargetID/01_FindTarget/MergedTargetList'
InterTarget_DIR='/home/ajl1213/Projects/PD/SciAdv_Anal/CombinedTargetID/01_FindTarget/RESULT'
dataTypeList=['Down','Up']
cellTypeList=['DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']
winSize=[1000000,'1mb']
minAbcScore=0

os.system('mkdir '+DIR+'/BED/')


def loadGeneInfo():
    geneDict={}
    tssDict={}
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
	    tssDict[ensembleID]=tssPos

            if (int(tssPos)-winSize[0]) < 0:
                pt1=0
            else:
                pt1=int(tssPos)-winSize[0]

            pt2=int(tssPos)+winSize[0]

            output1.write(chrID+'\t'+str(pt1)+'\t'+str(pt2)+'\t'+ensembleID+'\t'+geneID+'\t'+tssPos+'\n')
    input1.close()
    output1.close()

    return geneDict, tssDict


def loadGWAS(): 
    print 'loading GWAS'
    gwasDict={}
    input1=open(GWAS_file,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        snpPos=each[0].replace('.',':')
        if not gwasDict.has_key(snpPos):
            gwasDict[snpPos]=1
    input1.close()
    print 'loading complete'

    return gwasDict


def loadEQTL():
    geneDict, tssDict=loadGeneInfo()

    print 'loading eQTL'
    eqtlDict={}
    input1=open(sigEQTL_file,'r')
    all_input1=input1.readlines()
    for line in all_input1:
        each=line.strip().split('\t')
        eqtlPos=each[0]+':'+each[1]
	ensembleID=each[5]
	if geneDict.has_key(ensembleID):
	    if eqtlDict.has_key(eqtlPos):
		tmp=eqtlDict[eqtlPos]
		tmp.append(ensembleID)
		eqtlDict[eqtlPos]=tmp
	    else:
		eqtlDict[eqtlPos]=[ensembleID]
    input1.close()
    print 'loading complete'

    return eqtlDict


def match():
    geneDict, tssDict=loadGeneInfo()
    gwasDict=loadGWAS()
    eqtlDict=loadEQTL()

    commonList=set(gwasDict) & set(eqtlDict)
    output1=open(DIR+'/BED/GWAS_eQTL.common.bed','w')
    for i in commonList:
	chrID=i.split(':')[0]
	pt=int(i.split(':')[1])
	pt1=str(pt)
	pt2=str(pt+1)
	for ensembleID in eqtlDict[i]:
	    newLine=[chrID, pt1, pt2, ensembleID]
	    output1.write('\t'.join(newLine)+'\n')
    output1.close()

    runLine='intersectBed -a '+DIR+'/BED/GWAS_eQTL.common.bed -b '+totalPeak_file+' -wa -wb > '+DIR+'/BED/GWAS_eQTL.peakIntersect.bed'
    print runLine
    os.system(runLine)
    
    eqtlTargetDict={}
    eqtlInfoDict={}
    input1=open(DIR+'/BED/GWAS_eQTL.peakIntersect.bed','r')
    all_input1=input1.readlines()
    for line in all_input1:
	each=line.strip().split('\t')
	peakID=each[4]+':'+each[5]+'-'+each[6]
	eqtlID=each[0]+':'+each[1]
	ensembleID=each[3]

	if eqtlTargetDict.has_key(peakID):
	    tmp=eqtlTargetDict[peakID]
	    tmp.append(ensembleID)
	    eqtlTargetDict[peakID]=tmp
	else:
	    eqtlTargetDict[peakID]=[ensembleID]

	if eqtlInfoDict.has_key(peakID):
	    tmp=eqtlInfoDict[peakID]
	    if not eqtlID in tmp:
		tmp.append(eqtlID)
	    eqtlInfoDict[peakID]=tmp
	else:
	    eqtlInfoDict[peakID]=[eqtlID]

    input1.close()

    return eqtlTargetDict, eqtlInfoDict


def load_HiC():
    interDict={}
    interInfoDict={}
    for cellType in cellTypeList:
	input1=open(InterTarget_DIR+'/'+cellType+'.fullList.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    ensembleID=each[0]
	    geneID=each[1]
	    peakID=each[5]
	    peakLabel=each[7]
	    dist=each[8]
	    abcScore=each[11]
	    if peakLabel=='GwasPeak':
		if interDict.has_key(peakID):
		    tmp=interDict[peakID]
		    tmp.append(ensembleID)
		    interDict[peakID]=tmp
		else:
		    interDict[peakID]=[ensembleID]

		if interInfoDict.has_key(peakID):
		    tmp=interInfoDict[peakID]
		    if not cellType in tmp:
			tmp.append(cellType+';'+abcScore)
		    interInfoDict[peakID]=tmp
		else:
		    interInfoDict[peakID]=[cellType+';'+abcScore]

	input1.close()

    return interDict, interInfoDict


def getMatch():
    geneDict, tssDict=loadGeneInfo()
    eqtlTargetDict, eqtlInfoDict= match()
    interDict, interInfoDict=load_HiC()

    eqtlAssoList=[]
    interAssoList=[]
    commonPeakList = set(eqtlTargetDict) & set(interDict)
    for peakID in commonPeakList:

	for ensembleID in eqtlTargetDict[peakID]:
	    assoID=ensembleID+';'+peakID
	    eqtlAssoList.append(assoID)

	for ensembleID in interDict[peakID]:
	    assoID=ensembleID+';'+peakID
	    interAssoList.append(assoID)

    commonSet=set(eqtlAssoList) & set(interAssoList)
    eqtlSpec=set(eqtlAssoList)
    interSpec=set(interAssoList)

    output1=open(DIR+'/GWAS_eQTL_HiC.comparison.list.txt','w')
    output1.write('peakID\tGWAS_eQTL\tensembleID\tgeneID\tassoLabel\tdist\tabcScore\n')
    for assoID in eqtlSpec:
	ensembleID=assoID.split(';')[0]
	peakID=assoID.split(';')[1]
	label='eqtl_spec'
	dist=min(abs(int(tssDict[ensembleID]) - int(peakID.split(':')[1].split('-')[0])), abs(int(tssDict[ensembleID]) - int(peakID.split(':')[1].split('-')[1])))
	abcScore=','.join(interInfoDict[peakID])

	for eqtlID in eqtlInfoDict[peakID]:
	    newLine=[peakID, eqtlID, ensembleID, geneDict[ensembleID], label, str(dist), abcScore]
	    output1.write('\t'.join(newLine)+'\n')

    for assoID in commonSet:
	ensembleID=assoID.split(';')[0]
	peakID=assoID.split(';')[1]
	label='common'
	dist=min(abs(int(tssDict[ensembleID]) - int(peakID.split(':')[1].split('-')[0])), abs(int(tssDict[ensembleID]) - int(peakID.split(':')[1].split('-')[1])))
	abcScore=','.join(interInfoDict[peakID])

	for eqtlID in eqtlInfoDict[peakID]:
	    newLine=[peakID, eqtlID, ensembleID, geneDict[ensembleID], label, str(dist), abcScore]
	    output1.write('\t'.join(newLine)+'\n')

    for assoID in interSpec:
	ensembleID=assoID.split(';')[0]
	peakID=assoID.split(';')[1]
	label='hic_spec'
	dist=min(abs(int(tssDict[ensembleID]) - int(peakID.split(':')[1].split('-')[0])), abs(int(tssDict[ensembleID]) - int(peakID.split(':')[1].split('-')[1])))
	abcScore=','.join(interInfoDict[peakID])

	newLine=[peakID, '.', ensembleID, geneDict[ensembleID], label, str(dist), abcScore]
	output1.write('\t'.join(newLine)+'\n')

    output1.close()

    ## background
    output1=open(DIR+'/BED/GWAS_eQTL.commonPeaks.bed','w')
    for peakID in commonPeakList:
	chrID=peakID.split(':')[0]
	pt1=peakID.split(':')[1].split('-')[0]
	pt2=peakID.split(':')[1].split('-')[1]

	output1.write(chrID+'\t'+pt1+'\t'+pt2+'\n')
    output1.close()

    runLine='intersectBed -a '+DIR+'/BED/GeneTSS.'+winSize[1]+'.bed -b '+DIR+'/BED/GWAS_eQTL.commonPeaks.bed -wa -wb > '+DIR+'/BED/GWAS_eQTL.commonPeaks.GeneTSS.intersect.bed'
    print runLine
    os.system(runLine)
 
    input1=open(DIR+'/BED/GWAS_eQTL.commonPeaks.GeneTSS.intersect.bed', 'r')
    output1=open(DIR+'/BED/GWAS_eQTL.AllAsso.txt','w')
    all_input1=input1.readlines()
    for line in all_input1:
	each=line.strip().split('\t')
	ensembleID=each[3]
	geneID=each[4]
	peakID=each[6]+':'+each[7]+'-'+each[8]

	assoID=ensembleID+';'+peakID
	output1.write(assoID+'\n')
    input1.close()


def makeData():

    output1=open(DIR+'/GWAS_eQTL_HiC.comparison.count.txt','w')
    output1.write('Datatype\teqtl_spec\thic_spec\tcommon\tallAsso\n')
    #
    allcount=0
    input1=open(DIR+'/BED/GWAS_eQTL.AllAsso.txt','r')
    all_input1=input1.readlines()
    for line in all_input1:
	allcount+=1
    input1.close()

    #
    countDict={}
    input1=open(DIR+'/GWAS_eQTL_HiC.comparison.list.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	assoLabel=each[4]
	if countDict.has_key(assoLabel):
	    countDict[assoLabel]+=1
	else:
	    countDict[assoLabel]=1
    input1.close()

    # 
    newLine=['GWAS_eQTL', str(countDict['eqtl_spec']), str(countDict['hic_spec']), str(countDict['common']), str(allcount)]
    output1.write('\t'.join(newLine)+'\n')
    output1.close()




getMatch()
makeData()




