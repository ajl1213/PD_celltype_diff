#!/home/ajl1213/anaconda2/bin/python


import os
import numpy



DIR=os.getcwd()
resolution='5kb'
NearbyWindow=[5000,'5kb']
GTF='/home/ajl1213/genome.info/gencode/hg19/gencode.v19.ptn_coding.lincRNA.level_1_2'
DysPeakDIR='/home/ajl1213/Projects/PD/SciAdv_Anal/Diff_cRE/1_BulkH3K27ac/BedFiles'
SnpPeakFile='/home/ajl1213/Projects/PD/SciAdv_Anal/GWAS_Enrich/2_ClassifySNP/SnpTable/MergedPD.SnpTable.txt'
SigInterFile='/home/ajl1213/Projects/PD/SciAdv_Data/HiC/FitHiC/SigInter/MergedTotal.'+resolution+'.1e-2.FitHiC.ProCre.txt'
dataTypeList=['Down','Up']
cellTypeList=['DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']

os.system('mkdir '+DIR+'/GeneLocBED')
os.system('mkdir '+DIR+'/DysCreBED')
os.system('mkdir '+DIR+'/NearbyTargetList')
os.system('mkdir '+DIR+'/LongRangeTargetList')
os.system('mkdir '+DIR+'/MergedTargetList')



def loadGeneInfo():
    geneDict={}
    tssDict={}
    input1=open(GTF,'r')
    all_input1=input1.readlines()
    output1=open(DIR+'/GeneLocBED/GeneTss.bed','w')
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        ensembleID=each[0]
        geneID=each[1]
        chrID=each[2]
        tssPos=each[3]
        if each[6]=='protein_coding':
            geneDict[ensembleID]=geneID
            tssDict[ensembleID]=chrID+'.'+tssPos
            output1.write(chrID+'\t'+tssPos+'\t'+str(int(tssPos)+1)+'\t'+ensembleID+'\t'+geneID+'\n')
    input1.close()

    return geneDict, tssDict


def loadPeaks():
    # dys cRE
    dysPeakDict={}
    for cellType in cellTypeList:
	for dataType in dataTypeList:
	    tmpList=[]
            input1=open(DysPeakDIR+'/'+dataType+'.'+cellType+'.bed','r')
            all_input1=input1.readlines()
            for line in all_input1:
                each=line.strip().split('\t')
                chrID=each[0]
                pt1=each[1]
                pt2=each[2]
                peakID=each[3]
		tmpList.append(peakID)
            input1.close()

	    dysPeakDict[dataType+';'+cellType]=tmpList

    # GWAS cRE
    gwasPeakList=[]
    input1=open(SnpPeakFile,'r')
    all_input1=input1.readlines()
    for line in all_input1:
        each=line.strip().split('\t')
        peakID=each[2]
        snpPos=each[0]
        rsID=each[1]
	gwasPeakList.append(peakID)
    input1.close()

    return dysPeakDict, gwasPeakList


def getNearbyTarget():
    geneDict, tssDict = loadGeneInfo()
    dysPeakDict, gwasPeakList = loadPeaks()

    for cellType in cellTypeList:
	output1=open(DIR+'/DysCreBED/'+cellType+'.'+NearbyWindow[1]+'.bed','w')
	#
	for dataType in dataTypeList:
            for creID in dysPeakDict[dataType+';'+cellType]:
                chrID=creID.split(':')[0]
                pt1=int(creID.split(':')[1].split('-')[0])
                pt2=int(creID.split(':')[1].split('-')[1])

		newLine=[chrID, str(pt1-NearbyWindow[0]), str(pt2+NearbyWindow[0]), creID, dataType+'Peak']
		output1.write('\t'.join(newLine)+'\n')
	#
	for creID in gwasPeakList:
	    chrID=creID.split(':')[0]
	    pt1=int(creID.split(':')[1].split('-')[0])
	    pt2=int(creID.split(':')[1].split('-')[1])

	    newLine=[chrID, str(pt1-NearbyWindow[0]), str(pt2+NearbyWindow[0]), creID, 'GwasPeak']
	    output1.write('\t'.join(newLine)+'\n')
	output1.close()

	#
	runLine='intersectBed -a '+DIR+'/DysCreBED/'+cellType+'.'+NearbyWindow[1]+'.bed -b '+DIR+'/GeneLocBED/GeneTss.bed -wa -wb > '+DIR+'/DysCreBED/'+cellType+'.geneTssIntersect.'+NearbyWindow[1]+'.bed'
	print runLine
	os.system(runLine)

	#
    	elementList=[]
	input1=open(DIR+'/DysCreBED/'+cellType+'.geneTssIntersect.'+NearbyWindow[1]+'.bed','r')
	all_input1=input1.readlines()
	for line in all_input1:
	    each=line.strip().split('\t')
	    creID=each[3]
	    creType=each[4]
	    ensembleID=each[8]
	    valID=[ensembleID, creID, creType]
	    elementList.append(';'.join(valID))
	input1.close()

	#
	output1=open(DIR+'/NearbyTargetList/'+cellType+'.NearbyTargetList.'+NearbyWindow[1]+'.txt','w')
	output1.write('ensembleID\tgeneID\tcreID\tcreType\n')
	for i in set(elementList):
	    [ensembleID, creID, creType]=i.split(';')
	    geneID=geneDict[ensembleID]
	    newLine=[ensembleID, geneID, creID, creType]
	    output1.write('\t'.join(newLine)+'\n')
	output1.close()


def loadSigInter():
    interDict={}
    input1=open(SigInterFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        interLabel=each[0]
        if interLabel=='PC':
            frag1=each[1]
            geneList=each[2]
            frag2=each[3]
            creList=each[4]

            for ensembleID in geneList.split(';'):
                for creID in creList.split(';'):
                    if interDict.has_key(creID):
                        tmp=interDict[creID]
                        tmp.append(ensembleID)
                        interDict[creID]=tmp
                    else:
                        interDict[creID]=[ensembleID]
    input1.close()

    return interDict


def getLongRangeTarget():
    geneDict, tssDict = loadGeneInfo()
    dysPeakDict, gwasPeakList = loadPeaks()
    interDict = loadSigInter()

    for cellType in cellTypeList:
	elementList=[]
	#
	for dataType in dataTypeList:
            for creID in dysPeakDict[dataType+';'+cellType]:
                if interDict.has_key(creID):
                    for ensembleID in interDict[creID]:
			valID=[ensembleID, creID, dataType+'Peak']
			elementList.append(';'.join(valID))
	#
	for creID in gwasPeakList:
	    if interDict.has_key(creID):
		for ensembleID in interDict[creID]:
		    valID=[ensembleID, creID, 'GwasPeak']
		    elementList.append(';'.join(valID))

	#
	output1=open(DIR+'/LongRangeTargetList/'+cellType+'.LongRangeTargetList.txt','w')
	output1.write('ensembleID\tgeneID\tcreID\tcreType\n')
	for i in set(elementList):
	    [ensembleID, creID, creType]=i.split(';')
	    if geneDict.has_key(ensembleID):
		geneID=geneDict[ensembleID]
		newLine=[ensembleID, geneID, creID, creType]
		output1.write('\t'.join(newLine)+'\n')
	output1.close()


def mergeTargetList():
    geneDict, tssDict = loadGeneInfo()
    dysPeakDict, gwasPeakList = loadPeaks()

    for cellType in cellTypeList:
	elementList=[]
	#
	input1=open(DIR+'/NearbyTargetList/'+cellType+'.NearbyTargetList.'+NearbyWindow[1]+'.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    elementList.append(';'.join(each))
	input1.close()

	#
	input1=open(DIR+'/LongRangeTargetList/'+cellType+'.LongRangeTargetList.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    elementList.append(';'.join(each))
	input1.close()

	#
	output1=open(DIR+'/MergedTargetList/'+cellType+'.MergedTargetList.txt','w')
	output1.write('ensembleID\tgeneID\tcreID\tcreType\tdist\n')
	for i in set(elementList):
	    [ensembleID, geneID, creID, creType]=i.split(';')

	    #
	    tss=tssDict[ensembleID]
	    pt1=int(creID.split(':')[1].split('-')[0])
	    pt2=int(creID.split(':')[1].split('-')[1])
	    tssPt=int(tss.split('.')[1])

	    pt=(pt1+pt2)/2
	    dist=abs(pt-tssPt)

	    newLine=[ensembleID, geneID, creID, creType, str(dist)]
	    output1.write('\t'.join(newLine)+'\n')
	output1.close()



getNearbyTarget()
getLongRangeTarget()
mergeTargetList()



