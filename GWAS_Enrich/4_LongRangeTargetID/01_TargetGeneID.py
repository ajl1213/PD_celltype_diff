#!/home/ajl1213/anaconda2/bin/python


import os
import numpy



DIR=os.getcwd()
resolution='5kb'
gwasID='MergedPD'
GTF='/home/ajl1213/genome.info/gencode/hg19/gencode.v19.ptn_coding.lincRNA.level_1_2'
SnpCreFile='/home/ajl1213/Projects/PD/SciAdv_Anal/GWAS_Enrich/2_ClassifySNP/SnpTable/'+gwasID+'.SnpTable.txt'
SigInterFile='/home/ajl1213/Projects/PD/data/HiC/FitHiC/SigInter/MergedTotal.'+resolution+'.1e-2.FitHiC.ProCre.txt'
NearbyWindow=[5000,'5kb']

os.system('mkdir '+DIR+'/GeneLocBED')
os.system('mkdir '+DIR+'/GwasCreBED')
os.system('mkdir '+DIR+'/NearbyTargetList')
os.system('mkdir '+DIR+'/LongRangeTargetList')
os.system('mkdir '+DIR+'/MergedTargetList')


def loadGeneInfo():
    geneDict={}
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
	    output1.write(chrID+'\t'+tssPos+'\t'+str(int(tssPos)+1)+'\t'+ensembleID+'\t'+geneID+'\n')
    input1.close()

    return geneDict


def loadGwasCre():
    gwasCreDict={}

    input1=open(SnpCreFile,'r')
    all_input1=input1.readlines()
    for line in all_input1:
	each=line.strip().split('\t')
	creID=each[2]
	snpPos=each[0]
	rsID=each[1]

	if gwasCreDict.has_key(gwasID):
	    tmp=gwasCreDict[gwasID]
	    if not creID in tmp:
		tmp.append(creID)
	    gwasCreDict[gwasID]=tmp
	else:
	    gwasCreDict[gwasID]=[creID]
    input1.close()
    
    return gwasCreDict


def getNearbyTarget():
    geneDict=loadGeneInfo()
    gwasCreDict=loadGwasCre()

    print 'Number of gwas harboring cREs'
    print gwasID, len(gwasCreDict[gwasID])

    output1=open(DIR+'/GwasCreBED/'+gwasID+'.'+NearbyWindow[1]+'.bed','w')
    for creID in gwasCreDict[gwasID]:
	chrID=creID.split(':')[0]
	pt1=int(creID.split(':')[1].split('-')[0])
	pt2=int(creID.split(':')[1].split('-')[1])

	newLine=[chrID, str(pt1-NearbyWindow[0]), str(pt2+NearbyWindow[0]), creID]
	output1.write('\t'.join(newLine)+'\n')
    output1.close()

    runLine='intersectBed -a '+DIR+'/GwasCreBED/'+gwasID+'.'+NearbyWindow[1]+'.bed -b '+DIR+'/GeneLocBED/GeneTss.bed -wa -wb > '+DIR+'/GwasCreBED/'+gwasID+'.geneTssIntersect.'+NearbyWindow[1]+'.bed'
    os.system(runLine)

    tmpDict={}
    input1=open(DIR+'/GwasCreBED/'+gwasID+'.geneTssIntersect.'+NearbyWindow[1]+'.bed','r')
    all_input1=input1.readlines()
    for line in all_input1:
	each=line.strip().split('\t')
	creID=each[3]
	ensembleID=each[7]
	if tmpDict.has_key(ensembleID):
	    tmp=tmpDict[ensembleID]
	    tmp.append(creID)
	    tmpDict[ensembleID]=tmp
	else:
	    tmpDict[ensembleID]=[creID]
    input1.close()

    output1=open(DIR+'/NearbyTargetList/'+gwasID+'.NearbyTargetList.'+NearbyWindow[1]+'.txt','w')
    output1.write('ensembleID\tgeneID\tcreList\n')
    for ensembleID in tmpDict:
	if geneDict.has_key(ensembleID):
	    geneID=geneDict[ensembleID]
	    creList=','.join(set(tmpDict[ensembleID]))
	    output1.write(ensembleID+'\t'+geneID+'\t'+creList+'\n')
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

            for geneID in geneList.split(';'):
                for creID in creList.split(';'):
                    if interDict.has_key(creID):
                        tmp=interDict[creID]
                        tmp.append(geneID)
                        interDict[creID]=tmp
                    else:
                        interDict[creID]=[geneID]
    input1.close()

    return interDict


def getLongRangeTarget():
    geneDict=loadGeneInfo()
    gwasCreDict=loadGwasCre()
    interDict=loadSigInter()

    tmpDict={}
    for creID in gwasCreDict[gwasID]:
	if interDict.has_key(creID):
	    for geneID in interDict[creID]:

		if tmpDict.has_key(geneID):
		    tmp=tmpDict[geneID]
		    tmp.append(creID)
		    tmpDict[geneID]=tmp
		else:
		    tmpDict[geneID]=[creID]

    output1=open(DIR+'/LongRangeTargetList/'+gwasID+'.LongRangeTargetList.txt','w')
    output1.write('ensembleID\tgeneID\tcreList\n')
    for ensembleID in tmpDict:
	if geneDict.has_key(ensembleID):
	    geneID=geneDict[ensembleID]
	    creList=','.join(set(tmpDict[ensembleID]))
	    output1.write(ensembleID+'\t'+geneID+'\t'+creList+'\n')
    output1.close()


def mergeTargetList():
    geneDict=loadGeneInfo()
    gwasCreDict=loadGwasCre()

    targetDict={}

    ## Nearby Target
    input1=open(DIR+'/NearbyTargetList/'+gwasID+'.NearbyTargetList.'+NearbyWindow[1]+'.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	ensembleID=each[0]
	targetDict[ensembleID]='NearbyTarget'
    input1.close()

    ## Long-range Target
    input1=open(DIR+'/LongRangeTargetList/'+gwasID+'.LongRangeTargetList.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	ensembleID=each[0]
	if not targetDict.has_key(ensembleID):
	    targetDict[ensembleID]='LongRangeTarget'
    input1.close()

    ##
    output1=open(DIR+'/MergedTargetList/'+gwasID+'.MergedTargetList.txt','w')
    output1.write('ensembleID\tgeneID\tdataType\n')
    for ensembleID in targetDict:
	geneID=geneDict[ensembleID]
	dataType=targetDict[ensembleID]
	newLine=[ensembleID, geneID, dataType]
	output1.write('\t'.join(newLine)+'\n')
    output1.close()






getNearbyTarget()
getLongRangeTarget()
mergeTargetList()



