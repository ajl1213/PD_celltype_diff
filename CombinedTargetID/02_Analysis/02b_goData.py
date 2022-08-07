#!/home/ajl1213/anaconda2/bin/python


import os
import numpy
import operator



DIR=os.getcwd()
cellTypeList=['DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']

os.system('mkdir '+DIR+'/GOTABLESORTED')


## BP (enrichr)
for cellType in cellTypeList:
    nDict={}
    allDict={}
    input1=open(DIR+'/GOTABLE/'+cellType+'.BP.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	
	termID=each[0]
	nGene=float(each[1].split('/')[0])
	pval=each[2]
	adjPval=each[3]
	oddRatio=each[6]
	combinedScore=each[7]
	genes=each[8]
	cellLabel=each[9]
	dataType=each[10]

	nDict[termID]=nGene
	allDict[termID]=[pval, dataType, genes]

    input1.close()

    output1=open(DIR+'/GOTABLESORTED/'+cellType+'.enrichr.BP.txt','w')
    output1.write('termID\tdataType\tpval\tnGene\tgeneList\n')
    for key in sorted(nDict.items(), key=operator.itemgetter(1), reverse=True):
	termID=key[0]
	nGene=key[1]
	pval=allDict[termID][0]
	dataType=allDict[termID][1]
	genes=allDict[termID][2]
	newLine=[termID, dataType, pval, str(nGene), genes]
	output1.write('\t'.join(newLine)+'\n')
    output1.close()

## BP (metascape)
for cellType in ['DopaN','Oligo','Ast','Micro']:
    nDict={}
    allDict={}
    for dataType in ['Down','Up']:
	input1=open(DIR+'/GOTABLE/metascape_'+cellType+'.'+dataType+'.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    if each[0].count('Summary') > 0:
		termID = each[3]
		logPval=each[4]
		genes=each[8].replace('\"', '')
		nGene=len(genes.split(','))

	    nDict[termID]=nGene
	    allDict[termID]=[logPval, dataType, genes]
    input1.close()

    output1=open(DIR+'/GOTABLESORTED/'+cellType+'.metascape.BP.txt','w')
    output1.write('termID\tdataType\tpval\tnGene\tgeneList\n')
    for key in sorted(nDict.items(), key=operator.itemgetter(1), reverse=True):
	termID=key[0]
	nGene=key[1]
	pval=allDict[termID][0]
	dataType=allDict[termID][1]
	genes=allDict[termID][2]
	newLine=[termID, dataType, pval, str(nGene), genes]
	output1.write('\t'.join(newLine)+'\n')
    output1.close()


#=====================================================================================
## Oligo Down
# myelination     Down    -2.273680937    4       ASPA,DEGS1,MTMR2,ACER3 
# proteasomal protein catabolic process   Down    -2.37650727     11      PSMC6,UBE2G1,CDC16,TRIM13,RNF139,FBXL5,HECW2,RFFL,PDCD6IP,USP33,FBXO28
#
## Oligo Up
# autophagy       Up      -3.223652183    9       ATM,MARK2,ATP6V0D1,TGFBRAP1,CALCOCO2,ATG2A,TOMM7,VPS18,MTMR14
# regulation of microtubule cytoskeleton organization     Up      -3.251159046    9       MARK2,FKBP4,NUMA1,PARP3,KIF15,DST,C2CD3,TTL,TUBB2B
# protein folding Up      -2.568302762    5       FKBP4,PFDN6,PTGES3,RUVBL2,DNAJB4
#
## Ast Down
#brain morphogenesis (GO:0048854)        2/12    0.000779106356037233    0.119769713728356       0       0       58.5882352941176        419.337267088932        BBS2;PAFAH1B1   Ast     Down    GO_Biological_Process_2021
#protein ubiquitination (GO:0016567)     7/525   0.00237279961700305     0.170033071269187       0       0       4.16387816387816        25.1651669515274        RNF138;HLTF;TRIM13;PEX2;UBE2G1;SKP2;SH3RF1      Ast     Down    GO_Biological_Process_2021
# 
## Ast Up
# autophagy       Up      -3.31602839     13      MARK2,SUPT5H,CALCOCO2,ATG2A,ARFIP2,VPS13D,MTMR14,ZC3H12A,SBF2,KLHL22,SMARCC2,WDR1,PPP1R9B
# cellular response to lipid      Up      -2.468734976    9       ZFP36,PTGES3,RUVBL2,ABHD2,ADNP2,PAF1,ZC3H12A,PPP1R9B,SGMS1
# protein folding Up      -3.856389441    7       DNAJC7,PTGES3,RUVBL2,HSPA4L,DNAJB5,LMAN2L,FUT10
#
## DopaN Down
# learning or memory      Down    -2.691905812    3       PAFAH1B1,ATP8A1,ARL6IP5
# establishment or maintenance of cell polarity   Down    -4.191282063    9       PAFAH1B1,PDCD6IP,CLASP2,RICTOR,CDC16,TTC19,HECW2,CFL2,TTLL7
#
## DopaN Up
# regulation of phosphoprotein phosphatase activity       Up      -2.705233884    14      HSP90AB1,PPP1R11,BOD1,ARG2,N4BP1,RNF216,CUEDC2,ERLIN1,FBXO8,USP35,KLHL22,PPP2R5B,TIMM50,INPP5B
# autophagy       Up      -3.343848879    10      ITPR1,ATP6V0D1,HAX1,ATG2A,TOMM7,VPS18,SBF2,KLHL22,VTI1A,CACNG7
# protein folding Up      -5.393894982    9       FKBP4,DNAJA1,HSP90AB1,DNAJC7,PFDN6,CCT8,RUVBL2,HSPA4L,ABCA7
#=====================================================================================
    

## MP
nDict={}
allDict={}
input1=open(DIR+'/GOTABLE/Merged.MP.txt','r')
all_input1=input1.readlines()
for line in all_input1[1:]:
    each=line.strip().split('\t')
    termID=' '.join(each[0].split(' ')[:-1])
    nGene=float(each[1].split('/')[0])
    pval=each[2]
    adjPval=each[3]
    oddRatio=each[6]
    combinedScore=each[7]
    genes=each[8]

    nDict[termID]=nGene
    allDict[termID]=[pval, oddRatio, combinedScore, genes]
input1.close()

output1=open(DIR+'/GOTABLESORTED/Merged.MP.txt','w')
output1.write('termID\tpval\tOR\tcombinedScore\tnGene\tgeneList\n')
for key in sorted(nDict.items(), key=operator.itemgetter(1), reverse=True):
    termID=key[0]
    nGene=key[1]
    pval=allDict[termID][0]
    oddRatio=allDict[termID][1]
    combinedScore=allDict[termID][2]
    genes=allDict[termID][3]

    if not termID.count('lethality') > 0:
	if not termID.count('death') > 0:
	    if not termID.count('no abnormal phenotype') > 0:

		if int(nGene) >= 8:
		    if float(pval) < 0.05:
			    newLine=[termID, pval, oddRatio, combinedScore, str(nGene), genes]
			    output1.write('\t'.join(newLine)+'\n')
output1.close()

#
termList=[
'hyperactivity',
'hypoactivity',
'tremors',
'decreased grip strength',
'abnormal gait',
'limb grasping',
'impaired coordination',
'decreased startle reflex',
'decreased prepulse inhibition',
'ataxia',
'abnormal motor capabilities/coordination/movement',
'abnormal neural tube morphology',
'abnormal brain morphology',
'decreased vertical activity',
'decreased lymphocyte cell number',
'abnormal locomotor behavior',
'decreased CD8-positive, alpha-beta T cell number',
'decreased T cell number',
'abnormal axon morphology',
'abnormal myelination',
'increased monocyte cell number',
'abnormal locomotor activation',
'demyelination',
'astrocytosis',
'increased macrophage cell number',
'impaired righting response',
'increased neuron apoptosis',
'abnormal neural tube closure'
]


print len(termList)

geneList=[]
input1=open(DIR+'/GOTABLESORTED/Merged.MP.txt','r')
all_input1=input1.readlines()
for line in all_input1[1:]:
    each=line.strip().split('\t')
    termID=each[0]
    if termID in termList:
        for geneID in each[3].split(';'):
            geneList.append(geneID)
input1.close()

print len(set(geneList))
print set(geneList)







