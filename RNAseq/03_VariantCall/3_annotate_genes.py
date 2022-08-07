#!/python2


import os

# PMID: 30448284

DIR=os.getcwd()

sampleList=[
'X4870NOSN',
'X4689NOSN',
'X5628NOSN',
'X5130_1NOSN',
'X5130_2NOSN',
'X4996NOSN',
'X5006_1NOSN',
'X5006_2NOSN',
'X5244PDSN',
'X4831PDSN',
'X5215PDSN',
'X5627PDSN',
'X5778PDSN',
'X5742PDSN',
'X5649PDSN',
'X5591PDSN'
]

geneList=[
'APOE',
'SNCA',
'MAPT',
'LRRK2',
'PRKN', # Parkin
'PINK1',
'PARK7', #DJ1
'GBA'
]

annoDict={}
for sampleID in sampleList:
    print sampleID
    input1=open(DIR+'/annovar/'+sampleID+'.avinput.exonic_variant_function','r')
    all_input1=input1.readlines()

    tmpDict={}
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	varType=each[1]
	varAnno=each[2].split(',')[0]
	geneID=varAnno.split(':')[0]
	varInfo=each[20]
	genotype=varInfo.split(':')[0]
	rsID=each[13]
	#valID=rsID+'('+genotype+')'

	if not varType=='synonymous SNV':

	    if varType in ['nonsynonymous SNV','unknown']:
		varType_final='missense'
	    else:
		varType_final=varType
		
	    if (genotype=='0/1' or genotype=='1/1'):

		valID=rsID+'('+varType_final+','+genotype+')'

		if tmpDict.has_key(geneID):
		    tmp=tmpDict[geneID]
		    if not valID in tmp:
			tmp.append(valID)
		    tmpDict[geneID]=tmp
		else:
		    tmpDict[geneID]=[valID]
    input1.close()

    annoDict[sampleID]=tmpDict
    

output1=open(DIR+'/PD.genotyping.txt','w')
output1.write('SampleID\t'+'\t'.join(geneList)+'\n')
for sampleID in sampleList:
    tmpList=[]
    for geneID in geneList:
	if annoDict[sampleID].has_key(geneID):
	    a=','.join(annoDict[sampleID][geneID])
	else:
	    a='.'
	tmpList.append(a)

    output1.write(sampleID+'\t'+'\t'.join(tmpList)+'\n')
output1.close()




