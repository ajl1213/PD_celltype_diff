#!/python


import os


DIR=os.getcwd()
rawGWAS_file='/home/ajl1213/Projects/PD/SciAdv_Anal/GWAS_Enrich/1_CollectSNP/LdFinal/MergedPD.allchr.txt'
sampleList=[
'X4870_1NOSN',
'X4870_2NOSN',
'X4689NOSN',
'X5628_1NOSN',
'X5628_2NOSN',
'X5130NOSN',
'X4996_1NOSN',
'X4996_2NOSN',
'X5006NOSN',
'X5244PDSN',
'X4831PDSN',
'X5532PDSN',
'X5215PDSN',
'X5627PDSN',
'X5778PDSN',
'X5742PDSN',
'X5649PDSN',
'X5591PDSN'
]


def loadGWAS():
    gwasList=[]
    input1=open(rawGWAS_file,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	snpPos=each[0].replace('.',':')
	gwasList.append(snpPos)
    input1.close()

    return gwasList


def loadSNPs():
    snpDict={}
    for sampleID in sampleList:
	tmpList=[]
	input1=open(DIR+'/VCF/'+sampleID+'.step2.vcf','r')
	all_input1=input1.readlines()
	for line in all_input1:
	    if not line.count('#') > 0:
		each=line.strip().split('\t')

		chrID=each[0]
		pt=each[1]
		ref=each[3]
		alt=each[4]

		genotype=each[-1].split(':')[0]

		if (genotype=='1/1' or genotype=='0/1'):
		    tmpList.append(chrID+':'+pt)
	input1.close()

	snpDict[sampleID]=tmpList

    return snpDict


def compare():
    gwasList=loadGWAS()
    snpDict=loadSNPs()

    output1=open(DIR+'/PD_GWAS.match.txt','w')
    output1.write('SampleID\tGWAS_match\tGWAS_frac\n')

    output2=open(DIR+'/MatchedList.txt','w')
    for sampleID in sampleList:

	intersect= set(gwasList) & set(snpDict[sampleID])
	frac=float(len(intersect)) / float(len(snpDict[sampleID])) * 1000

	output1.write(sampleID+'\t'+str(len(intersect)) +'\t'+str(frac)+'\n')

	for i in intersect:
	    output2.write(sampleID+'\t'+i.replace(':','.')+'\n')

    output1.close()
    output2.close()




compare()



