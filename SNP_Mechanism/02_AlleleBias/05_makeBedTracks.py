#!/python2


import os


DIR=os.getcwd()
GWAS_file='/home/ajl1213/Projects/PD/SciAdv_Anal/GWAS_Enrich/1_CollectSNP/LdFinal/MergedPD.allchr.txt'
VCF_DIR='/home/ajl1213/Projects/PD/SciAdv_Anal/SNP_Mechanism/01_SNP_Distrib/1_LdSNP/VCF'
dysPeak_DIR='/home/ajl1213/Projects/PD/SciAdv_Anal/Diff_cRE/1_BulkH3K27ac/BedFiles'
sampleList=[
'X4870_1NOSN',
'X4870_2NOSN',
'X4689NOSN',
'X5628_1NOSN',
'X5628_2NOSN',
'X5130NOSN',
'X4996_1NOSN',
'X4996_2NOSN',
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

cellTypeList=['DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']

os.system('mkdir '+DIR+'/BED')



def loadGWAS():
    gwasList=[]
    input1=open(GWAS_file,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        snpPos=each[0].replace('.',':')
        gwasList.append(snpPos)
    input1.close()

    return gwasList


def loadSNPs():
    snpDict={}
    infoDict={}
    for sampleID in sampleList:
        tmpList=[]
        tmpDict={}
        input1=open(VCF_DIR+'/'+sampleID+'.step2.vcf','r')
        all_input1=input1.readlines()
        for line in all_input1:
            if not line.count('#') > 0:
                each=line.strip().split('\t')
                chrID=each[0]
                pt=each[1]
                ref=each[3]
                alt=each[4]
                genotype=each[-1].split(':')[0]
                depth=each[-1].split(':')[1]
                if not chrID=='chrY':
                    if (genotype=='0/1' or genotype=='1/1'):
                        tmpList.append(chrID+':'+pt)
                        tmpDict[chrID+':'+pt]=[ref, alt, genotype, depth]
        input1.close()

        snpDict[sampleID]=tmpList
        infoDict[sampleID]=tmpDict

    return snpDict, infoDict


def matchSNP():
    gwasList=loadGWAS()
    snpDict, infoDict=loadSNPs()

    for sampleID in sampleList:
	output1=open(DIR+'/BED/'+sampleID+'.GWAS_matched.bed','w')
        intersect= set(gwasList) & set(snpDict[sampleID])
	if len(intersect) > 0:
	    for i in intersect:
		chrID=i.split(':')[0]
		pt=i.split(':')[1]
		pt1=int(pt)
		pt2=int(pt)+1

		newLine=[chrID, str(pt1), str(pt2)]
		output1.write('\t'.join(newLine)+'\n')
	output1.close()


def intersect():

    for sampleID in sampleList:
	runLine='intersectBed -a '+DIR+'/BED/'+sampleID+'.GWAS_matched.bed -b '+dysPeak_DIR+'/Down.Merged.bed -wa -wb | cut -f1,2,3 | sort -u > '+DIR+'/BED/'+sampleID+'.DownPeak.intersect.bed'
	print runLine
	os.system(runLine)

	#
	tmpList=[]
	input1=open(DIR+'/BED/'+sampleID+'.DownPeak.intersect.bed','r')
	all_input1=input1.readlines()
	for line in all_input1:
	    each=line.strip().split('\t')
	    chrID=each[0]
	    pt1=each[1]
	    pt2=each[2]
	    snpPt=chrID+':'+pt1+'-'+pt2
	    tmpList.append(snpPt)
	input1.close()

	if len(tmpList) > 0:
	    output1=open(DIR+'/BED/'+sampleID+'.GWAS_matched.DownPeak.final.bed','w')
	    for i in tmpList:
		chrID=i.split(':')[0]
		pt1=i.split(':')[1].split('-')[0]
		pt2=i.split(':')[1].split('-')[1]
		newLine=[chrID, pt1, pt2]
		output1.write('\t'.join(newLine)+'\n')
	    output1.close()



matchSNP()
intersect()




