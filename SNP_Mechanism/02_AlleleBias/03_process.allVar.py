#!/python2


import os
import pysam
from scipy import stats
import numpy
import random


DIR=os.getcwd()
VCF_DIR='/home/ajl1213/Projects/PD/SciAdv_Anal/SNP_Mechanism/01_SNP_Distrib/1_LdSNP/VCF'
BAM_DIR='/home/ajl1213/Projects/PD/SciAdv_Data/ChIP/01_Mapping/SortedNodupFilteredBam'
FASTA = pysam.Fastafile('/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa')
FAI='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa.fai'
nSampling=5000
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

os.system('mkdir '+DIR+'/SNP_Reads')



def get_heteroSNPs():
    infoDict={}
    allList=[]
    for sampleID in sampleList:
	if sampleID.count('PD') > 0:
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

		    if not chrID == 'chrY':
			if genotype=='0/1':
			    allList.append(sampleID+':'+chrID+':'+pt)
			    tmpDict[chrID+':'+pt]=[ref, alt, genotype, depth]

	    input1.close()

	    infoDict[sampleID]=tmpDict

    #
    matchDict={}
    output1=open(DIR+'/hetero_Var.txt','w')
    output1.write('sampleID\tchrID\tpt\tref\talt\tgenotype\tdepth\n')
    for i in random.sample(allList, nSampling):
	sampleID=i.split(':')[0]
	chrID=i.split(':')[1]
	pt=i.split(':')[2]
	pt1=int(pt)-1
	pt2=int(pt)+1

	snpID=chrID+':'+str(pt1)+'-'+str(pt2)
	if matchDict.has_key(sampleID):
	    tmp=matchDict[sampleID]
	    tmp.append(snpID)
	    matchDict[sampleID]=tmp
	else:
	    matchDict[sampleID]=[snpID]

	[ref, alt, genotype, depth]=infoDict[sampleID][chrID+':'+pt]
	newLine=[sampleID, chrID, pt, ref, alt, genotype, depth]
	output1.write('\t'.join(newLine)+'\n')
    output1.close()

    return matchDict


def BamIndex():
    ncpu=2
    for sampleID in sampleList:
	output1=open('PBS_bamidx_'+sampleID+'.pbs','w')
	output1.write('#PBS -N bamidx_'+sampleID+'\n')
	output1.write('#PBS -q workq\n')
	output1.write('#PBS -l nodes=1:ppn='+str(ncpu)+'\n')
	output1.write('#PBS -j oe\n')
	output1.write('\n')
	output1.write('# go workdir\n')
	output1.write('cd $PBS_O_WORKDIR\n')
	output1.write('\n')
	output1.write('# run command \n')
	output1.write('sleep 5\n')
	output1.write('\n')
	output1.write('echo -n \"I am on: \"\n')
	output1.write('hostname;\n')    
	output1.write('echo finding ssh-accessible nodes:\n')
	output1.write('echo -n \"running on: \"\n')
	output1.write('\n')

	output1.write('samtools index '+BAM_DIR+'/'+sampleID+'.H3K27ac.sorted.nodup.filtered.bam\n')

	output1.write('\n')
	output1.write('sleep 30\n')
	output1.write('exit 0')
	output1.close()

	os.system('qsub PBS_bamidx_'+sampleID+'.pbs')


def collectReads():
    matchDict=get_heteroSNPs()

    for sampleID in sampleList:
	if sampleID.count('PD') > 0:
	    for snpID in matchDict[sampleID]:
		runLine='samtools view '+BAM_DIR+'/'+sampleID+'.H3K27ac.sorted.nodup.filtered.bam \''+snpID+'\' > '+DIR+'/SNP_Reads/'+sampleID+'.'+snpID.replace(':','_').replace('-','_')+'.allVar.txt'
		print runLine
		os.system(runLine)



def collectSeq():
    # load query allele
    queryDict={}
    input1=open(DIR+'/hetero_Var.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	sampleID=each[0]
	chrID=each[1]
	pt=each[2]
	ref=each[3]
	alt=each[4]
	queryDict[sampleID+';'+chrID+';'+pt]=[ref, alt]
    input1.close() 

    #
    output1=open(DIR+'/AllelicReads.allVar.txt','w')
    output1.write('queryID\tfastaSeq\trefSeq\taltSeq\trefReads\taltReads\totherReads\n')
    for i in queryDict:
	[sampleID, chrID, query_pt]=i.split(';')
	[ref, alt]=queryDict[i]

	query_pt=int(query_pt)

	# FASTA check
	fasta_seq = FASTA.fetch(chrID, query_pt-1, query_pt).upper()

	# read check
	pt1=str(query_pt-1)
	pt2=str(query_pt+1)
	input1=open('SNP_Reads/'+sampleID+'.'+chrID+'_'+pt1+'_'+pt2+'.allVar.txt','r')
	all_input1=input1.readlines()
	countDict={}
	for line in all_input1:
	    each=line.strip().split('\t')
	    readID=each[0]
	    chrID=each[2]
	    read_pt=int(each[3])

	    nPos=query_pt-read_pt

	    seqList=each[9]
	    if len(list(seqList)) > (nPos):
		querySeq=list(seqList)[nPos]

		if querySeq == ref:
		    seqType='ref'
		elif querySeq == alt:
		    seqType='alt'
		else:
		    seqType='other'

		if countDict.has_key(seqType):
		    countDict[seqType]+=1
		else:
		    countDict[seqType]=1	
	    else:
		'error'
	input1.close()

	if countDict.has_key('ref'):
	    ref_count=countDict['ref']
	else:
	    ref_count=0

	if countDict.has_key('alt'):
	    alt_count=countDict['alt']
	else:
	    alt_count=0

	if countDict.has_key('other'):
	    other_count=countDict['other']
	else:
	    other_count=0

	newLine=[sampleID+';'+chrID+';'+str(query_pt), fasta_seq, ref, alt, str(ref_count), str(alt_count), str(other_count)]
	output1.write('\t'.join(newLine)+'\n')

    output1.close()


def calcBinomP():
    input1=open(DIR+'/AllelicReads.allVar.txt', 'r')
    all_input1=input1.readlines()
    output1=open(DIR+'/AllelicReads.allVar.pval.txt','w')
    output1.write('queryID\tfastaSeq\trefSeq\taltSeq\trefReads\taltReads\totherReads\tallReads\taltExp\tlog2Enrich\tbinomP\tlabel\n')
    for line in all_input1[1:]:
	[queryID, fastaSeq, refSeq, altSeq, refReads, altReads, otherReads]=line.strip().split('\t')

	allReads=int(refReads)+int(altReads)
	altExp= float(allReads) / 2
	log2Enrich= numpy.log2((float(altReads)+1)/(altExp+1))
	binomP = stats.binom_test(x=altReads, n=allReads, p=0.5)

	if abs(log2Enrich) > 0.5:
	    if log2Enrich < 0:
		label= 'refEnrich'
	    elif log2Enrich >=0:
		label= 'altEnrich'

	else:
	    label='unchanged'

	newLine=[queryID, fastaSeq, refSeq, altSeq, refReads, altReads, otherReads, str(allReads), str(altExp), str(log2Enrich), str(binomP), label]
	output1.write('\t'.join(newLine)+'\n')

    input1.close()
    output1.close()




#get_heteroSNPs()
#BamIndex()
#collectReads()
collectSeq()
calcBinomP()


