#!/python2


import os
import operator
import gzip


DIR=os.getcwd()
code_DIR='/home/ajl1213/programs/ANNOVAR/annovar'
inputVCF_DIR=DIR+'/freebayes'
FASTA='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa'
dbSNP_file='/home/ajl1213/genome.info/dbSNP/dbsnp_b151_GRCh37p13.vcf.gz'

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

os.system('mkdir '+DIR+'/annovar')
os.system('mkdir '+DIR+'/dbSNP')


def processVCF():

    print 'step 1 / step 2'
    for sampleID in sampleList:
	print sampleID

	## command for downlading appropriate database into humandb directory
	# annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/

	## Process vcf file according to ANNOVAR suggestion
	## Step1: split vcf lines so that each line contains one and only one variant (originall vcf contains multiple variants in a single line)
	## Step2: peform left-normalization (A way to process indel)
   
	run_line='bcftools norm -m-both -o '+DIR+'/annovar/'+sampleID+'.step1.vcf '+inputVCF_DIR+'/'+sampleID+'.vcf'
	print run_line
	os.system(run_line)

	run_line='bcftools norm -f '+FASTA+' -o '+DIR+'/annovar/'+sampleID+'.step2.vcf '+DIR+'/annovar/'+sampleID+'.step1.vcf'
	print run_line
	os.system(run_line)


    ## Step3: conduct filter
    ## Since its "Minor" allele frequency, remove variants with more than 50% occurrence.
    print 'step 3'
    allDict={}
    for sampleID in sampleList:
	print sampleID
	input1=open(DIR+'/annovar/'+sampleID+'.step2.vcf','r')
	all_input1=input1.readlines()
	for line in all_input1:
	    if not line.count('#') > 0:
		each=line.strip().split('\t')

		chrID=each[0]
		pt=each[1]
		ref=each[3]
		alt=each[4]

		keyID=chrID+';'+pt+';'+ref+';'+alt
		if allDict.has_key(keyID):
		    allDict[keyID]+=1
		else:
		    allDict[keyID]=1
	input1.close()

    for sampleID in sampleList:
	input1=open(DIR+'/annovar/'+sampleID+'.step2.vcf','r')
	output1=open(DIR+'/annovar/'+sampleID+'.step3.vcf','w')
	all_input1=input1.readlines()
	for line in all_input1:
	    if not line.count('#') > 0:
		each=line.strip().split('\t')

		chrID=each[0]
		pt=each[1]
		ref=each[3]
		alt=each[4]

		keyID=chrID+';'+pt+';'+ref+';'+alt
		if not allDict[keyID] > len(sampleList)/2: ##
		    output1.write('\t'.join(each)+'\n')
	    else:
		output1.write(line)

	input1.close()
	output1.close()


def match_dbSNP():

    for sampleID in sampleList:
	input1=open(DIR+'/annovar/'+sampleID+'.step3.vcf','r')
	output1=open(DIR+'/annovar/'+sampleID+'.step3.mod.vcf','w')
	all_input1=input1.readlines()
	for line in all_input1:
	    if not line.count('#') > 0:
		each=line.strip().split('\t')

		each[0]=each[0].split('chr')[1]
		output1.write('\t'.join(each)+'\n')
	    else:
		output1.write(line)

	input1.close()
	output1.close()

	runLine='bgzip '+DIR+'/annovar/'+sampleID+'.step3.mod.vcf'
	print runLine
	os.system(runLine)

	runLine='tabix -p vcf '+DIR+'/annovar/'+sampleID+'.step3.mod.vcf.gz'
	print runLine
	os.system(runLine)


    #
    ncpu=2
    for sampleID in sampleList:
	output1=open('PBS_bcfisec_'+sampleID+'.pbs','w')
	output1.write('#PBS -N bcfisec_'+sampleID+'\n')
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

	# desired command lines
	output1.write('bcftools isec '+DIR+'/annovar/'+sampleID+'.step3.mod.vcf.gz '+dbSNP_file+' -p '+DIR+'/dbSNP/'+sampleID+'\n')

	output1.write('\n')
	output1.write('sleep 30\n')
	output1.write('exit 0')

	output1.close()

	os.system('qsub PBS_bcfisec_'+sampleID+'.pbs')


def cleanup():
    for sampleID in sampleList:
	varDict={}
	input1=open(DIR+'/dbSNP/'+sampleID+'/0003.vcf','r')
	for line in input1:
	    if not line.count('#') > 0:
		each=line.strip().split('\t')
		chrID='chr'+each[0]
		pt=each[1]
		rsID=each[2]
		ref=each[3]
		alt=each[4]

		keyID=chrID+';'+pt+';'+ref+';'+alt
		varDict[keyID]=rsID
	input1.close()

	input1=open(DIR+'/dbSNP/'+sampleID+'/0002.vcf','r')
	output1=open(DIR+'/annovar/'+sampleID+'.dbSNP.final.vcf','w')
	for line in input1:
	    if not line.count('#') > 0:
		each=line.strip().split('\t')
		chrID='chr'+each[0]
		pt=each[1]
		ref=each[3]
		alt=each[4]
	    
		keyID=chrID+';'+pt+';'+ref+';'+alt
		each[2]=varDict[keyID]
		output1.write('\t'.join(each)+'\n')

	    else:
		output1.write(line)
	input1.close()
	output1.close()

    ## remove dbSNP dir


def ANNOVAR():
    ## annotate by ANNOVAR
    print 'performing ANNOVAR'
    for sampleID in sampleList:
	print sampleID

	run_line='table_annovar.pl '+DIR+'/annovar/'+sampleID+'.dbSNP.final.vcf '+code_DIR+'/humandb/ -buildver hg19 -out '+DIR+'/annovar/'+sampleID+' -remove -protocol refGene -operation g -nastring . -vcfinput -polish'
	print run_line
	os.system(run_line)

	run_line='annotate_variation.pl -geneanno -dbtype refGene -buildver hg19 '+DIR+'/annovar/'+sampleID+'.avinput '+code_DIR+'/humandb/'
	print run_line
	os.system(run_line)


def extract_exon():


    for sampleID in sampleList:
	print sampleID

	input1=open(DIR+'/annovar/'+sampleID+'.hg19_multianno.txt','r')
	all_input1=input1.readlines()
	output1=open(DIR+'/annovar/'+sampleID+'.hg19_multianno.exonic.txt','w')
	output1.write('\t'.join(all_input1[0].strip().split('\t'))+'\n')
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    if each[5]=='exonic':
		output1.write('\t'.join(each)+'\n')

	input1.close()
	output1.close()





#processVCF()
#match_dbSNP()
cleanup()
ANNOVAR()
extract_exon()



