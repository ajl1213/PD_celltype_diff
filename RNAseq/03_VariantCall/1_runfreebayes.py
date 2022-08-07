#!/home/ajl1213/anaconda2/bin/python
## conducted in ajl2 account


import os


DIR=os.getcwd()
cpu=4
BamDIR='/home/ajl1213/Projects/PD/SciAdv_Data/RNA/01_Mapping/AlignedBam'
FASTA='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa'

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

os.system('mkdir '+DIR+'/freebayes')


def bringBam():
	for sampleID in sampleList:
		runLine='cp '+BamDIR+'/'+sampleID+'/Aligned.sortedByCoord.out.bam '+DIR+'/freebayes/'+sampleID+'.RNA.bam'
		print(runLine)
		os.system(runLine)

def makeScript():

	for sampleID in sampleList:

		output1=open('PBS_freebayes.'+sampleID+'.pbs','w')
		output1.write('#PBS -N freebayes.'+sampleID+'\n')
		output1.write('#PBS -q workq\n')
		output1.write('#PBS -l nodes=1:ppn='+str(cpu)+'\n')
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

		#
		output1.write('source activate souporcell\n')
		output1.write('freebayes -f '+FASTA+' -iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 20 --pooled-continuous --skip-coverage 100000 '+DIR+'/freebayes/'+sampleID+'.RNA.bam > '+DIR+'/freebayes/'+sampleID+'.vcf\n')

		output1.write('\n')
		output1.write('sleep 30\n')
		output1.write('exit 0')

		output1.close()

		os.system('qsub PBS_freebayes.'+sampleID+'.pbs')



bringBam()
makeScript()


