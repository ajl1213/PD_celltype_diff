#!/home/ajl1213/anaconda2/bin/python


import os


ncpu=8
DIR=os.getcwd()
inputDIR=DIR+'/SampleFastq/'

FASTA='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa'    # hg19
FAI='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa.fai'  # hg19
PICARDJAR='/home/ajl1213/programs/picard_tool/picard.jar'
sampleList=[
'X4870NOSN',
'X4689NOSN',
'X5628NOSN',
'X5130NOSN',
'X5006NOSN',
'X4996NOSN',
'X5244PDSN',
'X4831PDSN',
'X5215PDSN',
'X5649PDSN',
'X5591PDSN'
]

os.system('mkdir '+DIR+'/FlagStat')
os.system('mkdir '+DIR+'/SortedMergedBam')
os.system('mkdir '+DIR+'/SortedNodupBam')
os.system('mkdir '+DIR+'/PicardStat')
os.system('mkdir '+DIR+'/TransCis')


def makeScript():
    for sampleID in sampleList:
	output1=open('PBS_hicmapping_'+sampleID+'.pbs','w')
	output1.write('#PBS -N hic_mapping_'+sampleID+'\n')
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
	output1.write('\n')
	output1.write('zcat '+inputDIR+'/'+sampleID+'_1.fq.gz | bwa mem -M -t '+str(ncpu)+' '+FASTA+' - > '+sampleID+'_1.sam\n')
	output1.write('zcat '+inputDIR+'/'+sampleID+'_2.fq.gz | bwa mem -M -t '+str(ncpu)+' '+FASTA+' - > '+sampleID+'_2.sam\n')
	output1.write('samtools flagstat '+sampleID+'_1.sam > FlagStat/'+sampleID+'_1.flagstatResult.txt\n')
	output1.write('samtools flagstat '+sampleID+'_2.sam > FlagStat/'+sampleID+'_2.flagstatResult.txt\n')
	output1.write('python 2b_make_paired_sam.py '+sampleID+'_1.sam '+sampleID+'_2.sam | samtools view -Sb -t '+FAI+' -o SortedMergedBam/'+sampleID+'.merged.bam\n')
	output1.write('samtools sort  -@ 8 -m 4G SortedMergedBam/'+sampleID+'.merged.bam > SortedMergedBam/'+sampleID+'.sorted.merged.bam\n')
	output1.write('java -Xmx50g -jar '+PICARDJAR+' MarkDuplicates VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE I=SortedMergedBam/'+sampleID+'.sorted.merged.bam O=SortedNodupBam/'+sampleID+'.sorted.nodup.bam M=PicardStat/'+sampleID+'.METRIX.TXT TMP_DIR=/tmp\n')
	output1.write('samtools view -S SortedNodupBam/'+sampleID+'.sorted.nodup.bam | python check_trans_cis.py TransCis/'+sampleID+'.TransCis.txt\n')
	output1.write('samtools index SortedNodupBam/'+sampleID+'.sorted.nodup.bam\n')

	output1.write('\n')
	output1.write('sleep 30\n')
	output1.write('exit 0')
	output1.close()


	os.system('qsub PBS_hicmapping_'+sampleID+'.pbs')



makeScript()




