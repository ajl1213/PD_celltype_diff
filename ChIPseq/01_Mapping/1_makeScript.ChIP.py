#!/home/ajl1213/anaconda2/bin/python


import os


DIR=os.getcwd()
inputDIR=DIR+'/Fastq'
PICARDJAR='/home/ajl1213/programs/picard_tool/picard.jar'
FASTA='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa'    
FAI='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa.fai'  
jobTitle='ChIP'
ncpu=8

sampleList=[
'X4870_1NOSN.Input','X4870_1NOSN.H3K27ac',
'X4870_2NOSN.Input','X4870_2NOSN.H3K27ac',
'X4689NOSN.Input',  'X4689NOSN.H3K27ac',  
'X5628_1NOSN.Input','X5628_1NOSN.H3K27ac',
'X5628_2NOSN.Input','X5628_2NOSN.H3K27ac',
'X5130NOSN.Input',  'X5130NOSN.H3K27ac',
'X4996_1NOSN.Input','X4996_1NOSN.H3K27ac',
'X4996_2NOSN.Input','X4996_2NOSN.H3K27ac',
'X5006NOSN.Input',  'X5006NOSN.H3K27ac',
'X5244PDSN.Input',  'X5244PDSN.H3K27ac',
'X4831PDSN.Input',  'X4831PDSN.H3K27ac',
'X5532PDSN.Input',  'X5532PDSN.H3K27ac',
'X5215PDSN.Input',  'X5215PDSN.H3K27ac',
'X5627PDSN.Input',  'X5627PDSN.H3K27ac',
'X5778PDSN.Input',  'X5778PDSN.H3K27ac',
'X5742PDSN.Input',  'X5742PDSN.H3K27ac',
'X5649PDSN.Input',  'X5649PDSN.H3K27ac',
'X5591PDSN.Input',  'X5591PDSN.H3K27ac'
]


print len(sampleList)

os.system('mkdir '+DIR+'/SortedBam')
os.system('mkdir '+DIR+'/SortedNodupBam')
os.system('mkdir '+DIR+'/SortedNodupFilteredBam')
os.system('mkdir '+DIR+'/Flagstat')
os.system('mkdir '+DIR+'/PicardStat')


def makeScript(sampleID):
    output1=open('PBS_'+jobTitle+'.'+sampleID+'.pbs','w')
    output1.write('#PBS -N '+jobTitle+'.'+sampleID+'\n')
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
    output1.write('bwa mem -M -t '+str(ncpu)+' '+FASTA+' '+inputDIR+'/'+sampleID+'.R1.fq.gz '+inputDIR+'/'+sampleID+'.R2.fq.gz | samtools view -bS -t '+FAI+' -o '+DIR+'/SortedBam/'+sampleID+'.bam\n')
    output1.write('samtools flagstat '+DIR+'/SortedBam/'+sampleID+'.bam > '+DIR+'/Flagstat/'+sampleID+'.flagstat.txt\n')
    output1.write('samtools sort -m 4G -@ '+str(ncpu)+' '+DIR+'/SortedBam/'+sampleID+'.bam > '+DIR+'/SortedBam/'+sampleID+'.sorted.bam\n')
    output1.write('java -Xmx50g -jar '+PICARDJAR+' MarkDuplicates VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE I='+DIR+'/SortedBam/'+sampleID+'.sorted.bam O='+DIR+'/SortedNodupBam/'+sampleID+'.sorted.nodup.bam M='+DIR+'/PicardStat/'+sampleID+'.METRIX.TXT TMP_DIR=/tmp\n')
    output1.write('samtools view -bS -q 10 -t '+FAI+' '+DIR+'/SortedNodupBam/'+sampleID+'.sorted.nodup.bam > '+DIR+'/SortedNodupFilteredBam/'+sampleID+'.sorted.nodup.filtered.bam\n')
    output1.write('samtools flagstat '+DIR+'/SortedNodupFilteredBam/'+sampleID+'.sorted.nodup.filtered.bam > '+DIR+'/Flagstat/'+sampleID+'.flagstat.final.txt\n')
    output1.write('\n')
    output1.write('sleep 30\n')
    output1.write('exit 0')
    output1.close()



for sampleID in sampleList:
    makeScript(sampleID)
    os.system('qsub PBS_'+jobTitle+'.'+sampleID+'.pbs')






