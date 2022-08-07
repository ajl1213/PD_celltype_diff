#!/usr/bin/python


import os



DIR=os.getcwd()
inputDIR=DIR+'/Fastq'
ncpu=8
GTF='/home/ajl1213/genome.info/gencode/hg19.release38/GTF/gencode.v38lift37.annotation.gtf'
StarGenome_DIR='/home/ajl1213/genome.info/StarGenome/hg19'
RsemGenome_DIR='/home/ajl1213/genome.info/RsemRef/hg19/hg19'

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


print len(sampleList)
os.system('mkdir '+DIR+'/AlignedBam')
os.system('mkdir '+DIR+'/CountInfo')


def makeScript():
    output1=open('PBS_RNA_'+sampleID+'.pbs','w')
    output1.write('#PBS -N RNA_'+sampleID+'\n')
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

    output1.write('STAR --runThreadN '+str(ncpu)+' --genomeDir '+StarGenome_DIR+' --sjdbGTFfile '+GTF+' --readFilesIn '+inputDIR+'/'+sampleID+'_1.fq.gz '+inputDIR+'/'+sampleID+'_2.fq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix '+DIR+'/AlignedBam/'+sampleID+'/ --outWigType wiggle --quantMode TranscriptomeSAM\n')
    output1.write('rsem-calculate-expression --bam -p '+str(ncpu)+' --paired-end --strandedness reverse '+DIR+'/AlignedBam/'+sampleID+'/Aligned.toTranscriptome.out.bam '+RsemGenome_DIR+' '+DIR+'/CountInfo/'+sampleID+'\n')

    output1.write('\n')
    output1.write('sleep 30\n')
    output1.write('exit 0')
    output1.close()


for sampleID in sampleList:
    makeScript()
    os.system('qsub PBS_RNA_'+sampleID+'.pbs')




