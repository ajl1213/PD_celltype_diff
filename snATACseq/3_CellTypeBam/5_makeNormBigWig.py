#!/home/ajl1213/anaconda2/bin/python


import os
import numpy


ncpu=8
DIR=os.getcwd()
FASTA='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa'    
FAI='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa.fai'  
sampleList=['AllCells','DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']
dataTypeList=['NOSN','PDSN']


os.system('mkdir '+DIR+'/Bigwig')


def getTotalCount():
    totalCountDict={}
    input1=open(DIR+'/FRIP/FRIP.final.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	sampleID=each[0]
	dataType=each[1]
	peakCount=each[3]
	totalCountDict[sampleID+'.'+dataType]=peakCount
    input1.close()

    return totalCountDict


def makeScript():
    totalCountDict=getTotalCount()

    for sampleID in sampleList:
	for dataType in dataTypeList:
	    totalCount=totalCountDict[sampleID+'.'+dataType]

	    output1=open('PBS_ChipBigWig_'+sampleID+'.'+dataType+'.pbs','w')
	    output1.write('#PBS -N ChipBigWig_'+sampleID+'.'+dataType+'\n')
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
	    output1.write('genomeCoverageBed -bg -ibam '+DIR+'/MergedBam/'+sampleID+'.'+dataType+'.bam -g '+FAI+' -scale '+str((float(1000000)/float(totalCount)))+' > '+DIR+'/Bigwig/'+sampleID+'.'+dataType+'.bg\n')
	    output1.write('cat '+DIR+'/Bigwig/'+sampleID+'.'+dataType+'.bg | grep \'chr\' | sort -k1,1 -k2,2n - > '+DIR+'/Bigwig/'+sampleID+'.'+dataType+'.sorted.bg\n')
	    output1.write('bedGraphToBigWig '+DIR+'/Bigwig/'+sampleID+'.'+dataType+'.sorted.bg '+FAI+' '+DIR+'/Bigwig/'+sampleID+'.'+dataType+'.bw\n')
	    output1.write('\n')
	    output1.write('sleep 30\n')
	    output1.write('exit 0')
	    output1.close()

	    os.system('qsub PBS_ChipBigWig_'+sampleID+'.'+dataType+'.pbs')



makeScript()


