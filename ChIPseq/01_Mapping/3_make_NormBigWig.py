#!/usr/bin/python


import os
import numpy


DIR=os.getcwd()
FASTA='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa'    # hg19
FAI='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa.fai'  # hg19
threads=4

sampleList=[
'X4870_1NOSN;H3K27ac',
'X4870_2NOSN;H3K27ac',
'X4689NOSN;H3K27ac',  
'X5628_1NOSN;H3K27ac',
'X5628_2NOSN;H3K27ac',
'X5130NOSN;H3K27ac',
'X4996_1NOSN;H3K27ac',
'X4996_2NOSN;H3K27ac',
'X5006NOSN;H3K27ac',
'X5244PDSN;H3K27ac',
'X4831PDSN;H3K27ac',
'X5532PDSN;H3K27ac',
'X5215PDSN;H3K27ac',
'X5627PDSN;H3K27ac',
'X5778PDSN;H3K27ac',
'X5742PDSN;H3K27ac',
'X5649PDSN;H3K27ac',
'X5591PDSN;H3K27ac'
]

os.system('mkdir '+DIR+'/Bigwig')



def getStat():
    PeakCountDict={}
    for dataType in ['H3K27ac']:
	input1=open(DIR+'/FRIP/'+dataType+'.FRIP.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    [sampleID, peakNumber, peakCount, totalCount, FRIP]=line.strip().split('\t')
	    keyID=sampleID+';'+dataType
	    PeakCountDict[keyID]=float(peakCount)	    
	input1.close()

    return PeakCountDict


def makeScript():
    PeakCountDict = getStat()

    for i in sampleList:
	sampleID=i.split(';')[0]
	dataType=i.split(';')[1]

	output1=open('PBS_ChipBigWig_'+sampleID+'.'+dataType+'.pbs','w')
	output1.write('#PBS -N ChipBigWig_'+sampleID+'_'+dataType+'\n')
	output1.write('#PBS -q workq\n')
	output1.write('#PBS -l nodes=1:ppn='+str(threads)+'\n')
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
	### Total peak read normalized
	output1.write('genomeCoverageBed -bg -ibam '+DIR+'/SortedNodupFilteredBam/'+sampleID+'.'+dataType+'.sorted.nodup.filtered.bam -g '+FAI+' -scale '+str((float(1000000)/float(PeakCountDict[i])))+' > '+DIR+'/Bigwig/'+sampleID+'.'+dataType+'.bg\n')
	output1.write('sort -k1,1 -k2,2n '+DIR+'/Bigwig/'+sampleID+'.'+dataType+'.bg > '+DIR+'/Bigwig/'+sampleID+'.'+dataType+'.sorted.bg\n')
	output1.write('bedGraphToBigWig '+DIR+'/Bigwig/'+sampleID+'.'+dataType+'.sorted.bg '+FAI+' '+DIR+'/Bigwig/'+sampleID+'.'+dataType+'.bw\n')
	output1.write('\n')
	output1.write('sleep 30\n')
	output1.write('exit 0')
	output1.close()

	os.system('qsub PBS_ChipBigWig_'+sampleID+'.'+dataType+'.pbs')


makeScript()



