#!/usr/bin/python

import os


DIR=os.getcwd()
FAI='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa.fai'
BamDIR=DIR+'/SortedNodupFilteredBam'
MergedPeakBed='/home/ajl1213/Projects/PD/data/snATAC/03_CellTypeBam/PeakList/snATAC.peaks.bed'

##
histMark='H3K27ac'
sampleList=[
'X4870_1NOSN-X4870_2NOSN-X4689NOSN-X5628_1NOSN-X5628_2NOSN-X5130NOSN-X4996_1NOSN-X4996_2NOSN-X5006NOSN;MergedNOSN',
'X5244PDSN-X4831PDSN-X5532PDSN-X5215PDSN-X5627PDSN-X5778PDSN-X5742PDSN-X5649PDSN-X5591PDSN;MergedPDSN'
]


def merge():
    ncpu=8
    for i in sampleList:
	mergeList=i.split(';')[0].split('-')
	newLabel=i.split(';')[1]
	output1=open('PBS_ChIPMerge_'+newLabel+'.'+histMark+'.pbs','w')
	output1.write('#PBS -N ChIPMerge_'+newLabel+'.'+histMark+'\n')
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
	output1.write('samtools merge '+BamDIR+'/'+newLabel+'.'+histMark+'.sorted.nodup.filtered.bam')
	for j in mergeList:
	    output1.write(' '+BamDIR+'/'+j+'.'+histMark+'.sorted.nodup.filtered.bam')
	output1.write('\n')
	output1.write('samtools flagstat '+BamDIR+'/'+newLabel+'.'+histMark+'.sorted.nodup.filtered.bam > '+DIR+'/Flagstat/'+newLabel+'.'+histMark+'.flagstat.final.txt\n')

	output1.write('\n')
	output1.write('sleep 30\n')
	output1.write('exit 0')
	output1.close()

	os.system('qsub PBS_ChIPMerge_'+newLabel+'.'+histMark+'.pbs')


def peakCoverage():
    ncpu=4
    for i in sampleList:
	newLabel=i.split(';')[1]
        output1=open(DIR+'/PBS_peakCoverage.'+newLabel+'.'+histMark+'.pbs','w')
        output1.write('#PBS -N peakCoverage.'+newLabel+'.'+histMark+'\n')
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
        output1.write('coverageBed -counts -abam '+BamDIR+'/'+newLabel+'.'+histMark+'.sorted.nodup.filtered.bam -b '+MergedPeakBed+' > '+DIR+'/PeakCoverage/'+newLabel+'.'+histMark+'.peakCoverage.txt\n')
        output1.write('\n')
        output1.write('sleep 30\n')
        output1.write('exit 0')
        output1.close()

        os.system('qsub '+DIR+'/PBS_peakCoverage.'+newLabel+'.'+histMark+'.pbs')


def getTotalCount():
    totalCountDict={}

    for i in sampleList:
	newLabel=i.split(';')[1]

        input1=open(DIR+'/Flagstat/'+newLabel+'.'+histMark+'.flagstat.final.txt','r')
        all_input1=input1.readlines()
        for line in all_input1:
            if line.count('mapped (')==1:
                each=line.strip().split(' ')
                totalCount=each[0]

                totalCountDict[newLabel]=totalCount

        input1.close()

    return totalCountDict


def calcFRIP():
    totalCountDict=getTotalCount()

    nPeakDict={}
    peakCountDict={}

    ## Load data
    for i in sampleList:
	newLabel=i.split(';')[1]

        nPeak=0
        peakCount=0
        input1=open(DIR+'/PeakCoverage/'+newLabel+'.'+histMark+'.peakCoverage.txt','r')
        all_input1=input1.readlines()
        for line in all_input1:
            each=line.strip().split('\t')
            readCount=int(each[-1])
            peakCount+=readCount    
            nPeak+=1

        input1.close()

        nPeakDict[newLabel]=nPeak
        peakCountDict[newLabel]=peakCount

    ## Compute FRIP    
    output1=open(DIR+'/FRIP/'+histMark+'.FRIP.merged.txt','w')
    output1.write('sampleID\tnPeak\tpeakReads\ttotalReads\tFRIP(%)\n')
    for i in sampleList:
	newLabel=i.split(';')[1]

	nPeak=nPeakDict[newLabel]    
	peakCount=peakCountDict[newLabel]    
	totalCount=totalCountDict[newLabel]    

	fripVal=str(float(peakCount)/float(totalCount)*100)

	outputLine=[newLabel, str(nPeak), str(peakCount), str(totalCount), str(fripVal)]
	output1.write('\t'.join(outputLine)+'\n')
    output1.close()


def getStat():
    PeakCountDict={}
    input1=open(DIR+'/FRIP/'+histMark+'.FRIP.merged.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	[sampleID, peakNumber, peakCount, totalCount, FRIP]=line.strip().split('\t')
	PeakCountDict[sampleID]=float(peakCount)    
    input1.close()

    return PeakCountDict


def makeBigwig():
    PeakCountDict = getStat()
    ncpu=4

    for i in sampleList:
	newLabel=i.split(';')[1]

        output1=open('PBS_ChipBigWig_'+newLabel+'.'+histMark+'.pbs','w')
        output1.write('#PBS -N ChipBigWig_'+newLabel+'_'+histMark+'\n')
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
        ## Total peak read normalized
        output1.write('genomeCoverageBed -bg -ibam '+BamDIR+'/'+newLabel+'.'+histMark+'.sorted.nodup.filtered.bam -g '+FAI+' -scale '+str((float(1000000)/float(PeakCountDict[newLabel])))+' > '+DIR+'/Bigwig/'+newLabel+'.'+histMark+'.bg\n')
        output1.write('sort -k1,1 -k2,2n '+DIR+'/Bigwig/'+newLabel+'.'+histMark+'.bg > '+DIR+'/Bigwig/'+newLabel+'.'+histMark+'.sorted.bg\n')
        output1.write('bedGraphToBigWig '+DIR+'/Bigwig/'+newLabel+'.'+histMark+'.sorted.bg '+FAI+' '+DIR+'/Bigwig/'+newLabel+'.'+histMark+'.bw\n')
        output1.write('\n')
        output1.write('sleep 30\n')
        output1.write('exit 0')
        output1.close()

        os.system('qsub PBS_ChipBigWig_'+newLabel+'.'+histMark+'.pbs')



#merge()
#peakCoverage()
calcFRIP()
makeBigwig()


