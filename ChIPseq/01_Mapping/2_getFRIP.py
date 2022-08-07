#!/home/ajl1213/anaconda2/bin/python


import os


DIR=os.getcwd()
BamDIR=DIR+'/SortedNodupFilteredBam'
MergedPeakBed='/home/ajl1213/Projects/PD/data/snATAC/03_CellTypeBam/PeakList/snATAC.peaks.bed'
sampleList=[
'X4870_1NOSN.H3K27ac',
'X4870_2NOSN.H3K27ac',
'X4689NOSN.H3K27ac',  
'X5628_1NOSN.H3K27ac',
'X5628_2NOSN.H3K27ac',
'X5130NOSN.H3K27ac',
'X4996_1NOSN.H3K27ac',
'X4996_2NOSN.H3K27ac',
'X5006NOSN.H3K27ac',
'X5244PDSN.H3K27ac',
'X4831PDSN.H3K27ac',
'X5532PDSN.H3K27ac',
'X5215PDSN.H3K27ac',
'X5627PDSN.H3K27ac',
'X5778PDSN.H3K27ac',
'X5742PDSN.H3K27ac',
'X5649PDSN.H3K27ac',
'X5591PDSN.H3K27ac'
]


os.system('mkdir '+DIR+'/PeakCoverage')
os.system('mkdir '+DIR+'/FRIP')


def peakCoverage():
    cpu=4
    for i in sampleList:
	sampleID=i.split('.')[0]
	histMark=i.split('.')[1]
	output1=open(DIR+'/PBS_peakCoverage.'+sampleID+'.'+histMark+'.pbs','w')
	output1.write('#PBS -N peakCoverage.'+sampleID+'.'+histMark+'\n')
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
	output1.write('\n')
	output1.write('coverageBed -counts -abam '+BamDIR+'/'+sampleID+'.'+histMark+'.sorted.nodup.filtered.bam -b '+MergedPeakBed+' > '+DIR+'/PeakCoverage/'+sampleID+'.'+histMark+'.peakCoverage.txt\n')
	output1.write('\n')
	output1.write('sleep 30\n')
	output1.write('exit 0')
	output1.close()

	os.system('qsub '+DIR+'/PBS_peakCoverage.'+sampleID+'.'+histMark+'.pbs')



def getTotalCount():
    totalCountDict={}

    for i in sampleList:
	sampleID=i.split('.')[0]
	histMark=i.split('.')[1]

        input1=open(DIR+'/Flagstat/'+sampleID+'.'+histMark+'.flagstat.final.txt','r')
        all_input1=input1.readlines()
        for line in all_input1:
            if line.count('mapped (')==1:
                each=line.strip().split(' ')
                totalCount=each[0]
	    
		keyID=sampleID+'.'+histMark
		totalCountDict[keyID]=totalCount
	
        input1.close()

    return totalCountDict


def calcFRIP():
    totalCountDict=getTotalCount()

    nPeakDict={}
    peakCountDict={}

    ## Load data
    for i in sampleList:
	sampleID=i.split('.')[0]
	histMark=i.split('.')[1]

	nPeak=0
	peakCount=0
	input1=open(DIR+'/PeakCoverage/'+sampleID+'.'+histMark+'.peakCoverage.txt','r')
	all_input1=input1.readlines()
	for line in all_input1:
	    each=line.strip().split('\t')
	    readCount=int(each[-1])
	    peakCount+=readCount	    
	    nPeak+=1

	input1.close()

	keyID=sampleID+'.'+histMark
	nPeakDict[keyID]=nPeak
	peakCountDict[keyID]=peakCount

    ## Compute FRIP    
    for i in ['H3K27ac']:
	output1=open(DIR+'/FRIP/'+i+'.FRIP.txt','w')
	output1.write('sampleID\tnPeak\tpeakReads\ttotalReads\tFRIP(%)\n')
	for j in sampleList:
	    sampleID=j.split('.')[0]
	    histMark=j.split('.')[1]
	    keyID=sampleID+'.'+histMark

	    if i==histMark:
		nPeak=nPeakDict[keyID]		
		peakCount=peakCountDict[keyID]		
		totalCount=totalCountDict[keyID]		

		fripVal=str(float(peakCount)/float(totalCount)*100)

		outputLine=[sampleID, str(nPeak), str(peakCount), str(totalCount), str(fripVal)]
		output1.write('\t'.join(outputLine)+'\n')
	output1.close()


def makeSubjectInfo():
    infoDict={}
    input1=open(DIR+'/Metadata/metadata.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        sampleID=each[0]
        info=each[1:]
        infoDict[sampleID]=info
    input1.close()

    output1=open(DIR+'/Metadata/metadata.final.txt','w')
    output1.write('SampleID\tIndivID\tDiagnosis\tAGE\tGender\tBraakCombo\tMMSE\tDRS\n')
    for i in sampleList:
	i=i.split('.')[0]
        if i.count('_')>0:
            if i.count('NOSN')>0:
                tail='NOSN'
            if i.count('PDSN')>0:
                tail='PDSN'

            sampleID=i.split('_')[0]+tail

        else:
            sampleID=i

        sampleID=sampleID.split(';')[0]
        output1.write(i+'\t'+sampleID+'\t')
        output1.write('\t'.join(infoDict[sampleID])+'\n')
    output1.close()





#peakCoverage()
#calcFRIP()
makeSubjectInfo()




