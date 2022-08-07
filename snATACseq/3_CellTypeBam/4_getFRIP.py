#!/home/ajl1213/anaconda2/bin/python


import os


DIR=os.getcwd()
BamDIR=DIR+'/MergedBam'
peakFile='/home/ajl1213/Projects/PD/data/snATAC/03_CellTypeBam/PeakList/snATAC.peaks.bed'
sampleList=['AllCells','DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']
dataTypeList=['NOSN','PDSN']

os.system('mkdir '+DIR+'/PeakCoverage/')
os.system('mkdir '+DIR+'/FRIP/')


def peakCoverage():
    ncpu=4
    for sampleID in sampleList:
	for dataType in dataTypeList:
	    
	    output1=open(DIR+'/PBS_peakCoverage_'+sampleID+'.'+dataType+'.pbs','w')
	    output1.write('#PBS -N peakCoverage_'+sampleID+'.'+dataType+'\n')
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
	    #output1.write('coverageBed -counts -bed -a '+peakFile+' -b '+BamDIR+'/'+sampleID+'.'+dataType+'.bam > '+DIR+'/PeakCoverage/'+sampleID+'.'+dataType+'.peakCoverage.txt\n')
	    output1.write('cat '+peakFile+' | cut -f1,2,3 | coverageBed -counts -bed -sorted -a stdin -b '+BamDIR+'/'+sampleID+'.'+dataType+'.bam > '+DIR+'/PeakCoverage/'+sampleID+'.'+dataType+'.peakCoverage.txt\n')
	    output1.write('\n')
	    output1.write('sleep 30\n')
	    output1.write('exit 0')
	    output1.close()

	    #print 'qsub '+DIR+'/PBS_peakCoverage_'+sampleID+'.pbs'
	    os.system('qsub '+DIR+'/PBS_peakCoverage_'+sampleID+'.'+dataType+'.pbs')


def getTotalCount():
    totalCountDict={}

    for sampleID in sampleList:
	for dataType in dataTypeList:
	    input1=open(DIR+'/flagstat/'+sampleID+'.'+dataType+'.flagstat.txt','r')
	    all_input1=input1.readlines()
	    for line in all_input1:
		if line.count('mapped (')==1:
		    each=line.strip().split(' ')
		    totalCount=each[0]
		
		    totalCountDict[sampleID+'.'+dataType]=totalCount
	    
	    input1.close()

    return totalCountDict


def calcFRIP():
    totalCountDict=getTotalCount()

    peakNumberDict={}
    peakCountDict={}

    ## Load data
    for sampleID in sampleList:
	for dataType in dataTypeList:
	    peakNumber=0
	    peakCount=0
	    input1=open(DIR+'/PeakCoverage/'+sampleID+'.'+dataType+'.peakCoverage.txt','r')
	    all_input1=input1.readlines()
	    for line in all_input1:
		each=line.strip().split('\t')
		readCount=int(each[-1])
		peakCount+=readCount	    
		peakNumber+=1

	    input1.close()

	    peakNumberDict[sampleID+'.'+dataType]=peakNumber
	    peakCountDict[sampleID+'.'+dataType]=peakCount

    ## Compute FRIP    
    output1=open(DIR+'/FRIP/FRIP.final.txt','w')
    output1.write('sampleID\tdataType\tpeakNumber\tpeakCount\ttotalCount\tFRIP(%)\n')
    for sampleID in sampleList:
	for dataType in dataTypeList:
	    peakNumber=peakNumberDict[sampleID+'.'+dataType]
	    peakCount=peakCountDict[sampleID+'.'+dataType]		
	    totalCount=totalCountDict[sampleID+'.'+dataType]		

	    fripVal=str(float(peakCount)/float(totalCount)*100)

	    outputLine=[sampleID, dataType, str(peakNumber), str(peakCount), str(totalCount), str(fripVal)]
	    output1.write('\t'.join(outputLine)+'\n')
    output1.close()



#peakCoverage()
calcFRIP()



