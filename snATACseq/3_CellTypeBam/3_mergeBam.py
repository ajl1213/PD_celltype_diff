#!/home/ajl1213/anaconda2/bin/python


import os


ncpu=4
DIR=os.getcwd()

sampleList=[
'X4870NOSN',
'X4996NOSN',
'X5006NOSN',
'X5130NOSN',
'X5628NOSN',
'Public1NOSN',
'Public2NOSN',
'X4831PDSN',
'X5215PDSN',
'X5244PDSN',
'X5532PDSN',
'X5591PDSN',
'X5627PDSN',
'X5742PDSN',
'X5778PDSN'
]

cellTypeList=['AllCells','DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']
dataTypeList=['NOSN','PDSN']


os.system('mkdir '+DIR+'/MergedBam')
os.system('mkdir '+DIR+'/flagstat')


def make_script():

    ##
#    for cellType in cellTypeList:
#	output1=open('PBS_ChIPMerge_'+cellType+'.pbs','w')
#	output1.write('#PBS -N ChIPMerge_'+cellType+'\n')
#	output1.write('#PBS -q workq\n')
#	output1.write('#PBS -l nodes=1:ppn='+str(ncpu)+'\n')
#	output1.write('#PBS -j oe\n')
#	output1.write('\n')
#	output1.write('# go workdir\n')
#	output1.write('cd $PBS_O_WORKDIR\n')
#	output1.write('\n')
#	output1.write('# run command \n')
#	output1.write('sleep 5\n')
#	output1.write('\n')
#	output1.write('echo -n \"I am on: \"\n')
#	output1.write('hostname;\n')    
#	output1.write('echo finding ssh-accessible nodes:\n')
#	output1.write('echo -n \"running on: \"\n')
#	output1.write('\n')
#	output1.write('samtools merge '+DIR+'/MergedBam/'+cellType+'.bam')
#	for sampleID in sampleList:
#	    output1.write(' '+DIR+'/BamFiles/'+cellType+'.'+sampleID+'.filter.bam')
#	output1.write('\n')
#	output1.write('samtools flagstat '+DIR+'/MergedBam/'+cellType+'.bam > '+DIR+'/flagstat/'+cellType+'.flagstat.txt\n')
#
#	output1.write('\n')
#	output1.write('sleep 30\n')
#	output1.write('exit 0')
#	output1.close()
#
#	make_script()
#	os.system('qsub PBS_ChIPMerge_'+cellType+'.pbs')


    ##
    for cellType in cellTypeList:
	for dataType in dataTypeList:
	    print cellType, dataType

	    output1=open('PBS_ChIPMerge_'+cellType+'.'+dataType+'.pbs','w')
	    output1.write('#PBS -N ChIPMerge_'+cellType+'.'+dataType+'\n')
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
	    output1.write('samtools merge '+DIR+'/MergedBam/'+cellType+'.'+dataType+'.bam')
	    for sampleID in sampleList:
		if sampleID.count(dataType) > 0:
		    output1.write(' '+DIR+'/BamFiles/'+cellType+'.'+sampleID+'.filter.bam')
	    output1.write('\n')
	    output1.write('samtools flagstat '+DIR+'/MergedBam/'+cellType+'.'+dataType+'.bam > '+DIR+'/flagstat/'+cellType+'.'+dataType+'.flagstat.txt\n')

	    output1.write('\n')
	    output1.write('sleep 30\n')
	    output1.write('exit 0')
	    output1.close()

	    os.system('qsub PBS_ChIPMerge_'+cellType+'.'+dataType+'.pbs')



make_script()



