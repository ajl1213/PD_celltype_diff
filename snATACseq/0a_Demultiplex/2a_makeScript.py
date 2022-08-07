#!/home/ajl1213/anaconda2/bin/python


import os
import gzip


## After souporcell and cluster-patient assignment

DIR=os.getcwd()
ncpu=8
barcodeDIR='/home/ajl2/Projects/PD/Data/snATACseq/0a_Demultiplex.v2'
inFastqDIR='/home/ajl1213/Projects/PD/SciAdv_Data/snATACseq/0a_Demultiplex/Fastq'
outFastqDIR=DIR+'/FastqDemul'
readTypeList=['R1','R2','R3']
sampleDict={
'JW568_0':'X5627PDSN',
'JW568_1':'X5532PDSN',

'JW569_0':'X4996NOSN',
'JW569_1':'X5130NOSN',

'JW570_0':'X5591PDSN',
'JW570_1':'X5244PDSN'
}


laneDict={
'JW568':1,
'JW569':1,
'JW570':1,
}

os.system('mkdir '+DIR+'/Barcodes')
os.system('mkdir '+DIR+'/ReadPerCell')
os.system('mkdir '+DIR+'/FastqDemul')



def getBarcodes():
    for i in sampleDict:
	seqID=i.split('_')[0]
	pseudoID=i.split('_')[1]
	sampleID=sampleDict[i]

	input1=open(barcodeDIR+'/souporcell_'+seqID+'/clusters.tsv','r')
	output1=open(DIR+'/Barcodes/'+sampleID+'.barcodes.txt','w')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    barcodeID=each[0]
	    if each[2] == pseudoID: 
		output1.write(barcodeID+'\n')
	input1.close()
	output1.close()


def makeScript():
    for i in sampleDict:
	seqID=i.split('_')[0]
	sampleID=sampleDict[i]
	os.system('mkdir '+outFastqDIR+'/'+sampleID)

	for nLane in range(laneDict[seqID]):
	    nLane=nLane+1

	    for readType in readTypeList:

		output1=open('PBS_splitFastq.'+sampleID+'.lane'+str(nLane)+'.'+readType+'.pbs','w')
		output1.write('#PBS -N splitFastq.'+sampleID+'.lane'+str(nLane)+'.'+readType+'\n')
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

		output1.write('python 2b_splitFastq.py '+seqID+' '+sampleID+' '+str(nLane)+' '+readType+' '+DIR+'/Barcodes/'+sampleID+'.barcodes.txt '+DIR+'/'+seqID+'/outs/possorted_bam.bam '+inFastqDIR+'/'+seqID+' '+outFastqDIR+'/'+sampleID+'\n')

		output1.write('\n')
		output1.write('sleep 30\n')
		output1.write('exit 0')

		output1.close()

		os.system('qsub PBS_splitFastq.'+sampleID+'.lane'+str(nLane)+'.'+readType+'.pbs')



getBarcodes()
makeScript()



