#!/home/ajl1213/anaconda2/bin/python


import os
import glob


DIR=os.getcwd()
inputDIR=DIR+'/Fastq'
cpu=16
ref='/home/ambrosio/reference/hg19-3.0.0_premrna/'
cellCountDIR='/home/ajl1213/Projects/PD/SciAdv_Data/snRNAseq/00a_Demultiplex/Barcodes'

sampleList=[
'X5628NOSN',
'X4689NOSN',
'X4870NOSN',
'X4996NOSN',
'X5006NOSN',
'X5130NOSN',

'Public1NOSN',
'Public2_1NOSN',
'Public2_2NOSN',
'Public3_1NOSN',
'Public3_2NOSN',
'Public4NOSN',
'Public5NOSN',

'X5591PDSN',
'X4831PDSN',
'X5215PDSN',
'X5244PDSN',
'X5742PDSN',
'X5778PDSN'
]


def loadCellCount():
    cellCountDict={}
    for i in sorted(glob.glob(cellCountDIR+'/*.barcodes.txt')):
        cellCount=0
        sampleID=i.split('/')[-1].split('.')[0]
        input1=open(cellCountDIR+'/'+sampleID+'.barcodes.txt','r')
        all_input1=input1.readlines()
        for line in all_input1:
            cellCount+=1
        input1.close()

        cellCountDict[sampleID]=cellCount

    return cellCountDict


def makeScript():
    cellCountDict=loadCellCount()

    for sampleID in sampleList:
	output1=open('PBS_scRNA_'+sampleID+'.pbs','w')
	output1.write('#PBS -N scRNA.'+sampleID+'\n')
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
	
	if cellCountDict.has_key(sampleID):
	    output1.write('cellranger count --force-cells='+str(cellCountDict[sampleID])+' --id='+sampleID+' --transcriptome='+ref+' --fastqs='+inputDIR+'/'+sampleID+' --sample='+sampleID+'\n')
	else:
	    output1.write('cellranger count --id='+sampleID+' --transcriptome='+ref+' --fastqs='+inputDIR+'/'+sampleID+' --sample='+sampleID+'\n')


	output1.write('\n')
	output1.write('sleep 30\n')
	output1.write('exit 0')

	output1.close()

	os.system('qsub PBS_scRNA_'+sampleID+'.pbs')



makeScript()



