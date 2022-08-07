#!/home/ajl1213/anaconda2/bin/python


import os



ncpu=4
DIR=os.getcwd()
inputDIR='/home/ajl1213/Projects/PD/SciAdv_Data/snATACseq/1_Process'
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

os.system('mkdir '+DIR+'/BamFiles')
os.system('mkdir '+DIR+'/AlignStat')


def ExtractCellScript(sampleID, cellType):
    output1=open('PBS_CellExtract_'+sampleID+'.'+cellType+'.pbs','w')
    output1.write('#PBS -N CellExtract_'+sampleID+'.'+cellType+'\n')
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

    output1.write('python 2b_extractCells.py '+sampleID+' '+cellType+' '+DIR+'/Barcodes/'+sampleID+'.'+cellType+'.barcodes.txt '+inputDIR+'/'+sampleID+'/outs/possorted_bam.bam\n')

    output1.write('\n')
    output1.write('sleep 30\n')
    output1.write('exit 0')

    output1.close()


for sampleID in sampleList:
    for cellType in cellTypeList:
        ExtractCellScript(sampleID, cellType)
	os.system('qsub PBS_CellExtract_'+sampleID+'.'+cellType+'.pbs')





