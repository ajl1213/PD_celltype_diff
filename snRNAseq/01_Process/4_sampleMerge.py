#!/home/ajl1213/anaconda2/bin/python


import os


DIR=os.getcwd()
ncpu=16
ref='/home/ambrosio/reference/hg19-3.0.0_premrna/'
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

outputName='MergedAll'


def makeMergeInfo():
    output1=open(DIR+'/aggr.'+outputName+'.csv', 'w')
    output1.write('library_id,molecule_h5\n')
    for sampleID in sampleList:
	output1.write(sampleID+',')
	output1.write(DIR+'/'+sampleID+'/outs/molecule_info.h5\n')
    output1.close()


def makeScript():
    output1=open('PBS_scRNA_aggr.'+outputName+'.pbs','w')
    output1.write('#PBS -N scRNA_aggr.'+outputName+'\n')
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

    output1.write('cellranger aggr --id='+outputName+' --csv='+DIR+'/aggr.'+outputName+'.csv --normalize=mapped\n')
    
    output1.write('\n')
    output1.write('sleep 30\n')
    output1.write('exit 0')

    output1.close()



makeMergeInfo()
makeScript()
os.system('qsub PBS_scRNA_aggr.'+outputName+'.pbs')



