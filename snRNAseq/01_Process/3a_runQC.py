#!/home/ajl1213/anaconda2/bin/python


import os


DIR=os.getcwd()
ncpu=4
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

os.system('mkdir '+DIR+'/Plots')


def makeScript():

    for sampleID in sampleList:

	output1=open('PBS_scRNA.QC.'+sampleID+'.pbs','w')
	output1.write('#PBS -N scRNA.QC.'+sampleID+'\n')
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

	output1.write('source activate scEnv\n')
	output1.write('R CMD BATCH --no-save --no-restore \'--args '+sampleID+'\' '+DIR+'/3b_sampleQC.R '+DIR+'/Plots/'+sampleID+'.sampleQC.Rout\n')

	output1.write('\n')
	output1.write('sleep 30\n')
	output1.write('exit 0')

	output1.close()

	os.system('qsub PBS_scRNA.QC.'+sampleID+'.pbs')




makeScript()



