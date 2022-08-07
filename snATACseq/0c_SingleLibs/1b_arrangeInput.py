#!/usr/bin/python

import os


ncpu=4
DIR=os.getcwd()

sampleList=[
'/home/CBFG/20210803_MGIseq;SW228;X5006NOSN',
'/home/CBFG/20210803_MGIseq;SW227;X5215PDSN',
'/home/CBFG/20210803_MGIseq;SW229;X4831PDSN'
]

os.system('mkdir '+DIR+'/Fastq')


def makeScript():
    for i in sampleList:
	inputDIR=i.split(';')[0]
	seqID=i.split(';')[1]
	sampleID=i.split(';')[2]
	
	os.system('mkdir '+DIR+'/Fastq/'+seqID)

	## R1
	output1=open('PBS_'+seqID+'.R1.pbs','w')
	output1.write('#PBS -N '+seqID+'.R1\n')
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
	output1.write('zcat '+inputDIR+'/Sample_'+seqID+'/'+seqID+'_1.fq.gz | awk -F \'\' \'{if($1=="@"&& $2=="V"){gsub("/1$",":1:X:1:3:4:1 3:N:0:0",$0); print $0} else{print $0}}\' | gzip --stdout --fast > '+DIR+'/Fastq/'+sampleID+'/'+sampleID+'_S1_L001_R1_001.fastq.gz')
	output1.write('\n')
	output1.write('sleep 30\n')
	output1.write('exit 0')
	output1.close()

	os.system('qsub PBS_'+seqID+'.R1.pbs')

	## R2
	output1=open('PBS_'+seqID+'.R2.pbs','w')
	output1.write('#PBS -N '+seqID+'.R2\n')
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
	output1.write('zcat '+inputDIR+'/Sample_'+seqID+'/'+seqID+'_3.fq.gz | awk -F \'\' \'{if($1=="@"&& $2=="V"){gsub("/2$",":1:X:1:3:4:1 3:N:0:0",$0); print $0} else{print $0}}\' | gzip --stdout --fast > '+DIR+'/Fastq/'+sampleID+'/'+sampleID+'_S1_L001_R2_001.fastq.gz')
	output1.write('\n')
	output1.write('sleep 30\n')
	output1.write('exit 0')
	output1.close()

	os.system('qsub PBS_'+seqID+'.R2.pbs')

	## R3
	output1=open('PBS_'+seqID+'.R3.pbs','w')
	output1.write('#PBS -N '+seqID+'.R3\n')
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
	output1.write('zcat '+inputDIR+'/Sample_'+seqID+'/'+seqID+'_2.fq.gz | awk -F \'\' \'{if($1=="@"&& $2=="V"){gsub("/2$",":1:X:1:3:4:1 3:N:0:0",$0); print $0} else{print $0}}\' | gzip --stdout --fast > '+DIR+'/Fastq/'+sampleID+'/'+sampleID+'_S1_L001_R3_001.fastq.gz')
	output1.write('\n')
	output1.write('sleep 30\n')
	output1.write('exit 0')
	output1.close()

	os.system('qsub PBS_'+seqID+'.R3.pbs')



makeScript()



