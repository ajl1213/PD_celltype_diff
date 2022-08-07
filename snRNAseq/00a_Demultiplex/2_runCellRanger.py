#!/home/ajl1213/anaconda2/bin/python


import os
import glob


DIR=os.getcwd()
inputDIR=DIR+'/Fastq/'
ref='/home/ambrosio/reference/hg19-3.0.0_premrna/'
ncpu=16

sampleList=['SW211', 'SW212', 'BG40BG41', 'BG42BG43']



def makeScript():

    for sampleID in sampleList:
        output1=open('PBS_scRNA_'+sampleID+'.pbs','w')
        output1.write('#PBS -N scRNA.'+sampleID+'\n')
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

	output1.write('cellranger count --id='+sampleID+' --transcriptome='+ref+' --fastqs='+inputDIR+'/'+sampleID+' --sample='+sampleID+'\n')

        output1.write('\n')
        output1.write('sleep 30\n')
        output1.write('exit 0')

        output1.close()

        os.system('qsub PBS_scRNA_'+sampleID+'.pbs')



makeScript()



