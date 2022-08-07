#!/home/ajl1213/anaconda2/bin/python


import os


DIR=os.getcwd()
inputDIR=DIR+'/Fastq/'
ncpu=12
ref='/home/ambrosio/reference/refdata-cellranger-atac-hg19-1.2.0/'
sampleList=[
'JW568',
'JW569',
'JW570'
]


def make_script(sampleID):
    output1=open('PBS_scATAC_'+sampleID+'.pbs','w')
    output1.write('#PBS -N scATAC.'+sampleID+'\n')
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

    output1.write('cellranger-atac count --id='+sampleID+' --reference='+ref+' --fastqs='+inputDIR+'/'+sampleID+'/ --sample='+sampleID+'\n')
    
    output1.write('\n')
    output1.write('sleep 30\n')
    output1.write('exit 0')

    output1.close()



for sampleID in sampleList:
    make_script(sampleID)
    os.system('qsub PBS_scATAC_'+sampleID+'.pbs')




