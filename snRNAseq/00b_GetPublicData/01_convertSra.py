#!/home/ajl1213/anaconda2/bin/python


import os


ncpu=8
DIR=os.getcwd()
sampleList=[
'SRR10433204',
'SRR10433205',
'SRR10433206',
'SRR10433207',
'SRR10433224',
'SRR10433208',
'SRR10433228',
'SRR10433212',
'SRR10433213',
'SRR10433214',
'SRR10433218',
'SRR10433219',
'SRR10433220'
]

def convertSRA():
    ## pbs submission
    for srrID in sampleList:

        ## qsub
	runScript='fastq-dump --gzip -I --split-files '+DIR+'/SRA/'+srrID+'.1 --outdir '+DIR+'/SRA/'
        #print run_script

        output1=open('PBS_fastqdump_'+srrID+'.pbs','w')
        output1.write('#PBS -N fastqdump_'+srrID+'\n')
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
        output1.write(runScript+'\n')
        output1.write('\n')
        output1.write('sleep 30\n')
        output1.write('exit 0')
        output1.close()

        os.system('qsub PBS_fastqdump_'+srrID+'.pbs')



convertSRA()


