#!/home/ajl1213/anaconda2/bin/python


import os



ncpu=8
DIR=os.getcwd()
sampleList=[
'SRR11442505.1;Public1NOSN',
'SRR11442506.1;Public2NOSN',
]


def convert_sra():
    for sampleID in sampleList:
	srrID=sampleID.split(';')[0]
        newLabel=sampleID.split(';')[1]

	##### make directory
	os.system('mkdir '+DIR+'/'+newLabel)

        ##### qsub
	runScript='fastq-dump --gzip -I --split-files '+DIR+'/'+srrID+' --outdir '+DIR+'/'+newLabel
        #print run_script

        output1=open('PBS_fastqdump_'+newLabel+'.pbs','w')
        output1.write('#PBS -N fastqdump_'+newLabel+'\n')
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

        os.system('qsub PBS_fastqdump_'+newLabel+'.pbs')




convert_sra()


