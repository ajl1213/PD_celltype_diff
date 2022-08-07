#!/home/ajl1213/anaconda2/bin/python


import os



DIR=os.getcwd()
ncpu=4


def makeScript():

    ### Public1NOSN
    for sampleID in ['SRR11442505.1_1','SRR11442505.1_2','SRR11442505.1_3','SRR11442505.1_4']:
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
	output1.write('zcat Public1NOSN/'+sampleID+'.fastq.gz | awk -F \' \' \'{if($1~/SRR/){split($1,a,\".\"); print a[1]\".\"a[2]\".\"a[3]\" \"$2\" \"$3} else{print $0}}\' | gzip > Public1NOSN/'+sampleID+'.mod.fastq.gz\n')
	output1.write('\n')
	output1.write('sleep 30\n')
	output1.write('exit 0')
	output1.close()

	os.system('qsub PBS_scATAC_'+sampleID+'.pbs')

    ### Public2NOSN
    for sampleID in ['SRR11442506.1_1','SRR11442506.1_2','SRR11442506.1_3','SRR11442506.1_4']:
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
	output1.write('zcat Public2NOSN/'+sampleID+'.fastq.gz | awk -F \' \' \'{if($1~/SRR/){split($1,a,\".\"); print a[1]\".\"a[2]\".\"a[3]\" \"$2\" \"$3} else{print $0}}\' | gzip > Public2NOSN/'+sampleID+'.mod.fastq.gz\n')
	output1.write('\n')
	output1.write('sleep 30\n')
	output1.write('exit 0')

	output1.close()

	os.system('qsub PBS_scATAC_'+sampleID+'.pbs')



makeScript()




