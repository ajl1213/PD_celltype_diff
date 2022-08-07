#!/home/ajl1213/anaconda2/bin/python


import os


DIR=os.getcwd()
PICARDJAR='/home/ajl1213/programs/picard_tool/picard.jar'
cpu=8

sampleList=[
'X4870NOSN-X4689NOSN-X5628NOSN-X5130NOSN-X5006NOSN-X4996NOSN;MergedNOSN',
'X5244PDSN-X4831PDSN-X5215PDSN-X5649PDSN-X5591PDSN;MergedPDSN'
]

def makeScript(i):
    sampleList=i.split(';')[0].split('-')
    newLabel=i.split(';')[1]
    output1=open('PBS_hicMerge_'+newLabel+'.pbs','w')
    output1.write('#PBS -N hicMerge_'+newLabel+'\n')
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
    output1.write('samtools merge SortedNodupBam/'+newLabel+'.sorted.nodup.bam')
    for sampleID in sampleList:
	output1.write(' SortedNodupBam/'+sampleID+'.sorted.nodup.bam')
    output1.write('\n')
    output1.write('samtools view -S SortedNodupBam/'+newLabel+'.sorted.nodup.bam | python check_trans_cis.py TransCis/'+newLabel+'.TransCis.txt\n')
    output1.write('samtools index SortedNodupBam/'+newLabel+'.sorted.nodup.bam\n')

    output1.write('\n')
    output1.write('\n')
    output1.write('sleep 30\n')
    output1.write('exit 0')
    output1.close()


for i in sampleList:
    makeScript(i)
    newLabel=i.split(';')[1]
    os.system('qsub PBS_hicMerge_'+newLabel+'.pbs')



