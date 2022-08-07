#!/home/ajl1213/anaconda2/bin/python


import os
import gzip


ncpu=4
DIR=os.getcwd()
fithicDIR='/home/ajl1213/programs/fithic/fithic'
pythonDIR='/home/ajl1213/anaconda2/bin'
coverageDIR='/home/ajl1213/Projects/PD/SciAdv_Data/HiC/CovNorm/Coverage'
inputDIR='/home/ajl1213/Projects/PD/SciAdv_Data/HiC/CovNorm/CisOnlyBam'

sampleList=['MergedNOSN','MergedPDSN']
chrList=['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',  'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16','chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']

maxDist=1000000
resolution='5kb'
#resolution='10kb'


os.system('mkdir '+DIR+'/Intersfile')
os.system('mkdir '+DIR+'/Fragsfile')
os.system('mkdir '+DIR+'/OUTPUT')



def runFithic(sampleID, chrID):
    output1=open('runFithic.'+sampleID+'.'+chrID+'.'+resolution+'.pbs','w')
    output1.write('#PBS -N fithic.'+sampleID+'.'+chrID+'.'+resolution+'\n')
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
    output1.write('\n')
    output1.write('samtools view '+inputDIR+'/'+sampleID+'.'+chrID+'.bam | python bam2inter.py '+str(maxDist)+' '+resolution+' | gzip -  > '+DIR+'/Intersfile/intersfile.'+sampleID+'.'+chrID+'.'+resolution+'.gz\n')
    output1.write('python getFragsfile.py '+coverageDIR+' '+sampleID+' '+chrID+' '+resolution+'\n')
    output1.write(pythonDIR+'/python '+fithicDIR+'/fithic.py -i '+DIR+'/Intersfile/intersfile.'+sampleID+'.'+chrID+'.'+resolution+'.gz -f '+DIR+'/Fragsfile/fragsfile.'+sampleID+'.'+chrID+'.'+resolution+'.gz -o '+DIR+'/OUTPUT/'+sampleID+'.'+chrID+'.'+resolution+' -U '+str(maxDist)+' -r 0 -p 2\n')
    output1.write('\n')
    output1.write('sleep 30\n')
    output1.write('exit 0')
    output1.close()


for sampleID in sampleList:
    for chrID in chrList:
	runFithic(sampleID, chrID)
	os.system('qsub runFithic.'+sampleID+'.'+chrID+'.'+resolution+'.pbs')



