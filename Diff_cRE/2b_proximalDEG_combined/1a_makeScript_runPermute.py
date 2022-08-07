#!/python2


import os
import sys 
import numpy
import random



DIR=os.getcwd()
ncpu=1
fai='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa.fai'
gtf='/home/ajl1213/genome.info/gencode/hg19.release38/GTF/gencode.v38.ptncoding.lncRNA.level_1_2.txt'
PeakDIR='/home/ajl1213/Projects/PD/SciAdv_Anal/Diff_cRE/1_BulkH3K27ac/BedFiles'
MergedDegDIR='/home/ajl1213/Projects/PD/NeuronJ/DiffAnal/01_RNA/4_DEG_overlap'
dataTypeList=['Down','Up']
cellTypeList=['DopaN','GabaN','Oligo','OPC','Ast','Micro','Endo','Peri']
windowDict={
'100kb':100000
}
maxIter=10000

os.system('mkdir '+DIR+'/DegList')
os.system('mkdir '+DIR+'/RandPeakPermute')



def makeDegList():
    # Load TSS
    tssDict={}
    input1=open(gtf,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        geneID=each[1]
        chrID=each[2]
        tss=int(each[3])
        geneLoc=chrID+':'+str(tss)
        tssDict[geneID]=geneLoc
    input1.close()

    # DEG
    for dataType in dataTypeList:
	for cellType in cellTypeList:
	    for winDist in windowDict:
		input1=open(MergedDegDIR+'/'+dataType+'.'+cellType+'.txt','r')
		all_input1=input1.readlines()
		output1=open(DIR+'/DegList/'+dataType+'.'+cellType+'.'+winDist+'.MergedDEG.bed','w')
		for line in all_input1[1:]:
		    each=line.strip().split('\t')
		    geneID=each[0]

		    if tssDict.has_key(geneID):
			tss=tssDict[geneID]
			chrID=tss.split(':')[0]
			pt=int(tss.split(':')[1])

			pt1=pt-windowDict[winDist]
			if pt1 < 0:
			    pt1=0
			pt2=pt+windowDict[winDist]

			newLine=[chrID, str(pt1), str(pt2), geneID]
			output1.write('\t'.join(newLine)+'\n')
		input1.close()


def makeScript():

    for dataType in dataTypeList:
	for cellType in cellTypeList:
	    for winDist in windowDict:

		output1=open('PBS_randPeakPermute.'+dataType+'.'+cellType+'.'+winDist+'.pbs','w')
		output1.write('#PBS -N randPeakPermute.'+dataType+'.'+cellType+'.'+winDist+'\n')
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

		output1.write('python '+DIR+'/1b_permuteDiffPeak.py '+dataType+' '+cellType+' '+winDist+' '+PeakDIR+' '+str(maxIter)+' '+fai+'\n')

		output1.write('\n')
		output1.write('sleep 30\n')
		output1.write('exit 0')

		output1.close()


		os.system('qsub PBS_randPeakPermute.'+dataType+'.'+cellType+'.'+winDist+'.pbs')


makeDegList()
makeScript()




