#!/home/ajl1213/anaconda2/bin/python


## Individual the Normal and PD Hi-C interactions

import os
import numpy
import gzip


DIR=os.getcwd()
gtf='/home/ajl1213/genome.info/gencode/hg19/gencode.v19.ptn_coding.lincRNA.level_1_2'
inputDIR='/home/ajl1213/Projects/PD/SciAdv_Data/HiC/FitHiC/SigInter'
sampleList=['MergedNOSN', 'MergedPDSN']
resolution='5kb'
qvalCut='1e-3'

os.system('mkdir '+DIR+'/GenomeTrack')


def loadGeneInfo():
    geneDict={}
    input1=open(gtf,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	ensembleID=each[0]
	geneID=each[1]
	geneDict[ensembleID]=geneID
    input1.close()

    return geneDict



def getInter():
    geneDict=loadGeneInfo()

    ##
    for sampleID in sampleList:
        input1=open(inputDIR+'/'+sampleID+'.'+resolution+'.'+qvalCut+'.FitHiC.all.txt','r')
        all_input1=input1.readlines()
        output1=open(DIR+'/GenomeTrack/'+sampleID+'.Inter.washu.txt','w')
        index1=1
        for line in all_input1[1:]:
            each=line.strip().split('\t')

	    frag1=each[0]
	    frag2=each[1]

	    chrID=frag1.split('.')[0]
	    frag1Pos1=frag1.split('.')[1]
	    frag1Pos2=frag1.split('.')[2]
	    frag2Pos1=frag2.split('.')[1]
	    frag2Pos2=frag2.split('.')[2]

	    qval=each[4]
	    logQval=-numpy.log10(float(qval))

	    ## 
	    output1.write(chrID+'\t'+str(frag1Pos1)+'\t'+str(frag1Pos2)+'\t'+chrID+':'+str(frag2Pos1)+'-'+str(frag2Pos2)+','+str(logQval)+'\t'+str(index1)+'\t.\n')
	    index1+=1
	    output1.write(chrID+'\t'+str(frag2Pos1)+'\t'+str(frag2Pos2)+'\t'+chrID+':'+str(frag1Pos1)+'-'+str(frag1Pos2)+','+str(logQval)+'\t'+str(index1)+'\t.\n')
	    index1+=1

        input1.close()
        output1.close()

        ##
        runLine1='cat '+DIR+'/GenomeTrack/'+sampleID+'.Inter.washu.txt | sort -k1,1 -k2,2n - > '+DIR+'/GenomeTrack/'+sampleID+'.Inter.washu.sorted.txt'
        os.system(runLine1)
        runLine2='bgzip '+DIR+'/GenomeTrack/'+sampleID+'.Inter.washu.sorted.txt'
        os.system(runLine2)
        runLine3='tabix -p bed '+DIR+'/GenomeTrack/'+sampleID+'.Inter.washu.sorted.txt.gz'
        os.system(runLine3)




getInter()



