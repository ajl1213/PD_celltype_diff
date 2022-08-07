#!/home/ajl1213/anaconda2/bin/python


import os
import numpy


DIR=os.getcwd()

inputDIR='/home/ajl1213/Projects/PD/SciAdv_Anal/CombinedTargetID/01_FindTarget/RESULT'
GTF='/home/ajl1213/genome.info/gencode/hg19/gencode.v19.ptn_coding.lincRNA.level_1_2'
AbcScoreDIR='/home/ajl1213/Projects/PD/data/HiC/ActivityByContact/ActivityByContactScore'
TargetDIR=DIR+'/MergedTargetList'
cellTypeList=['DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']

minAbcScore=10
minCor=0.3



def loadGeneInfo():
    geneDict={}
    input1=open(GTF,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        ensembleID=each[0]
        geneID=each[1]
        geneDict[ensembleID]=geneID
    input1.close()

    return geneDict


def loadContactScore(cellType):
    contactValDict={}
    input1=open(AbcScoreDIR+'/'+cellType+'.ActivityByContactScore.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        ensembleID=each[0]
        creID=each[3]
        contactVal=float(each[8])
        abcScore=float(each[9])

        keyID=ensembleID+';'+creID
        valID=[contactVal, abcScore]

        contactValDict[keyID]=valID
    input1.close()

    return contactValDict


def select():
    for cellType in cellTypeList:
        print cellType

        input1=open(inputDIR+'/'+cellType+'.fullList.txt','r')
        all_input1=input1.readlines()
        output1=open(inputDIR+'/'+cellType+'.selectList.txt','w')
        output1.write('ensembleID\tgeneID\tgeneLabel\tPdGeneLabel\tfuncGeneLabel\tcreID\tfuncPeakLabel\tcreType\tdist\tcorVal\tcontactVal\tabcScore\n')
        for line in all_input1[1:]:
            [ensembleID, geneID, geneLabel, PdGeneLabel, funcGeneLabel, creID, funcPeakLabel, creType, dist, corVal, contactVal, abcScore]=line.strip().split('\t')
            if (cellType in funcGeneLabel.split(',') and cellType in funcPeakLabel.split(',')):
                if float(corVal) > minCor:
#                   if float(abcScore) > minAbcScore:
                        newLine=[ensembleID, geneID, geneLabel, PdGeneLabel, funcGeneLabel, creID, funcPeakLabel, creType, dist, corVal, contactVal, abcScore]
                        output1.write('\t'.join(newLine)+'\n')

                        if float(abcScore) > minAbcScore:
                            if PdGeneLabel == 'knownPdGene':
                                print '\t'.join(newLine)    

        input1.close()
        output1.close()



select()




