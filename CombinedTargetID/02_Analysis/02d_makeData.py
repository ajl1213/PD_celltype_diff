#!/home/ajl1213/anaconda2/bin/python


import os
import numpy
import operator



DIR=os.getcwd()
inputDIR=DIR+'/GOTABLE'
cellTypeList=['DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']

termList=[
'hyperactivity',
'hypoactivity',
'tremors',
'decreased grip strength',
'abnormal gait',
'limb grasping',
'impaired coordination',
'decreased startle reflex',
'decreased prepulse inhibition',
'ataxia',
'abnormal motor capabilities/coordination/movement',
'abnormal neural tube morphology',
'abnormal brain morphology',
'decreased vertical activity',
'decreased lymphocyte cell number',
'abnormal locomotor behavior',
'decreased CD8-positive, alpha-beta T cell number',
'decreased T cell number',
'abnormal axon morphology',
'abnormal myelination',
'increased monocyte cell number',
'abnormal locomotor activation',
'demyelination',
'astrocytosis',
'increased macrophage cell number',
'impaired righting response',
'increased neuron apoptosis',
'abnormal neural tube closure'
]



def makeData():

    ## Load data
    allDict={}
    for cellType in cellTypeList:
	input1=open(inputDIR+'/'+cellType+'.MP.txt', 'r')
	all_input1=input1.readlines()
        tmpDict={}
        for line in all_input1[1:]:
            each=line.strip().split('\t')
	    termID=' '.join(each[0].split(' ')[:-1])

            pval=each[2]
            adjPval=each[3]
            oddRatio=each[6]
            combinedScore=each[7]

	    tmpDict[termID]=[pval, oddRatio, combinedScore]

        input1.close()

        allDict[cellType]=tmpDict

    ## pval
    output1=open(DIR+'/GOTABLESORTED/MP.pvalList.txt','w')
    output1.write('termID\tcellType\tpval\n')
    for termID in termList:
        for cellType in cellTypeList:
            pval=str(allDict[cellType][termID][0]) if allDict[cellType].has_key(termID) else '1.0'
            newLine=[termID, cellType, pval]
            output1.write('\t'.join(newLine)+'\n')
    output1.close()

    ## plot matrix
    output1=open(DIR+'/GOTABLESORTED/MP.celltype.txt','w')
    output1.write('termID\t'+'\t'.join(cellTypeList)+'\n')
    for termID in termList:
        valList=[]
        for cellType in cellTypeList:
            if allDict[cellType].has_key(termID):
                [pval, oddRatio, combinedScore]=allDict[cellType][termID]
                logPval=-numpy.log10(float(pval))
                valList.append(str(logPval))
            else:
                valList.append('0.0')
        output1.write(termID+'\t'+'\t'.join(valList)+'\n')

    output1.close()


makeData()


