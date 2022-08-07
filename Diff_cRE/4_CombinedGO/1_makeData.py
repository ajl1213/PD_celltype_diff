#!/python2


import os
import numpy
import operator


DIR=os.getcwd()
RnaGO_DIR='/home/ajl1213/Projects/PD/SciAdv_Anal/Diff_RNA/5_runGO/GOBP'
GREAT_DIR='/home/ajl1213/Projects/PD/SciAdv_Anal/Diff_cRE/3_Great/GOBP'

dataTypeList=['Down','Up']
cellTypeList=['DopaN','GabaN','Oligo','OPC','Ast','Micro','Endo','Peri']
maxTerm=5

os.system('mkdir '+DIR+'/PvalMat')



def load_RNA_GO(dataType, cellType):
    RNA_pvalDict={}
    RNA_goDict={}
    RNA_geneDict={}

    input1=open(RnaGO_DIR+'/'+dataType+'DEGs.GO.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	termID='_'.join(each[0].split(' ')[:-1])
	goID=each[0].split(' ')[-1].replace('(','').replace(')','')

	geneList=each[8]
	cellLabel=each[10]
	logPval=-numpy.log10(float(each[2]))
	if cellType == cellLabel:
	    RNA_pvalDict[goID]=logPval
	    RNA_goDict[goID]=[termID, geneList]
	    RNA_geneDict[goID]=geneList.split(';')

    input1.close()

    return RNA_pvalDict, RNA_goDict, RNA_geneDict


def load_GREAT_GO(dataType, cellType):
    great_pvalDict={}
    great_goDict={}
    great_geneDict={}

    input1=open(GREAT_DIR+'/'+dataType+'.'+cellType+'.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	termID='_'.join(each[2].split(' '))
	goID=each[1]

	geneList=each[9]
	logPval=-numpy.log10(float(each[3]))

	great_pvalDict[goID]=logPval
	great_goDict[goID]=[termID, geneList]
	great_geneDict[goID]=geneList.split(',')

    input1.close()

    return great_pvalDict, great_goDict, great_geneDict


def makePvalMat():
    for dataType in dataTypeList:
	output1=open(DIR+'/PvalMat/'+dataType+'.logPval.txt','w')
	output1.write('cellType\ttermID\tsharedGenes\tMergedDEGs\tdiffPeaks\n')
	for cellType in cellTypeList:

	    # load GO
	    RNA_pvalDict, RNA_goDict, RNA_geneDict = load_RNA_GO(dataType, cellType)
	    great_pvalDict, great_goDict, great_geneDict = load_GREAT_GO(dataType, cellType)

	    # select GO
	    GoList=[]
	    for key in sorted(RNA_pvalDict.items(), key=operator.itemgetter(1), reverse=True):

		goID=key[0]
		pval=key[1]

		if great_pvalDict.has_key(goID):
		    if len(GoList) < maxTerm:
			GoList.append(goID)

	    # create pval matrix
	    for i in GoList:
		termID=RNA_goDict[i][0]

		RNA_pval=RNA_pvalDict[i]
		great_pval=great_pvalDict[i]

		commonGeneList=set(RNA_geneDict[i]) & set(great_geneDict[i])
		if len(commonGeneList) > 0:
		    finalGeneList=';'.join(commonGeneList)
		else:
		    finalGeneList='none'


		newLine=[cellType, termID, finalGeneList, str(RNA_pval), str(great_pval)]
		output1.write('\t'.join(newLine)+'\n')


makePvalMat()



