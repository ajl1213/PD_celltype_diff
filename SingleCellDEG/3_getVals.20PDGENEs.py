#!python2


import os
import numpy
import operator


DIR=os.getcwd()
snRNA_file='/home/ajl1213/Projects/PD/SciAdv_Data/snRNAseq/02_SeuratProcess/CellTypeDEGs/snRNA.cellTypeDEGs.txt'
cellTypeList=['DopaN','GabaN','Oligo', 'OPC', 'Ast', 'Micro','Endo','Peri']
PDGeneList = [
'SNCA',
'MAPT',
'UCHL1',
'PINK1',
'PARK7',
'LRRK2',
'GBA',
'SCARB2',
'ATP13A2',
'FBXO7',
'VPS35',
'DNAJC6',
'SYNJ1',
'CHCHD2',
'VPS13C',
'GAK',
'SMPD1',
'SLC17A5',
'ATP6V0A1',
'CTSB'
]


def load_rna():
    valDict={}
    input1=open(snRNA_file,'r')
    all_input1=input1.readlines()
    for cellType in cellTypeList:
	tmpDict={}
	for line in all_input1[1:]:
	    each=line.strip().split()

	    pct1=float(each[2])
	    pct2=float(each[3])
	    log2fc=float(each[1])
	    adjPval=float(each[4])
	    cluster=each[6]
	    geneID=each[7]

	    if cluster==cellType:
		tmpDict[geneID]=[log2fc, adjPval, pct1, pct2]

	valDict[cellType]=tmpDict

    input1.close()

    return valDict


def makeValMat():
    valDict=load_rna()

    # Value list
    output1=open(DIR+'/PDGenes.valList.txt','w')
    output1.write('GeneID\tCellType\tLog2FC\tAdjPval\tpct1\tpct2\n')
    for geneID in PDGeneList:
	for cellType in cellTypeList:

	    log2fc = valDict[cellType][geneID][0]
	    adjPval = valDict[cellType][geneID][1]
	    pct1 = valDict[cellType][geneID][2]
	    pct2 = valDict[cellType][geneID][3]

	    newLine=[geneID, cellType, str(log2fc), str(adjPval), str(pct1), str(pct2)]
	    output1.write('\t'.join(newLine)+'\n')

    output1.close()

    # Log2FC Matrix
    output1=open(DIR+'/PDGenes.valMat.txt','w')
    output1.write('GeneID\t'+'\t'.join(cellTypeList)+'\n')
    for geneID in PDGeneList:
	tmpList=[]
	for cellType in cellTypeList:

	    log2fc = valDict[cellType][geneID][0]
	    adjPval = valDict[cellType][geneID][1]

	    tmpList.append(str(log2fc))

	output1.write(geneID+'\t'+'\t'.join(tmpList)+'\n')
    output1.close()


makeValMat()





