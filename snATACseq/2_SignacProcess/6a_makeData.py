#!python2



import os
import numpy



DIR=os.getcwd()
cellTypeList=['DopaN','GabaN','Oligo','OPC','Ast','Micro','Endo','Peri']
inFile=DIR+'/CellTypeDEGs/snATAC.cellTypeMarkers.txt'
outFile=DIR+'/CellTypeDEGs/snATAC.cellTypeMarkers.inputPlot.txt'
nGene=5

input1=open(inFile,'r')
all_input1=input1.readlines()
geneList=[]
valDict={}
for cell in cellTypeList:
    index=0
    tmpDict={}
    for line in all_input1[1:]:
	each=line.strip().split('\t')

	avgLog2FC=float(each[1])
	adjPval=float(each[4])
	if adjPval == 0:
	    adjPval = 2.225074e-308 
	    logAdjPval = -numpy.log10(adjPval)
	else:
	    logAdjPval = -numpy.log10(adjPval)
	    
	celltype=each[5]
	geneID=each[6]

	if cell == celltype:
	    tmpDict[geneID] = [avgLog2FC, logAdjPval]

	    if index < nGene:
		if not geneID in geneList:
		    geneList.append(geneID)
		    index+=1	    

    valDict[cell] = tmpDict
    input1.close()

##
output1=open(outFile,'w')
output1.write('geneID\tcelltype\tavgLogFC\tlogAdjPval\n')
for geneID in geneList:
    for cellType in cellTypeList:
	if valDict[cellType].has_key(geneID):
	    avgLog2FC=valDict[cellType][geneID][0]
	    logAdjPval=valDict[cellType][geneID][1]
	#else:
	#    avgLog2FC=0
	#    logAdjPval=0

	    newLine=[geneID, cellType, str(avgLog2FC), str(logAdjPval)]
	    output1.write('\t'.join(newLine)+'\n')
output1.close()



