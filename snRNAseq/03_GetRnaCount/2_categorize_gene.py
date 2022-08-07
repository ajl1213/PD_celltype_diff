#!/python2


import os
import numpy


DIR=os.getcwd()
inFile=DIR+'/CountByCellType/ExpRatio.byCellType.txt'

os.system('mkdir '+DIR+'/CellLabel')
os.system('mkdir '+DIR+'/Plots')


def assignCellType():
    input1=open(inFile,'r')
    all_input1=input1.readlines()
    cellTypeList=all_input1[0].strip().split('\t')[1:]
    output1=open(DIR+'/CellLabel/GeneExp.CellLabel.txt','w')
    output1.write('geneID\tfunctionalCell\n')
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        geneID=each[0]
        vals=numpy.array(each[1:], dtype='float')

        if max(vals) > 10:
            index=0
            cellList=[]
            for val in vals:
                if val > 10:
                    assignedCell = cellTypeList[index]
                    cellList.append(assignedCell)
                index+=1

	    if len(cellList) > 0:
		output1.write(geneID+'\t'+','.join(cellList)+'\n')
	    else:
		output1.write(geneID+'\tnone\n')

    input1.close()
    output1.close()

    ## make separate list
    for cellType in cellTypeList:
        input1=open(DIR+'/CellLabel/GeneExp.CellLabel.txt','r')
        all_input1=input1.readlines()
        output1=open(DIR+'/CellLabel/'+cellType+'.GeneExp.txt','w')
        for line in all_input1[1:]:
            each=line.strip().split('\t')
            geneID=each[0]
            assignedCell=each[1]
            if assignedCell.count(cellType) > 0:
                output1.write(geneID+'\t'+assignedCell+'\n')
        input1.close()
        output1.close()


def makeChowRuskeyInput():
    a=[]
    b=[]
    c=[]
    d=[]
    e=[]

    input1=open(DIR+'/CellLabel/DopaN.GeneExp.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        geneID=each[0]
        a.append(geneID)
    input1.close()

    input1=open(DIR+'/CellLabel/Oligo.GeneExp.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        geneID=each[0]
        b.append(geneID)
    input1.close()

    input1=open(DIR+'/CellLabel/Ast.GeneExp.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        geneID=each[0]
        c.append(geneID)
    input1.close()

    input1=open(DIR+'/CellLabel/Micro.GeneExp.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        geneID=each[0]
        d.append(geneID)
    input1.close()

    input1=open(DIR+'/CellLabel/Endo.GeneExp.txt','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
        each=line.strip().split('\t')
        geneID=each[0]
        e.append(geneID)
    input1.close()

    a=set(a)
    b=set(b)
    c=set(c)
    d=set(d)
    e=set(e)

    ## compute combinations
    # 00000
    c00000=0
    # 10000
    c10000=len(a-b-c-d-e)
    # 01000
    c01000=len(b-a-c-d-e)
    # 11000
    c11000=len(a&b-c-d-e)
    # 00100
    c00100=len(c-a-b-d-e)
    # 10100
    c10100=len(a&c-b-d-e)
    # 01100
    c01100=len(b&c-a-d-e)
    # 11100
    c11100=len(a&b&c-d-e)
    # 00010
    c00010=len(d-a-b-c-e)
    # 10010
    c10010=len(a&d-b-c-e)
    # 01010
    c01010=len(b&d-a-c-e)
    # 11010
    c11010=len(a&b&d-c-e)
    # 00110 
    c00110=len(c&d-a-b-e)
    # 10110
    c10110=len(a&c&d-b-e)
    # 01110
    c01110=len(b&c&d-a-e)
    # 11110
    c11110=len(a&b&c&d-e)
    # 00001
    c00001=len(e-a-b-c-d)
    # 10001
    c10001=len(a&e-b-c-d)
    # 01001
    c01001=len(b&e-a-c-d)
    # 11001
    c11001=len(a&b&e-c-d)
    # 00101
    c00101=len(c&e-a-b-d)
    # 10101
    c10101=len(a&c&e-b-d)
    # 01101
    c01101=len(b&c&e-a-d)
    # 11101
    c11101=len(a&b&c&e-d)
    # 00011
    c00011=len(d&e-a-b-c)
    # 10011 
    c10011=len(a&d&e-b-c)
    # 01011
    c01011=len(b&d&e-a-c)
    # 11011
    c11011=len(a&b&d&e-c)
    # 00111
    c00111=len(c&d&e-a-b)
    # 10111
    c10111=len(a&c&d&e-b)
    # 01111
    c01111=len(b&c&d&e-a)
    # 11111 
    c11111=len(a&b&c&d&e)

    ##
    labelList=['c00000','c10000','c01000','c11000','c00100','c10100','c01100','c11100','c00010','c10010','c01010','c11010','c00110','c10110','c01110','c11110','c00001','c10001','c01001','c11001','c00101','c10101','c01101','c11101','c00011','c10011','c01011','c11011','c00111','c10111','c01111','c11111']
    valList=[c00000,c10000,c01000,c11000,c00100,c10100,c01100,c11100,c00010,c10010,c01010,c11010,c00110,c10110,c01110,c11110,c00001,c10001,c01001,c11001,c00101,c10101,c01101,c11101,c00011,c10011,c01011,c11011,c00111,c10111,c01111,c11111]

    output1=open(DIR+'/Plots/VennInput.txt','w')
    output1.write('comb\tnVal\n')
    for i in range(len(valList)):
        output1.write(str(labelList[i])+'\t'+str(valList[i])+'\n')

    output1.close()



assignCellType()
makeChowRuskeyInput()



