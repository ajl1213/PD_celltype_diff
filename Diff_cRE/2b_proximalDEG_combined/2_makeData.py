#!/python2


import os
import numpy
import random



DIR=os.getcwd()
PeakDIR='/home/ajl1213/Projects/PD/NeuronJ/DiffAnal/02_H3K27ac/1_BulkH3K27ac/BedFiles'
dataTypeList=['Down','Up']
cellTypeList=['DopaN','GabaN','Oligo','OPC','Ast','Micro','Endo','Peri']
windowDict={
'100kb':100000
}

os.system('mkdir '+DIR+'/ObservedData')
os.system('mkdir '+DIR+'/RESULT')


def getActualData():

    output1=open(DIR+'/RESULT/EnrichStatFinal.txt','w')
    output1.write('dataType\tcellType\twinDist\tDEG_type\tnPeakDEG\tpeakDegFrac\tmeanRandomFrac\tpval\n')
    for dataType in dataTypeList:
        for cellType in cellTypeList:
            for winDist in windowDict:
		print dataType, cellType, winDist

		## Merged DEGs
		# intersect
		runScript='intersectBed -a '+PeakDIR+'/'+dataType+'.'+cellType+'.bed -b '+DIR+'/DegList/'+dataType+'.'+cellType+'.'+winDist+'.MergedDEG.bed -wa -wb | cut -f8 | sort -u | wc -l > '+DIR+'/ObservedData/'+dataType+'.'+cellType+'.'+winDist+'.MergedDEG.Intersect.bed'
		os.system(runScript)
		runScript='intersectBed -a '+PeakDIR+'/'+dataType+'.'+cellType+'.bed -b '+DIR+'/DegList/'+dataType+'.'+cellType+'.'+winDist+'.MergedDEG.bed -wa -wb > '+DIR+'/ObservedData/'+dataType+'.'+cellType+'.'+winDist+'.MergedDEG.Intersect.full.bed'
		os.system(runScript)

		# get DEG number
		nTotalDEG=0
		input1=open(DIR+'/DegList/'+dataType+'.'+cellType+'.'+winDist+'.MergedDEG.bed','r')
		all_input1=input1.readlines()
		for line in all_input1:
		    nTotalDEG+=1
		input1.close()

		# load actual overlap
		input1=open(DIR+'/ObservedData/'+dataType+'.'+cellType+'.'+winDist+'.MergedDEG.Intersect.bed', 'r')
		all_input1=input1.readlines()
		nDEG=all_input1[0].strip().split('\t')[0]
		DegFrac=float(nDEG)/float(nTotalDEG)*100 if not nTotalDEG==0 else 0
		input1.close()

		# compute pvalue
		RandValList=[]
		input1=open(DIR+'/RandPeakPermute/'+dataType+'.'+cellType+'.'+winDist+'.RandPeakPermuteRecord.txt','r')
		all_input1=input1.readlines()
		for line in all_input1[1:]:
		    each=line.strip().split('\t')
		    DEG_type=each[1]
		    fracVal=float(each[3])
		    if DEG_type=='MergedDEG':
			RandValList.append(fracVal)
		input1.close()	

		TmpList=[]
		for i in RandValList:
		    if i > DegFrac:
			TmpList.append(i) 

		pval=float(len(TmpList))/float(len(RandValList)) if not len(RandValList)==0 else 1
		meanBulkRandVal=numpy.mean(RandValList) if not len(RandValList)==0 else 0

		#
		newLine=[dataType, cellType, winDist, 'MergedDEG', str(nDEG), str(DegFrac), str(meanBulkRandVal), str(pval)]
		output1.write('\t'.join(newLine)+'\n')


    output1.close()




getActualData()



