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

		## Bulk DEG
		# intersect
		runScript='intersectBed -a '+DIR+'/DegList/'+dataType+'.'+cellType+'.'+winDist+'.bulkDEG.bed -b '+PeakDIR+'/'+dataType+'.'+cellType+'.bed -wa -wb | cut -f4 | sort -u | wc -l > '+DIR+'/ObservedData/'+dataType+'.'+cellType+'.'+winDist+'.bulkDEG.Intersect.bed'
		os.system(runScript)
		runScript='intersectBed -a '+DIR+'/DegList/'+dataType+'.'+cellType+'.'+winDist+'.bulkDEG.bed -b '+PeakDIR+'/'+dataType+'.'+cellType+'.bed -wa -wb > '+DIR+'/ObservedData/'+dataType+'.'+cellType+'.'+winDist+'.bulkDEG.Intersect.full.bed'
		os.system(runScript)

		# get DEG number
		nTotalBulkDEG=0
		input1=open(DIR+'/DegList/'+dataType+'.'+cellType+'.'+winDist+'.bulkDEG.bed','r')
		all_input1=input1.readlines()
		for line in all_input1:
		    nTotalBulkDEG+=1
		input1.close()

		# load actual overlap
		input1=open(DIR+'/ObservedData/'+dataType+'.'+cellType+'.'+winDist+'.bulkDEG.Intersect.bed', 'r')
		all_input1=input1.readlines()
		nBulkDEG=all_input1[0].strip().split('\t')[0]
		bulkDegFrac=float(nBulkDEG)/float(nTotalBulkDEG)*100 if not nTotalBulkDEG==0 else 0
		input1.close()

		# compute pvalue
		bulkRandValList=[]
		input1=open(DIR+'/RandDegPermute/'+dataType+'.'+cellType+'.'+winDist+'.RandDegPermuteRecord.txt','r')
		all_input1=input1.readlines()
		for line in all_input1[1:]:
		    each=line.strip().split('\t')
		    DEG_type=each[1]
		    fracVal=float(each[3])
		    if DEG_type=='bulkDEG':
			bulkRandValList.append(fracVal)
		input1.close()	

		bulkTmpList=[]
		for i in bulkRandValList:
		    if i > bulkDegFrac:
			bulkTmpList.append(i) 

		bulkPval=float(len(bulkTmpList))/float(len(bulkRandValList)) if not len(bulkRandValList)==0 else 1
		meanBulkRandVal=numpy.mean(bulkRandValList) if not len(bulkRandValList)==0 else 0

		#
		newLine=[dataType, cellType, winDist, 'bulkDEG', str(nBulkDEG), str(bulkDegFrac), str(meanBulkRandVal), str(bulkPval)]
		output1.write('\t'.join(newLine)+'\n')


		## snRNA DEG
		# intersect
		runScript='intersectBed -a '+DIR+'/DegList/'+dataType+'.'+cellType+'.'+winDist+'.snDEG.bed -b '+PeakDIR+'/'+dataType+'.'+cellType+'.bed -wa -wb | cut -f4 | sort -u | wc -l > '+DIR+'/ObservedData/'+dataType+'.'+cellType+'.'+winDist+'.snDEG.Intersect.bed'
		os.system(runScript)
		runScript='intersectBed -a '+DIR+'/DegList/'+dataType+'.'+cellType+'.'+winDist+'.snDEG.bed -b '+PeakDIR+'/'+dataType+'.'+cellType+'.bed -wa -wb > '+DIR+'/ObservedData/'+dataType+'.'+cellType+'.'+winDist+'.snDEG.Intersect.full.bed'
		os.system(runScript)

		# get DEG number
		totalSnDEG=0
		input1=open(DIR+'/DegList/'+dataType+'.'+cellType+'.'+winDist+'.snDEG.bed','r')
		all_input1=input1.readlines()
		for line in all_input1:
		    totalSnDEG+=1
		input1.close()

		# load actual overlap
		input1=open(DIR+'/ObservedData/'+dataType+'.'+cellType+'.'+winDist+'.snDEG.Intersect.bed', 'r')
		all_input1=input1.readlines()
		nSnDEG=all_input1[0].strip().split('\t')[0]
		snDegFrac=float(nSnDEG)/float(totalSnDEG)*100 if not totalSnDEG==0 else 0
		input1.close()


		# compute pvalue
		snRandValList=[]
		input1=open(DIR+'/RandDegPermute/'+dataType+'.'+cellType+'.'+winDist+'.RandDegPermuteRecord.txt','r')
		all_input1=input1.readlines()
		for line in all_input1[1:]:
		    each=line.strip().split('\t')
		    DEG_type=each[1]
		    fracVal=float(each[3])
		    if DEG_type=='snDEG':
			snRandValList.append(fracVal)
		input1.close()	

		snTmpList=[]
		for i in snRandValList:
		    if i > snDegFrac:
			snTmpList.append(i)

		snPval=float(len(snTmpList))/float(len(snRandValList)) if not len(snRandValList)==0 else 1
		meanSnRandVal=numpy.mean(snRandValList) if not len(bulkRandValList)==0 else 0

		#
		newLine=[dataType, cellType, winDist, 'snDEG', str(nSnDEG), str(snDegFrac), str(meanSnRandVal), str(snPval)]
		output1.write('\t'.join(newLine)+'\n')

    output1.close()




getActualData()



