#!/home/ajl1213/anaconda2/bin/python



import os
import sys
import pysam


DIR=os.getcwd()
sampleID=sys.argv[1]
cellType=sys.argv[2]
BarcodeFile=sys.argv[3]
BamInput=sys.argv[4]


###
barcodeDict = {}
input1 = open(BarcodeFile,'r')
all_input1=input1.readlines()
for line in all_input1[1:]:
    each=line.strip().split('\t')
    code_transform = each[0].split('-')[0] + '-1'
    barcodeDict[code_transform]=1
input1.close()

###
countDict={}
f = pysam.AlignmentFile(BamInput, "rb")
out = pysam.AlignmentFile(DIR+'/BamFiles/'+cellType+'.'+sampleID+'.filter.bam','wb',template=f)
for line in f:
    if line.has_tag('CB'):
	if barcodeDict.has_key(line.get_tag('CB')):

	    tag=line.get_tag('CB')
	    out.write(line)

	    if countDict.has_key(tag):
		countDict[tag]+=1
	    else:
		countDict[tag]=1

f.close()
out.close()

###
output1=open(DIR+'/AlignStat/'+cellType+'.'+sampleID+'.cellList.txt','w')
output1.write('cellID\treadCount\n')
for cellID in countDict:
    output1.write(cellID+'\t'+str(countDict[cellID])+'\n')
output1.close()



