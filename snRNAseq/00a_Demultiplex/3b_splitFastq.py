#!/home/ajl1213/anaconda2/bin/python


import os
import gzip
import sys 
import pysam



DIR=os.getcwd()
seqID=sys.argv[1]
sampleID=sys.argv[2]
nLane=sys.argv[3]
readType=sys.argv[4]
barcodeFile=sys.argv[5]
bamFile=sys.argv[6]
inFastqDIR=sys.argv[7]
outFastqDIR=sys.argv[8]


def loadBarcodes():

    barcodeDict={}
    input1=open(barcodeFile,'r')
    all_input1=input1.readlines()
    for line in all_input1:
	each=line.strip().split('\t')
	barcodeID=each[0]
	barcodeDict[barcodeID]=1
    input1.close()

    return barcodeDict


def getReadIDs():
    barcodeDict=loadBarcodes()

    readDict={}
    countDict={}
    input1 = pysam.AlignmentFile(bamFile,'rb')
    for line in input1:
	if line.has_tag('CB'):
	    tag=line.get_tag('CB')
	    if barcodeDict.has_key(tag):
		readID=line.query_name.split('/')[0]
		readDict[readID]=1

		if countDict.has_key(tag):
		    tmp=countDict[tag]
		    tmp.append(readID)
		    countDict[tag]=tmp
		else:
		    countDict[tag]=[readID]
    input1.close()

    output1=open(DIR+'/ReadPerCell/'+sampleID+'.'+readType+'.txt','w')
    output1.write('BarcodeID\treadCount\n')
    for tag in countDict:
	readCount=len(set(countDict[tag]))
	newLine=[tag, str(readCount)]
	output1.write('\t'.join(newLine)+'\n')
    output1.close()


    return readDict


def splitFastq():
    readDict=getReadIDs()

    inFastq = gzip.open(inFastqDIR+'/'+seqID+'_S1_L00'+nLane+'_'+readType+'_001.fastq.gz','rb')
    outFastq = gzip.open(outFastqDIR+'/'+sampleID+'_S1_L00'+nLane+'_'+readType+'_001.fastq.gz','wb')

    count1 = 1 
    count2 = 0 
    for read in inFastq:
	    if count1 == 1:
		    query = read.decode('utf-8').split(' ')[0].split('@')[1].split('/')[0]
		    if readDict.has_key(query):
			    outFastq.write(read)
			    count2 = 1 
			    del readDict[query]
		    count1 = 2 

	    elif count1 ==2:
		    if count2 ==1:
			    outFastq.write(read)
		    count1=3

	    elif count1 ==3:
		    if count2 ==1:
			    outFastq.write(read)
		    count1=4

	    elif count1==4:
		    if count2 ==1:
			    outFastq.write(read)
		    count1=1
		    count2=0

    inFastq.close()
    outFastq.close()



splitFastq()


