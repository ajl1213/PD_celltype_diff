#!/home/ajl1213/anaconda2/bin/python



import os
import numpy



DIR=os.getcwd()
coverageDIR='/home/ajl1213/Projects/PD/SciAdv_Data/ChIP/01_Mapping/PeakCoverage'
peakBED='/home/ajl1213/Projects/PD/data/snATAC/03_CellTypeBam/PeakList/snATAC.peaks.bed'
fripDIR='/home/ajl1213/Projects/PD/SciAdv_Data/ChIP/01_Mapping/FRIP'
fai='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa.fai'

histMarkList=['H3K27ac']
sampleList=[
'X4870_1NOSN',
'X4870_2NOSN',
'X4689NOSN',  
'X5628_1NOSN',
'X5628_2NOSN',
'X5130NOSN',
'X4996_1NOSN',
'X4996_2NOSN',
'X5006NOSN',
'X5244PDSN',
'X4831PDSN',
'X5532PDSN',
'X5215PDSN',
'X5627PDSN',
'X5778PDSN',
'X5742PDSN',
'X5649PDSN',
'X5591PDSN'
]

print len(sampleList)
os.system('mkdir '+DIR+'/CountMatrix')


def makeCountMatrix():
    for histMark in histMarkList:
	countDict={}
	for sampleID in sampleList:
	    input1=open(coverageDIR+'/'+sampleID+'.'+histMark+'.peakCoverage.txt','r')
	    all_input1=input1.readlines()
	    tmpDict={}
	    for line in all_input1:
		[chrID,pos1,pos2,cellTypes, readCount]=line.strip().split('\t')
		#[chrID,pos1,pos2, readCount]=line.strip().split('\t')
		peakID=chrID+':'+pos1+'-'+pos2	
		tmpDict[peakID]=readCount
	    countDict[sampleID]=tmpDict
	    input1.close()

	output1=open(DIR+'/CountMatrix/'+histMark+'.readCount.txt','w')
	output1.write('peakID'+'\t'+'\t'.join(sampleList)+'\n')
	for peakID in countDict[sampleList[0]]:
	    valList=[]
	    for sampleID in sampleList:
		val=countDict[sampleID][peakID]
		valList.append(val)
	    output1.write(peakID+'\t'+'\t'.join(valList)+'\n')
	output1.close()


def computeNoiseLevel():
    genomeSize=0
    input1=open(fai, 'r')
    all_input1=input1.readlines()
    for line in all_input1:
	each=line.strip().split('\t')
	size=float(each[1])
	genomeSize+=size
    input1.close()

    peakRegionSize=0
    input1=open(peakBED, 'r')
    all_input1=input1.readlines()
    for line in all_input1:
	each=line.strip().split('\t')
	size=float(each[2]) - float(each[1])
	peakRegionSize+=size
    input1.close()

    nonPeakRegion=genomeSize-peakRegionSize

    print 'genome size: ' +str(genomeSize)
    print 'peak region: '+str(peakRegionSize)+' ('+str(peakRegionSize/genomeSize*100)+'%)'
    print 'nonpeak region: '+str(nonPeakRegion)+' ('+str(nonPeakRegion/genomeSize*100)+'%)'


    noiseDict={}
    for histMark in histMarkList:
	input1=open(fripDIR+'/'+histMark+'.FRIP.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    sampleID=each[0]
	    nonPeakReads=float(each[3])-float(each[2])
	    noisePer100bp = nonPeakReads/nonPeakRegion*100

	    print histMark, sampleID
	    print 'noise per 100bp : ' +str(noisePer100bp)

	    noiseDict[histMark+';'+sampleID]=noisePer100bp

	input1.close()

    return noiseDict


def removeNoise():
    noiseDict=computeNoiseLevel()

    for histMark in histMarkList:
	valDict={}

	input1=open(DIR+'/CountMatrix/'+histMark+'.readCount.txt','r')
	all_input1=input1.readlines()
	sampleList=all_input1[0].strip().split('\t')[1:] 
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
	    peakID=each[0]
	    vals=numpy.array(each[1:], dtype='float')
	    tmpDict={}
	    index=0
	    for val in vals:
		tmpDict[sampleList[index]]=val
		index+=1
	    valDict[peakID]=tmpDict
	input1.close()

	output1=open(DIR+'/CountMatrix/'+histMark+'.readCount.denoised.txt','w')
	output1.write('peakID'+'\t'+'\t'.join(sampleList)+'\n')
	for peakID in valDict:
	    peakSize=float(peakID.split(':')[1].split('-')[1]) - float(peakID.split(':')[1].split('-')[0])
	    valList=[]
	    for sampleID in sampleList:
		val=valDict[peakID][sampleID]
		noiseLevel=noiseDict[histMark+';'+sampleID]
		noiseVal=peakSize/100*noiseLevel
		
		finalValue = val-noiseVal if not val-noiseVal < 0 else 0
		valList.append(str(finalValue))
	    output1.write(peakID+'\t'+'\t'.join(valList)+'\n')
	output1.close()



makeCountMatrix()
removeNoise()



