#!/python2


import os
import numpy
from scipy import stats


DIR=os.getcwd()
matchedGWAS_file='/home/ajl1213/Projects/PD/SciAdv_Anal/SNP_Mechanism/01_SNP_Distrib/1_LdSNP/MatchedList.txt'
totalPeakBED='/home/ajl1213/Projects/PD/data/snATAC/04_GetAtacCount/CellLabel/TOTAL.atacPeak.bed'
dysPeak_DIR='/home/ajl1213/Projects/PD/NeuronJ/DiffAnal/02_H3K27ac/1_BulkH3K27ac/BedFiles'


os.system('mkdir '+DIR+'/BED')


def makeBED():
    input1=open(matchedGWAS_file, 'r')
    output1=open(DIR+'/BED/matchedGWAS.bed','w')
    all_input1=input1.readlines()
    for line in all_input1:
	each=line.strip().split('\t')
	sampleID=each[0]

	chrID=each[1].split('.')[0]
	pt=each[1].split('.')[1]

	pt1=int(pt)
	pt2=int(pt)+1

	if sampleID.count('PD') > 0:
	    output1.write(chrID+'\t'+str(pt1)+'\t'+str(pt2)+'\t'+sampleID+'\n')

    input1.close()
    output1.close()


def intersect():
    run_line='intersectBed -a '+DIR+'/BED/matchedGWAS.bed -b '+totalPeakBED+' -wa -wb > '+DIR+'/BED/matchedGWAS.peakIntersect.bed'
    print run_line
    os.system(run_line)


def loadDysPeak():
    dysPeakDict={}

    for i in ['Down','Up']:
        input1=open(dysPeak_DIR+'/'+i+'.Merged.bed','r')
        all_input1=input1.readlines()
        for line in all_input1:
            each=line.strip().split('\t')
            peakID=each[3]
            dysPeakDict[peakID]=i
        input1.close()

    return dysPeakDict


def makeData():
    dysPeakDict=loadDysPeak()

    # peak size
    peakSizeDict={}
    totalSize=0
    input1=open(totalPeakBED,'r')
    all_input1=input1.readlines()
    for line in all_input1:
        each=line.strip().split('\t')
        peakID=each[0]+':'+each[1]+'-'+each[2]
        size=int(each[2])-int(each[1])
        totalSize+=size

        if dysPeakDict.has_key(peakID):
            peakType=dysPeakDict[peakID]
        else:
            peakType='nonDys'

        if peakSizeDict.has_key(peakType):
            peakSizeDict[peakType]+=size
        else:
            peakSizeDict[peakType]=size
    input1.close()

    peakFracDict={}
    for peakType in peakSizeDict:
        peakFracDict[peakType] = float(peakSizeDict[peakType])/float(totalSize)


    # variant
    output1=open(DIR+'/matchedGWAS.cRE_asso.txt','w')
    output1.write('PeakType\tnVar\texpected\tlog2Enrich\tPval\n')
    varList=[]
    varDict={}
    input1=open(DIR+'/BED/matchedGWAS.peakIntersect.bed','r')
    all_input1=input1.readlines()
    for line in all_input1:
	each=line.strip().split('\t')
	chrID=each[0]
	pt=each[1]
	sampleID=each[3]
	varID=chrID+'.'+pt+';'+sampleID
	peakID=each[4]+':'+each[5]+'-'+each[6]

	if dysPeakDict.has_key(peakID):
	    peakType=dysPeakDict[peakID]
	else:
	    peakType='nonDys'

	if varDict.has_key(peakType):
	    tmp=varDict[peakType]
	    if not varID in tmp:
		tmp.append(varID)
	    varDict[peakType]=tmp
	else:
	    varDict[peakType]=[varID]

	if not varID in varList:
	    varList.append(varID)

    input1.close()

    # 
    for peakType in ['Down','Up','nonDys']:

	nVar=float(len(varDict[peakType]))
	nExp=float(len(varList)) * peakFracDict[peakType]

	log2Enrich=numpy.log2(nVar/nExp)
	binomP = stats.binom_test(x=nVar, n=float(len(varList)), p=peakFracDict[peakType])
	newLine=[peakType, str(nVar), str(nExp), str(log2Enrich), str(binomP)]
	#print peakType, nVar, nExp, log2Enrich, binomP
	output1.write('\t'.join(newLine)+'\n')

    output1.close()




makeBED()
intersect()
makeData()




