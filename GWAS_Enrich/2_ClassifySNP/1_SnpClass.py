#!/home/ajl1213/anaconda2/bin/python


import os



DIR=os.getcwd()
inputGwasFile='/home/ajl1213/Projects/PD/SciAdv_Anal/GWAS_Enrich/1_CollectSNP/LdFinal/MergedPD.allchr.txt'
allPeakBED='/home/ajl1213/Projects/PD/data/snATAC/04_GetAtacCount/CellLabel/AtacPeak.TOTAL.bed'
diffPeakFile='/home/ajl1213/Projects/PD/SciAdv_Data/ChIP/02_Processing/ValTable/H3K27ac.valTable.filter.txt'
gwasID='MergedPD'

os.system('mkdir '+DIR+'/BedFiles/')
os.system('mkdir '+DIR+'/RESULT/')
os.system('mkdir '+DIR+'/SnpTable/')



def makeBedFiles():
    ##
    input1=open(inputGwasFile,'r')
    output1=open(DIR+'/BedFiles/'+gwasID+'.bed','w')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	[LdSnpPos, ParentalSnpPos, ParentalSnpID]=line.strip().split('\t') 
	chrID=LdSnpPos.split('.')[0]
	pt1=LdSnpPos.split('.')[1]
	pt2=str(int(pt1)+1)
	newLine=[chrID, pt1, pt2, ParentalSnpPos, ParentalSnpID]
	output1.write('\t'.join(newLine)+'\n')
    input1.close()
    output1.close()


def loadSigPeakInfo():
    sigPeakDict={}
    input1=open(diffPeakFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	peakID=each[0]
	label=each[7]

        sigPeakDict[peakID]=label
    input1.close()

    return sigPeakDict


def labelGwasSNP():
    makeBedFiles()
    sigPeakDict=loadSigPeakInfo()

    runLine='intersectBed -a '+allPeakBED+' -b '+DIR+'/BedFiles/'+gwasID+'.bed -wa -wb | sort -u > '+DIR+'/BedFiles/'+gwasID+'.GWAS.AllPeak.intersect.bed'
    os.system(runLine)

    input1=open(DIR+'/BedFiles/'+gwasID+'.GWAS.AllPeak.intersect.bed','r')
    output1=open(DIR+'/RESULT/'+gwasID+'.peakIntersect.txt','w')
    output1.write('LdSnpPos\tParentalSnpPos\tParentalSnpID\tpeakID\tpeakLabel\tmaxCell\n')
    all_input1=input1.readlines()
    for line in all_input1:
	each=line.strip().split('\t')
	peakID=each[0]+':'+each[1]+'-'+each[2]
	if sigPeakDict.has_key(peakID):
	    peakLabel='dysPeak'
	else:
	    peakLabel='nonDysPeak'

	maxCell=each[3]
	LdSnpPos=each[4]+'.'+each[5]
	ParentalSnpPos=each[7]
	ParentalSnpID=each[8]
	newLine=[LdSnpPos, ParentalSnpPos, ParentalSnpID, peakID, peakLabel, maxCell]
	output1.write('\t'.join(newLine)+'\n')
    input1.close()
    output1.close()

    ## SigPeak
    input1=open(DIR+'/BedFiles/'+gwasID+'.GWAS.AllPeak.intersect.bed','r')
    output1=open(DIR+'/RESULT/'+gwasID+'.peakIntersect.sigPeak.txt','w')
    output1.write('LdSnpPos\tParentalSnpPos\tParentalSnpID\tpeakID\tpeakLabel\tmaxCell\n')
    all_input1=input1.readlines()
    for line in all_input1:
	each=line.strip().split('\t')
	peakID=each[0]+':'+each[1]+'-'+each[2]
	if sigPeakDict.has_key(peakID):
	    peakLabel='dysPeak'
	else:
	    peakLabel='nonDysPeak'
	maxCell=each[3]
	LdSnpPos=each[4]+'.'+each[5]
	ParentalSnpPos=each[7]
	ParentalSnpID=each[8]

	if sigPeakDict.has_key(peakID):
	    peakLabel=sigPeakDict[peakID]

	    newLine=[LdSnpPos, ParentalSnpPos, ParentalSnpID, peakID, peakLabel, maxCell]
	    output1.write('\t'.join(newLine)+'\n')
    input1.close()
    output1.close()




labelGwasSNP()



