#!/home/ajl1213/anaconda2/bin/python


import os


DIR=os.getcwd()
gwasList=['PD_Nalls_2014', 'PD_Chang_2017', 'PD_Foo_2017', 'PD_Nalls_2019']
rawGwasDIR='/home/ajl1213/Projects/PD/SciAdv_Anal/GWAS_Enrich/1_CollectSNP/rawGWAS'


def getSnpInfo():
    snpDict={}
    for gwasID in gwasList:
	input1=open(rawGwasDIR+'/'+gwasID+'.rawGWAS.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
	    [rsID, chrID, pos, pval]=line.strip().split('\t')
	    snpPt=chrID+'.'+pos
	    valID=gwasID+';'+pval
	    if snpDict.has_key(snpPt):
		tmp=snpDict[snpPt]
		tmp.append(valID)
		snpDict[snpPt]=tmp
	    else:
		snpDict[snpPt]=[valID]
	input1.close()

    return snpDict


def makeSnpTable():
    runScript='cat RESULT/MergedPD.peakIntersect.txt | grep -w -v \'LdSnpPos\' | cut -f2,3,4,5,6 | sort -u > '+DIR+'/SnpTable/MergedPD.SnpTable.txt'
    print runScript
    os.system(runScript)

    runScript='cat RESULT/MergedPD.peakIntersect.sigPeak.txt | grep -w -v \'LdSnpPos\' | cut -f2,3,4,5,6 | sort -u > '+DIR+'/SnpTable/MergedPD.SigPeak.SnpTable.txt'
    print runScript
    os.system(runScript)

    peakSnpList=[]
    peakDict={}
    input1=open(DIR+'/SnpTable/MergedPD.SnpTable.txt','r')
    all_input1=input1.readlines()
    for line in all_input1:
	each=line.strip().split('\t')
	snpPt=each[0]
	rsID=each[1]
	peakID=each[2]
	peakLabel=each[3]
	peakSnpList.append(snpPt)
	if peakDict.has_key(snpPt):
	    tmp=peakDict[snpPt]
	    tmp.append(peakLabel)
	    peakDict[snpPt]=tmp
	else:
	    peakDict[snpPt]=[peakLabel]
    input1.close()

    return peakSnpList, peakDict


def makeLabel():
    snpDict = getSnpInfo()
    peakSnpList, peakDict = makeSnpTable()

    input1=open(DIR+'/SnpTable/MergedPD.SigPeak.SnpTable.txt','r')
    output1=open(DIR+'/SnpTable/MergedPD.SigPeak.SnpTable.label.txt','w')
    output1.write('ParentalSnpPos\tParentalSnpID\tGwasType\tPval\tpeakID\tpeakLabel\tmaxCell\n')
    all_input1=input1.readlines()
    for line in all_input1:
	[ParentalSnpPos, ParentalSnpID, peakID, peakLabel, maxCell]=line.strip().split('\t')

	gwasList=[]
	pvalList=[]
	for i in snpDict[ParentalSnpPos]:
	    a=i.split(';')[0]
	    b=i.split(';')[1]

	    gwasList.append(a)
	    pvalList.append(b)
	
	gwasID=','.join(gwasList)
	pval=','.join(pvalList)

	newLine=[ParentalSnpPos, ParentalSnpID, gwasID, pval, peakID, peakLabel, maxCell]
	output1.write('\t'.join(newLine)+'\n')
    output1.close()

    output1=open(DIR+'/SnpTable/MergedPD.PeakAsso.txt','w')
    output1.write('SnpID\tpeakLabel\tpeakType\n')
    for snpPt in snpDict:
	if snpPt in peakSnpList:
	    peakLabel='peakAsso'
	else:
	    peakLabel='none'

	if peakDict.has_key(snpPt):
	    if 'dysPeak' in peakDict[snpPt]:
		peakType='dysPeak'
	    else:
		peakType='nonDysPeak'
	else:
	    peakType='none'

	newLine=[snpPt, peakLabel, peakType]
	output1.write('\t'.join(newLine)+'\n')
    output1.close()



makeLabel()



