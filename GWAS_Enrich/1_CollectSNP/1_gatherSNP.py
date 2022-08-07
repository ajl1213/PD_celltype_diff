#!/home/ajl1213/anaconda2/bin/python



import os



DIR=os.getcwd()
gwasInputDIR='/home/ajl1213/genome.info/GWAS_sumstat'
gwasList= ['PD_Nalls_2014', 'PD_Chang_2017', 'PD_Foo_2017', 'PD_Nalls_2019']
newID='MergedPD'
PvalCut=5e-8

os.system('mkdir '+DIR+'/rawGWAS')


def collectSNP():

    allDict={}
    for gwasID in gwasList:
	input1=open(gwasInputDIR+'/'+gwasID+'/sumStat.txt','r')
	all_input1=input1.readlines()
	colnames=all_input1[0].strip().split('\t')

	snp_idx=colnames.index('SNP')
	chr_idx=colnames.index('CHR')
	pos_idx=colnames.index('POS')
	pval_idx=colnames.index('P')

	output1=open(DIR+'/rawGWAS/'+gwasID+'.rawGWAS.txt','w')
	output1.write('rsID\tchrID\tpos\tpval\n')
	for line in all_input1[1:]:
	    each=line.strip().split('\t')
            chrID='chr'+each[chr_idx]
            pt=each[pos_idx]
            rsID=each[snp_idx]
            pval=each[pval_idx]

	    if float(pval) < PvalCut:
		newLine=[rsID, chrID, pt, pval]
		output1.write('\t'.join(newLine)+'\n')

		#
		keyID=chrID+'.'+pt
		if allDict.has_key(keyID):
		    tmp=allDict[keyID]
		    tmp.append(rsID)
		    allDict[keyID]=tmp
		else:
		    allDict[keyID]=[rsID]
	input1.close()

    return allDict


def mergeSNP():
    allDict=collectSNP()

    output1=open(DIR+'/rawGWAS/'+newID+'.rawGWAS.txt','w')
    output1.write('rsID\tchrID\tpos\n')
    for keyID in allDict:

	chrID=keyID.split('.')[0]
	pt=keyID.split('.')[1]
	if len(allDict[keyID])>1:
	    tmpList=[]
	    for i in allDict[keyID]:
		rsID=i.split(';')[0]
		tmpList.append(rsID)
	    rsID=';'.join(set(tmpList))

	else:
	    rsID=allDict[keyID][0]
	newLine=[rsID, chrID, pt]
	output1.write('\t'.join(newLine)+'\n')
    output1.close()




collectSNP()
mergeSNP()



