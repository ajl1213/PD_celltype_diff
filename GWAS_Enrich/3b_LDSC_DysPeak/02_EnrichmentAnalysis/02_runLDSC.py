#!/home/ajl1213/anaconda2/bin/python


import os
import numpy


DIR=os.getcwd()
pythonDIR='/home/ajl1213/anaconda2/bin'
scriptDIR='/home/ajl1213/programs/ldsc'
LdScoreDIR='/home/ajl1213/Projects/PD/NeuronJ/GwasEnrich/3b_LDSC_DysPeak/01_PartitionedLdScores/AnnoLdFile'
GwasList=[
'PD_Nalls_2014',
'PD_Chang_2017',
'PD_Foo_2017',
'PD_Nalls_2019',
'AD_Jansen_2019',
'AD_Kunkle_2019',
'ALS_VanRheenen_2016',
'ASD_Grove_2019',
'SCZ_Pardinas_2018'
]

dataTypeList=['Down','Up']
cellTypeList=['Merged', 'DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']


os.system('mkdir '+DIR+'/RES_LDSC')
os.system('mkdir '+DIR+'/ldcts')


def make_LDref():
    for GwasID in GwasList:
	output1=open(DIR+'/ldcts/'+GwasID+'.ldcts', 'w')
	for dataType in dataTypeList:
	    for cellType in cellTypeList:
		output1.write(dataType+'_'+cellType+'\t'+LdScoreDIR+'/'+dataType+'.'+cellType+'.\n')
	output1.close()


def runRegression():
    make_LDref()

    for GwasID in GwasList: 
	runScript=pythonDIR+'/python '+scriptDIR+'/ldsc.py --h2-cts '+DIR+'/GwasSummaryStat/'+GwasID+'.sumstats.gz --ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. --out '+DIR+'/RES_LDSC/'+GwasID+' --ref-ld-chr-cts '+DIR+'/ldcts/'+GwasID+'.ldcts --w-ld-chr weights_hm3_no_hla/weights.'
	print runScript
	os.system(runScript)


def makePvalMat():
    pvalDict={}
    for GwasID in GwasList:
        input1=open(DIR+'/RES_LDSC/'+GwasID+'.cell_type_results.txt','r')
        all_input1=input1.readlines()
        tmpDict={}
        for line in all_input1[1:]:
            each=line.strip().split('\t')
            keyID=each[0]
            pval=float(each[3])
            tmpDict[keyID]=pval

        pvalDict[GwasID]=tmpDict
        input1.close

    keyList=[]
    for dataType in dataTypeList:
	for cellType in cellTypeList:
	    keyID=dataType+'_'+cellType
	    keyList.append(keyID)

    ## pval
    output1=open(DIR+'/pvalList.txt','w')
    output1.write('GwasID\tdataType\tpval\n')
    for GwasID in GwasList:
        for keyID in keyList:
            pval=str(pvalDict[GwasID][keyID])
	    newLine=[GwasID, keyID, pval]
	    output1.write('\t'.join(newLine)+'\n')
    output1.close()

    ## logPval
    output1=open(DIR+'/logPvalMat.txt','w')
    output1.write('GwasID\t'+'\t'.join(keyList)+'\n')
    for GwasID in GwasList:
        valList=[]
        for keyID in keyList:
            logPval=str(-numpy.log10(pvalDict[GwasID][keyID]))
            valList.append(logPval)
        output1.write(GwasID+'\t'+'\t'.join(valList)+'\n')
    output1.close()





#runRegression()
makePvalMat()




