#!/home/ajl1213/anaconda2/bin/python


import os


DIR=os.getcwd()
PythonDIR='/home/ajl1213/anaconda2/bin'
ScriptDIR='/home/ajl1213/programs/ldsc'
GwasDIR='/home/ajl1213/genome.info/GWAS_sumstat'

# ctrl / case
nCaseDict={
'PD_Nalls_2014':[95282,13708],
'PD_Chang_2017':[403190,26035],
'PD_Foo_2017':[17604,5125],
'PD_Nalls_2019':[1400000,56306],
'AD_Jansen_2019':[455258],
'AD_Kunkle_2019':[21982,41944],
'ALS_VanRheenen_2016':[23475,12577],
'ASD_Grove_2019':[27969,18382],
'SCZ_Pardinas_2018':[64643,40675]
}

os.system('mkdir '+DIR+'/GwasSummaryStat')



def collectInfo():
    # Download baseline LD scores
    runScript='wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_baseline_ldscores.tgz'
    print runScript
    os.system(runScript)

    runScript='tar -xvzf 1000G_Phase3_baseline_ldscores.tgz'
    print runScript
    os.system(runScript)

    # Download SNP weight
    runScript='wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/weights_hm3_no_hla.tgz'
    print runScript
    os.system(runScript)

    runScript='tar -xvzf weights_hm3_no_hla.tgz'
    print runScript
    os.system(runScript)

    # Download HapMap3 SNP list
    runScript='wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/w_hm3.snplist.bz2'
    print runScript
    os.system(runScript)

    runScript='bunzip2 w_hm3.snplist.bz2'
    print runScript
    os.system(runScript)


def collectGWAS():

    for GwasID in nCaseDict:

	if len(nCaseDict[GwasID])==1:

	    nCase=str(nCaseDict[GwasID][0])
	    script=PythonDIR+'/python '+ScriptDIR+'/munge_sumstats.py --N '+str(nCase)+' --sumstats '+GwasDIR+'/'+GwasID+'/sumStat.txt --merge-alleles w_hm3.snplist --out '+DIR+'/GwasSummaryStat/'+GwasID
	    print script
	    os.system(script)


	if len(nCaseDict[GwasID])==2:
	    nCtrl=str(nCaseDict[GwasID][0])
	    nCase=str(nCaseDict[GwasID][1])

	    script=PythonDIR+'/python '+ScriptDIR+'/munge_sumstats.py --N-con '+str(nCtrl)+' --N-cas '+str(nCase)+' --a1-inc --sumstats '+GwasDIR+'/'+GwasID+'/sumStat.txt --merge-alleles w_hm3.snplist --out '+DIR+'/GwasSummaryStat/'+GwasID
	    print script
	    os.system(script)



collectInfo()
collectGWAS()



