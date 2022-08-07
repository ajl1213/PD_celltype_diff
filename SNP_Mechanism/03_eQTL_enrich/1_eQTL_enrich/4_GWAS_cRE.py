#!/python2


import os
import numpy
from scipy import stats


DIR=os.getcwd()
GWAS_file='/home/ajl1213/Projects/PD/SciAdv_Anal/GWAS_Enrich/2_ClassifySNP/BedFiles/MergedPD.bed'
totalPeak_file='/home/ajl1213/Projects/PD/data/snATAC/04_GetAtacCount/CellLabel/AtacPeak.TOTAL.bed'
sigEQTL_file=DIR+'/BED/BrainSN.sigEQTL.bed'
fai='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa.fai'


def loadGWAS(): 
    runLine='intersectBed -a '+GWAS_file+' -b '+totalPeak_file+' -wa -wb > '+DIR+'/BED/gwas_totalPeak.intersect.bed'
    print runLine
    os.system(runLine)

    gwasDict={}
    input1=open(DIR+'/BED/gwas_totalPeak.intersect.bed','r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	gwasPos=each[0]+':'+each[1]
	if not gwasDict.has_key(gwasPos):
	    gwasDict[gwasPos]=1
    input1.close()

    return gwasDict


def loadEQTL():
    runLine='intersectBed -a '+sigEQTL_file+' -b '+totalPeak_file+' -wa -wb > '+DIR+'/BED/eqtl_totalPeak.intersect.bed'
    print runLine
    os.system(runLine)

    eqtlDict={}
    input1=open(DIR+'/BED/eqtl_totalPeak.intersect.bed','r')
    all_input1=input1.readlines()
    for line in all_input1:
	each=line.strip().split('\t')
	eqtlPos=each[0]+':'+each[1]
	if not eqtlDict.has_key(eqtlPos):
	    eqtlDict[eqtlPos]=1
    input1.close()

    return eqtlDict


def getMatch():
    gwasDict=loadGWAS()
    eqtlDict=loadEQTL()

    commonList=set(gwasDict) & set(eqtlDict)

    ##
    genomeSize=0
    input1=open(fai, 'r')
    all_input1=input1.readlines()
    for line in all_input1:
        each=line.strip().split('\t')
        chrID=each[0]
        size=int(each[1])

        if not (chrID=='chrY' or chrID=='chrM'):
            genomeSize+=size
    input1.close()


    ##
    output1=open(DIR+'/GWAS_eQTL.asso.txt','w')
    output1.write('dataType\toverlapGWAS\ttotalGWAS\texpected\tlog2Enrich\tPval\n')

    overlapGWAS = len(commonList)
    totalGWAS = len(gwasDict)

    expected_frac = float(len(eqtlDict)) / float(genomeSize)
    expected_count = float(len(gwasDict)) * expected_frac
    log2Enrich = numpy.log2(float(overlapGWAS)/expected_count)
    fisherP = stats.fisher_exact([[float(overlapGWAS), float(totalGWAS-overlapGWAS)], [expected_count, float(totalGWAS-expected_count)]])[1]

    newLine=['GWAS', str(overlapGWAS), str(totalGWAS), str(expected_count), str(log2Enrich), str(fisherP)]
    output1.write('\t'.join(newLine)+'\n')
    output1.close()




getMatch()



