#!/python2


import os
import gzip
import numpy
from scipy import stats


DIR=os.getcwd()
sigEQTL_file='/home/ajl1213/Projects/PD/SciAdv_Anal/SNP_Mechanism/03_eQTL_enrich/0_Download/GTEx_Analysis_v7_eQTL/Brain_Substantia_nigra.v7.signif_variant_gene_pairs.txt.gz'
totalPeak_DIR='/home/ajl1213/Projects/PD/data/snATAC/04_GetAtacCount/CellLabel'
fai='/home/ajl1213/genome.info/genome_index/hg19/hg19.hg19.fa.fai'
cellTypeList=['DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']


os.system('mkdir '+DIR+'/BED')


def makeBED():
    input1=gzip.open(sigEQTL_file,'r')
    output1=open(DIR+'/BED/BrainSN.sigEQTL.bed','w')

    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	varInfo=each[0]
	ensembleID=each[1]
	pval=each[6]

	chrID='chr'+varInfo.split('_')[0]
	pt=varInfo.split('_')[1]
	ref=varInfo.split('_')[2]
	alt=varInfo.split('_')[3]

	pt1=int(pt)
	pt2=int(pt)+1

	newLine=[chrID, str(pt1), str(pt2), ref, alt, ensembleID, pval]
	output1.write('\t'.join(newLine)+'\n')
	
    output1.close()    
    input1.close()


def intersect():
    for cellType in cellTypeList:
        run_line='intersectBed -a '+DIR+'/BED/BrainSN.sigEQTL.bed -b '+totalPeak_DIR+'/'+cellType+'.atacPeak.bed -wa -wb > '+DIR+'/BED/'+cellType+'.peakIntersect.bed'
        print run_line
        os.system(run_line)


def makeData():
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

    # peak size
    peakFracDict={}
    for cellType in cellTypeList:
	totalSize=0
	input1=open(totalPeak_DIR+'/'+cellType+'.atacPeak.bed','r')
	all_input1=input1.readlines()
	for line in all_input1:
	    each=line.strip().split('\t')
	    peakID=each[0]+':'+each[1]+'-'+each[2]
	    size=int(each[2])-int(each[1])
	    totalSize+=size
	input1.close()

	celltype_frac=float(totalSize)/float(genomeSize)
	peakFracDict[cellType]=celltype_frac

    # eQTL
    totalEqtlList=[]
    input1=open(DIR+'/BED/BrainSN.sigEQTL.bed','r')
    all_input1=input1.readlines()
    for line in all_input1:
	each=line.strip().split('\t')
	eqtlID=each[0]+':'+each[1]
	if not eqtlID in totalEqtlList:
	    totalEqtlList.append(eqtlID)
    input1.close()

    output1=open(DIR+'/sigEQTL.cRE_asso.txt','w')
    output1.write('CellType\tnEQTL\ttotalEqtl\texpected\tlog2Enrich\tPval\n')
    for cellType in cellTypeList:
        eqtlList=[]
        input1=open(DIR+'/BED/'+cellType+'.peakIntersect.bed','r')
        all_input1=input1.readlines()
        for line in all_input1:
            each=line.strip().split('\t')
            chrID=each[0]
            pt=each[1]
            eqtlID=chrID+':'+pt
	    eqtlList.append(eqtlID)

	    if not eqtlID in eqtlList:
		eqtlList.append(eqtlID)
        input1.close()

	#
	nEqtl=float(len(eqtlList))
	nExp=float(len(totalEqtlList)) * peakFracDict[cellType]

	log2Enrich=numpy.log2(nEqtl/nExp)
	#binomP=stats.binom_test(x=nEqtl, n=float(len(totalEqtlList)), p=peakFracDict[cellType])
	fisherOR=stats.fisher_exact([[nEqtl, float(len(totalEqtlList)-nEqtl)], [nExp, float(len(totalEqtlList)-nExp)]])[0]
	fisherP=stats.fisher_exact([[nEqtl, float(len(totalEqtlList)-nEqtl)], [nExp, float(len(totalEqtlList)-nExp)]])[1]
	print stats.fisher_exact([[nEqtl, float(len(totalEqtlList)-nEqtl)], [nExp, float(len(totalEqtlList)-nExp)]])
	newLine=[cellType, str(nEqtl), str(len(totalEqtlList)), str(nExp), str(log2Enrich), str(fisherP)]
	output1.write('\t'.join(newLine)+'\n')

    output1.close()
    




#makeBED()
#intersect()
makeData()




