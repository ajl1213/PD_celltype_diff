#!/home/ajl1213/anaconda2/bin/python


import os
import numpy


DIR=os.getcwd()
allPairsEQTL_file='/home/ajl1213/Projects/PD/SciAdv_Anal/SNP_Mechanism/03_eQTL_enrich/0_Download/Brain_Substantia_nigra.allpairs.txt.gz'
GTF='/home/ajl1213/genome.info/gencode/hg19/gencode.v19.ptn_coding.lincRNA.level_1_2'

#targetGene='NUPL2'
#queryCoord='chr7:22800000-23250000'

#targetGene='GNA14'
#queryCoord='chr9:79975000-80680000'

targetGene='ITGB6'
queryCoord='chr2:160795000-161190000'


os.system('mkdir '+DIR+'/'+targetGene)


def loadGenes():
    geneDict={}
    input1=open(GTF,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	ensembleID=each[0]
	geneID=each[1]
	geneDict[geneID]=ensembleID
    input1.close()

    return geneDict


def extractGene():
    geneDict=loadGenes()

    if geneDict.has_key(targetGene):
	ensembleID=geneDict[targetGene]
	runLine='zcat '+allPairsEQTL_file+' | grep -w \''+ensembleID+'\' > '+DIR+'/'+targetGene+'/'+targetGene+'.all.txt'
	print runLine
	os.system(runLine)

    else:
	print 'matched gene not found'


def collectVal():
    extractGene()

    #q_chrID=queryCoord.split(':')[0]
    #q_pt1=queryCoord.split(':')[1].split('-')[0]
    #q_pt2=queryCoord.split(':')[1].split('-')[1]

    input1=open(DIR+'/'+targetGene+'/'+targetGene+'.all.txt','r')
    all_input1=input1.readlines()
    #output1=open(DIR+'/'+targetGene+'/'+q_chrID+'_'+q_pt1+'_'+q_pt2+'_'+targetGene+'.txt','w')
    output1=open(DIR+'/'+targetGene+'/'+targetGene+'.pval.txt','w')
    output1.write('genomicPt\tpval\n')
    for line in all_input1:
	each=line.strip().split('\t')
	pt=int(each[1].split('_')[1])
	pval=each[6]
	if float(pval) > 0:
	    #if (pt > int(q_pt1) and pt < int(q_pt2)):
	    newLine=[str(pt), str(pval)]
	    output1.write('\t'.join(newLine)+'\n')

    output1.close()
    input1.close()


def makePlot():
    q_chrID=queryCoord.split(':')[0]
    q_pt1=queryCoord.split(':')[1].split('-')[0]
    q_pt2=queryCoord.split(':')[1].split('-')[1]


    newLine='R CMD BATCH --no-save --no-restore \'--args '+targetGene+' '+q_chrID+' '+q_pt1+' '+q_pt2+'\' 01b_make_eQTL_track.R'
    print newLine
    os.system(newLine)


collectVal()
makePlot()



