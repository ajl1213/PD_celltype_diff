#!/home/ajl1213/anaconda2/bin/python


import os 
import gzip



DIR=os.getcwd()
sampleList=['MergedNOSN','MergedPDSN']
chrList=['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',  'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16','chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
resolution='5kb'

#qvalCut=[0.01,'1e-2']
qvalCut=[0.001,'1e-3']


def mergeInter():
    for sampleID in sampleList:
        print sampleID
        output1=open(DIR+'/OUTPUT/'+sampleID+'.'+resolution+'.'+str(qvalCut[1])+'.SigInter','w')
        output1.write('frag1Chr\tfrag1Mid\tfrag2Chr\tfrag2Mid\tcontactCount\tp_value\tq_value\tbias1\tbias2\n')
        i=1
        while i<=23:
            if i==23:
                chrID='chrX'
            else:
                chrID='chr'+str(i)

            print chrID
            input1=gzip.open(DIR+'/OUTPUT/'+sampleID+'.'+chrID+'.'+resolution+'/FitHiC.spline_pass2.significances.txt.gz','r')
            all_input1=input1.readlines()
            for j in all_input1[1:]:
                each=j.split()
                if float(each[6]) < qvalCut[0]:
                    output1.write(j)
            i+=1
        output1.close()






mergeInter()



