#!/python


import os


DIR=os.getcwd()
VCF_DIR='/home/ajl1213/Projects/PD/SciAdv_Data/ChIP/03_VariantCall'
FASTA='/home/ajl2/genome.info/refGenome/hg19/hg19.hg19.fa'

sampleList=[
'X4870_1NOSN',
'X4870_2NOSN',
'X4689NOSN',
'X5628_1NOSN',
'X5628_2NOSN',
'X5130NOSN',
'X4996_1NOSN',
'X4996_2NOSN',
'X5006NOSN',
'X5244PDSN',
'X4831PDSN',
'X5532PDSN',
'X5215PDSN',
'X5627PDSN',
'X5778PDSN',
'X5742PDSN',
'X5649PDSN',
'X5591PDSN'
]


os.system('mkdir '+DIR+'/VCF')


def processVCF():

        print('step 1 / step 2')
        for sampleID in sampleList:
                print(sampleID)

                ## command for downlading appropriate database into humandb directory
                # annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/

                ## Process vcf file according to ANNOVAR suggestion
                ## Step1: split vcf lines so that each line contains one and only one variant (originall vcf contains multiple variants in a single line)
                ## Step2: peform left-normalization (A way to process indel)

                run_line='bcftools norm -m-both -o '+DIR+'/VCF/'+sampleID+'.step1.vcf '+VCF_DIR+'/'+sampleID+'.vcf'
                print(run_line)
                os.system(run_line)

                run_line='bcftools norm -f '+FASTA+' -o '+DIR+'/VCF/'+sampleID+'.step2.vcf '+DIR+'/VCF/'+sampleID+'.step1.vcf'
                print(run_line)
                os.system(run_line)


        ## Step3: conduct filter
        ## Since its "Minor" allele frequency, remove variants with more than 50% occurrence.
        print('step 3')
        allDict={}
        for sampleID in sampleList:
                print(sampleID)
                input1=open(DIR+'/VCF/'+sampleID+'.step2.vcf','r')
                all_input1=input1.readlines()
                for line in all_input1:
                        if not line.count('#') > 0:
                                each=line.strip().split('\t')

                                chrID=each[0]
                                pt=each[1]
                                ref=each[3]
                                alt=each[4]

                                keyID=chrID+';'+pt+';'+ref+';'+alt
                                if keyID in allDict:
                                        allDict[keyID]+=1
                                else:
                                        allDict[keyID]=1
    
                input1.close()

        for sampleID in sampleList:
                input1=open(DIR+'/VCF/'+sampleID+'.step2.vcf','r')
                output1=open(DIR+'/VCF/'+sampleID+'.step3.vcf','w')
                all_input1=input1.readlines()
                for line in all_input1:

                        if not line.count('#') > 0:
                                each=line.strip().split('\t')

                                chrID=each[0]
                                pt=each[1]
                                ref=each[3]
                                alt=each[4]

                                keyID=chrID+';'+pt+';'+ref+';'+alt
                                if not allDict[keyID] > len(sampleList)/2:
                                        output1.write('\t'.join(each)+'\n')

                        else:
                                output1.write(line)

                input1.close()
                output1.close()




processVCF()





