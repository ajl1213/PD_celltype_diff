#!/python2


import os


DIR=os.getcwd()
inputFile=DIR+'/AllelicReads.GWAS_matched.pval.txt'
dysPeak_DIR='/home/ajl1213/Projects/PD/SciAdv_Anal/Diff_cRE/1_BulkH3K27ac/BedFiles'
sampleList=[
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

os.system('mkdir '+DIR+'/BED.v2')


for i in sampleList:
    tmpList=[]
    input1=open(inputFile,'r')
    all_input1=input1.readlines()
    for line in all_input1[1:]:
	each=line.strip().split('\t')
	queryID=each[0]

	sampleID=queryID.split(';')[0]
	chrID=queryID.split(';')[1]
	pt=queryID.split(';')[2]

	label=each[11]

	if i==sampleID:
	    if label=='refEnrich':
		snpID=chrID+':'+pt
		tmpList.append(snpID)
    input1.close()

    if len(tmpList) > 0:
	output1=open(DIR+'/BED.v2/'+i+'.GWAS_matched.bed','w')
	for snpID in tmpList:
	    chrID=snpID.split(':')[0]
	    pt=int(snpID.split(':')[1])

	    pt1=pt
	    pt2=pt+1
	    newLine=[chrID, str(pt1), str(pt2)]
	    output1.write('\t'.join(newLine)+'\n')
	output1.close()

	##
	runLine='intersectBed -a '+DIR+'/BED.v2/'+i+'.GWAS_matched.bed -b '+dysPeak_DIR+'/Down.Merged.bed -wa -wb | cut -f1,2,3 | sort -u > '+DIR+'/BED.v2/'+i+'.DownPeak.intersect.bed'
        print runLine
        os.system(runLine)

        #   
        tmpList=[]
        input1=open(DIR+'/BED.v2/'+i+'.DownPeak.intersect.bed','r')
        all_input1=input1.readlines()
        for line in all_input1:
            each=line.strip().split('\t')
            chrID=each[0]
            pt1=each[1]
            pt2=each[2]
            snpPt=chrID+':'+pt1+'-'+pt2
            tmpList.append(snpPt)
        input1.close()

        if len(tmpList) > 0:
            output1=open(DIR+'/BED.v2/'+i+'.GWAS_matched.DownPeak.final.bed','w')
            for i in tmpList:
                chrID=i.split(':')[0]
                pt1=i.split(':')[1].split('-')[0]
                pt2=i.split(':')[1].split('-')[1]
                newLine=[chrID, pt1, pt2]
                output1.write('\t'.join(newLine)+'\n')
            output1.close()



