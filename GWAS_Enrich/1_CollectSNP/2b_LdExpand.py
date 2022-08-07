#!/home/ajl1213/anaconda2/bin/python



import os
import fileinput
import sys



DIR=os.getcwd()
inputDIR=sys.argv[1]
LD_DIR=sys.argv[2]
gwasID=sys.argv[3]
chrID=sys.argv[4]


def load_LD(chrID,targetDict):
    LD_block={}
    print('loading LD')
    for i in ['AFR','AMR','EAS','EUR','SAS']:
	print(i)
	input1=open(LD_DIR+'/'+chrID+'_'+i+'_ld_all.hap.ld','rt')
	for line in input1:
	    if not line.startswith('CHR'):
		each=line.split()
		pos1=each[1]
		pos2=each[2]
		r_sqr=each[4]

		if float(r_sqr) > 0.8:
		    if targetDict.has_key(pos1):
			if LD_block.has_key(pos1):
			    tmp=LD_block[pos1]
			    tmp.append(pos2)
			    LD_block[pos1]=tmp
			else:
			    LD_block[pos1]=[pos2]

		    if targetDict.has_key(pos2):
			if LD_block.has_key(pos2):
			    tmp=LD_block[pos2]
			    tmp.append(pos1)
			    LD_block[pos2]=tmp
			else:
			    LD_block[pos2]=[pos1]
	input1.close()

    return LD_block


def expandData(gwasID, chrID):
    input1=open(inputDIR+'/'+gwasID+'.rawGWAS.txt','rt')
    all_input1=input1.readlines()
    targetDict={}
    targetList=[]
    for i in all_input1[1:]:
	each=i.strip().split('\t')
	chr=each[1]
	if chr==chrID:
	    pos=each[2]
	    snpID=each[0]
	    #pval=each[3]
	    targetDict[pos]=[snpID, chr]
	    targetList.append([snpID, chr, pos])
    input1.close()

    LD_block=load_LD(chrID, targetDict)

    output1=open(DIR+'/LdExpand/'+gwasID+'.'+chrID+'.LD.txt','w')
    for i in targetList:
	snpID=i[0]
	chr=i[1]

	if chr==chrID:
	    snpPos=i[2]
	    if snpPos in LD_block:
		output1.write(snpID+'\t'+chr+'\t'+snpPos+'\t')
		tmp_ld_list=LD_block[snpPos]
		for j in set(tmp_ld_list):
		    output1.write(j+':'+str(tmp_ld_list.count(j))+',')
		output1.write('\n')
	    else:
		output1.write(snpID+'\t'+chr+'\t'+snpPos+'\tno\n')
    output1.close()




expandData(gwasID, chrID)



