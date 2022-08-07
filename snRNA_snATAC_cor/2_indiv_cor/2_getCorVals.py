#!/python2


import os


DIR=os.getcwd()


cellTypeList=['DopaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']
atacPooledLibs=[
'X4996NOSN',
'X5130NOSN',
'X5627PDSN',
'X5532PDSN',
'X5591PDSN',
'X5244PDSN'
]


output1=open(DIR+'/Merged.cor.label.txt','w')
output1.write('cellType\trnaSample\tatacSample\tlabel\tatacLibType\tpcc\n')
for cellType in cellTypeList:

    valDict={}
    input1=open(DIR+'/CorMat/'+cellType+'.corMat.txt','r')
    all_input1=input1.readlines()
    atacSamples=all_input1[0].strip().split('\t')[1:]

    for line in all_input1[1:]:
	each=line.strip().split('\t')
	rnaSample=each[0]
	vals=each[1:]

	idx=0
	tmpDict={}
	for i in vals:
	    tmpDict[atacSamples[idx]]=i
	    idx+=1
	valDict[rnaSample]=tmpDict

    input1.close()

    for rnaSample in valDict:
	for atacSample in atacSamples:

	    if rnaSample.split('_')[0] == atacSample.split('_')[0]:
		label='match'
	    else:
		label='unmatch'

	    if atacSample.split('_')[0] in atacPooledLibs:
		libType='pooledLib'
	    else:
		libType='singleLib'

	    val=valDict[rnaSample][atacSample]

	    newLine=[cellType, rnaSample, atacSample, label, libType, val]
	    output1.write('\t'.join(newLine)+'\n')
output1.close()




