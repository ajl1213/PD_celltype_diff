#!/home/ajl1213/anaconda2/bin/python


import os


DIR=os.getcwd()
inputDIR=DIR+'/LdExpand'
gwasList=['MergedPD', 'PD_Nalls_2014', 'PD_Chang_2017', 'PD_Foo_2017', 'PD_Nalls_2019']

os.system('mkdir '+DIR+'/LdFinal')


def PostProcess():

    for gwasID in gwasList:
	snpDict={}
	for i in range(23):
	    if i==22:
		chrID='chrX'
	    else:
		chrID='chr'+str(i+1)

	    input1=open(DIR+'/LdExpand/'+gwasID+'.'+chrID+'.LD.txt','r')
	    all_input1=input1.readlines()
	    for line in all_input1:
		each=line.strip().split('\t')
		rsID=each[0]
		chrID=each[1]
		parentalPos=each[2]
		ldSnpList=each[3]
		snpDict[chrID+'.'+parentalPos]=[rsID+'.'+chrID+'.'+parentalPos]

		if not ldSnpList=='no':
		    for i in ldSnpList.split(',')[:-1]:
			if int(i.split(':')[1]) > 2:   ### Either 2, 3, 4 or 5
			    ldSnp=chrID+'.'+i.split(':')[0]
			    if snpDict.has_key(ldSnp):
				tmp=snpDict[ldSnp]
				tmp.append(rsID+'.'+chrID+'.'+parentalPos)
				snpDict[ldSnp]=tmp
			    else:
				snpDict[ldSnp]=[rsID+'.'+chrID+'.'+parentalPos]
	    input1.close()
		    
	output1=open(DIR+'/LdFinal/'+gwasID+'.allchr.txt','w')
	output1.write('LdSnpPos\tParentalSnpPos\tParantalSnpID\n')
	for snpID in snpDict:
	    for i in snpDict[snpID]:
		chrID=i.split('.')[1]
		parentalSnpID=i.split('.')[0]
		parentalSnpPos=i.split('.')[2]
		parentalSNP=chrID+'.'+parentalSnpPos

		newLine=[snpID, parentalSNP, parentalSnpID]
		output1.write('\t'.join(newLine)+'\n')
	output1.close()
	    


PostProcess()


