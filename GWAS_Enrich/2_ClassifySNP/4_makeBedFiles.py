#!/home/ajl1213/anaconda2/bin/python


import os




#=====
# target peak: chr7:23144778-23146888
# target snp: rs2072368, rs2072369

#=====
# target peak: chr4:90750339-90751275
#

snpList=[
'rs10008964',
'rs10019540','rs1372516',
'rs11188877','rs2583968',
'rs11447716','rs983361',
'rs11946349','rs2737004',
'rs12528978','rs2583965',
'rs12616630','rs2736988',
'rs12910772','rs2737024',
'rs12947890','rs2619366',
'rs13142587','rs4285703',
'rs13590404','rs2245801',
'rs13870302','rs2737012',
'rs14339790','rs2619370',
'rs1442144','rs4797451',
'rs1442146','rs5539728',
'rs1442147','rs16796700',
'rs1442148','rs9560798',
'rs15101855',
'rs15116944','rs2737034',
'rs15158283','rs2583971',
'rs15481903','rs1442149',
'rs15582908','rs2572320',
'rs15848575','rs2583989',
'rs16828787','rs2619360',
'rs1837891','rs6527334',
'rs2028535','rs5114204',
'rs2197120','rs14073149',
'rs2572318','rs12693090',
'rs2583958','rs10084431',
'rs2583959','rs7240758',
'rs2583969','rs2905582',
'rs2583973','rs1642125',
'rs2583977','rs14949701',
'rs2583978','rs13046577',
'rs2583979','rs12837215',
'rs2583981','rs5495230',
'rs2583983','rs8834808',
'rs2619339','rs11121980',
'rs2619341','rs8079234',
'rs2619342','rs5495060',
'rs2619343','rs7627598',
'rs2619346',
'rs2619347','rs3626430',
'rs2619349','rs4749813',
'rs2619351','rs6198233',
'rs2619352','rs3042787',
'rs2619354','rs2607869',
'rs2619355','rs3284580',
'rs2619356','rs11520957',
'rs2619359','rs4521184',
'rs2619364','rs1252228',
'rs2619373','rs1529860',
'rs2736996','rs10594219',
'rs2736997','rs7526933',
'rs2737001','rs12172894',
'rs2737002','rs15589307',
'rs2737003','rs1680809',
'rs2737008','rs1017621',
'rs2737014','rs11205752',
'rs2737015','rs15325353',
'rs2737021','rs1698756',
'rs2737025','rs16233739',
'rs2974588','rs2619362',
'rs3055063','rs2736995',
'rs3113355','rs13679906',
'rs315502','rs356184',
'rs3497219','rs2737018',
'rs356186','rs8037141',
'rs356187','rs2748330',
'rs356191','rs10135063',
'rs356192',
'rs356195',
'rs356197',
'rs356198',
'rs3662777','rs356162',
'rs4931916','rs2736993',
'rs5235655','rs990086',
'rs6044643','rs1372518',
'rs609100','rs356188',
'rs6163618','rs33965306',
'rs6184008','rs2737023',
'rs6587182','rs2737028',
'rs67876224','rs11947452',
'rs7018400','rs2619344',
'rs7284835','rs2583964',
'rs7450824','rs986609',
'rs748849','rs4638027',
'rs7539711','rs356185',
'rs7686527','rs2737000',
'rs7774127','rs2619358',
'rs7863432','rs2737009',
'rs8722590','rs2737005',
'rs8732156','rs2583988',
'rs8993173','rs2737017',
'rs9164108','rs2583985',
'rs9225949','rs2583975',
'rs9630024','rs990087',
'rs986610','rs1408496',
'rs990085','rs4001644',
'rs9935229','rs2619350'
] 




DIR=os.getcwd()
TagSnpFile='/home/ajl1213/Projects/PD/SciAdv_Anal/GWAS_Enrich/1_CollectSNP/rawGWAS/MergedPD.rawGWAS.txt'
LdSnpFile='/home/ajl1213/Projects/PD/SciAdv_Anal/GWAS_Enrich/1_CollectSNP/LdFinal/MergedPD.allchr.txt'
#snpList=['rs4951254', 'rs14345696', 'rs2153904', 'rs13627427', 'rs200752153', 'rs79716096', 'rs1418549', 'rs451724', 'rs6689008', 'rs61824662', 'rs823106', 'rs10551380']

os.system('mkdir '+DIR+'/BedFiles')


## Tag SNPs
input1=open(TagSnpFile,'r')
all_input1=input1.readlines()
output1=open(DIR+'/BedFiles/MergedPD.TargetTagSnp.bed','w')
for line in all_input1[1:]:
    each=line.strip().split('\t')

    rsIDs=each[0]
    chrID=each[1]
    pt=int(each[2])

    for snpID in rsIDs.split(';'):
	if snpID in snpList:
	    output1.write(chrID+'\t'+str(pt)+'\t'+str(pt+1)+'\t'+snpID+'\n')

input1.close()
output1.close()

## LD SNPs
input1=open(LdSnpFile,'r')
all_input1=input1.readlines()
output1=open(DIR+'/BedFiles/MergedPD.TargetLdSnp.bed','w')
for line in all_input1[1:]:
    each=line.strip().split('\t')

    LdSnpPos=each[0]
    ParentalSnpPos=each[1]
    ParentalSnpID=each[2]

    for snpID in ParentalSnpID.split(';'):
	if snpID in snpList:

	    chrID=LdSnpPos.split('.')[0]
	    pt=int(LdSnpPos.split('.')[1])

	    output1.write(chrID+'\t'+str(pt)+'\t'+str(pt+1)+'\n')

input1.close()
output1.close()




