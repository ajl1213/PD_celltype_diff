#!/home/ajl1213/anaconda2/bin/python


import os


DIR=os.getcwd()
PythonDIR='/home/ajl1213/anaconda2/bin'
ScriptDIR='/home/ajl1213/programs/ldsc'
PeakDIR='/home/ajl1213/Projects/PD/NeuronJ/DiffAnal/02_H3K27ac/1_BulkH3K27ac/BedFiles'
chrList=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']
dataTypeList=['Down','Up']
cellTypeList=['Merged', 'DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']

os.system('mkdir '+DIR+'/AnnoLdFile')


def getGenotypeDataFromPlink():
    runScript='wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_plinkfiles.tgz'
    print runScript
    os.system(runScript)

    runScript='tar -zxvf 1000G_Phase3_plinkfiles.tgz'
    print runScript
    os.system(runScript)


def getHapMap3Snps():
    runScript='wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/hapmap3_snps.tgz'
    print runScript
    os.system(runScript)

    runScript='tar -zxvf hapmap3_snps.tgz'
    print runScript
    os.system(runScript)


def makeAnnoFile():
    ncpu=2

    for dataType in dataTypeList:
	for cellType in cellTypeList:
	    for chrID in chrList:
		output1=open('PBS_makeAnno.'+dataType+'.'+cellType+'.'+chrID+'.pbs','w')
		output1.write('#PBS -N makeAnno.'+dataType+'.'+cellType+'.'+chrID+'\n')
		output1.write('#PBS -q workq\n')
		output1.write('#PBS -l nodes=1:ppn='+str(ncpu)+'\n')
		output1.write('#PBS -j oe\n')
		output1.write('\n')
		output1.write('# go workdir\n')
		output1.write('cd $PBS_O_WORKDIR\n')
		output1.write('\n')
		output1.write('# run command \n')
		output1.write('sleep 5\n')
		output1.write('\n')
		output1.write('echo -n \"I am on: \"\n')
		output1.write('hostname;\n')    
		output1.write('echo finding ssh-accessible nodes:\n')
		output1.write('echo -n \"running on: \"\n')
		output1.write('\n')

		output1.write(PythonDIR+'/python '+ScriptDIR+'/make_annot.py --bed-file '+PeakDIR+'/'+dataType+'.'+cellType+'.bed --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.'+chrID+'.bim --annot-file '+DIR+'/AnnoLdFile/'+dataType+'.'+cellType+'.'+chrID+'.anno.gz\n')
		output1.write('sleep 30\n')
		output1.write('exit 0')
		output1.close()

		os.system('qsub PBS_makeAnno.'+dataType+'.'+cellType+'.'+chrID+'.pbs')


def ComputeLD():
    ncpu=4
    for dataType in dataTypeList:
	for cellType in cellTypeList:
	    for chrID in chrList:
		output1=open('PBS_LdCompute.'+dataType+'.'+cellType+'.'+chrID+'.pbs','w')
		output1.write('#PBS -N LdCompute.'+dataType+'.'+cellType+'.'+chrID+'\n')
		output1.write('#PBS -q workq\n')
		output1.write('#PBS -l nodes=1:ppn='+str(ncpu)+'\n')
		output1.write('#PBS -j oe\n')
		output1.write('\n')
		output1.write('# go workdir\n')
		output1.write('cd $PBS_O_WORKDIR\n')
		output1.write('\n')
		output1.write('# run command \n')
		output1.write('sleep 5\n')
		output1.write('\n')
		output1.write('echo -n \"I am on: \"\n')
		output1.write('hostname;\n')    
		output1.write('echo finding ssh-accessible nodes:\n')
		output1.write('echo -n \"running on: \"\n')
		output1.write('\n')

		output1.write(PythonDIR+'/python '+ScriptDIR+'/ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.'+chrID+' --ld-wind-cm 1 --annot '+DIR+'/AnnoLdFile/'+dataType+'.'+cellType+'.'+chrID+'.anno.gz --thin-annot --out '+DIR+'/AnnoLdFile/'+dataType+'.'+cellType+'.'+chrID+' --print-snps hapmap3_snps/hm.'+chrID+'.snp\n')

		output1.write('\n')
		output1.write('sleep 30\n')
		output1.write('exit 0')
		output1.close()

		os.system('qsub PBS_LdCompute.'+dataType+'.'+cellType+'.'+chrID+'.pbs')



#getGenotypeDataFromPlink()
#getHapMap3Snps()

#makeAnnoFile()
ComputeLD()






