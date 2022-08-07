#!/home/ajl1213/anaconda2/bin/python


import os


ncpu=4
DIR=os.getcwd()
snRNA_Data='/home/ajl1213/Projects/PD/SciAdv_Data/snRNAseq/02_SeuratProcess/SeuratObjects/PD.SN.snRNA.postAlign.Anno.rds'
cellTypeList=['DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri']

os.system('mkdir '+DIR+'/ZmatByCells')


def makeScript():
    for cellType in cellTypeList:
	output1=open('PBS_snRNA_zval.'+cellType+'.pbs','w')
	output1.write('#PBS -N snRNA_zval.'+cellType+'\n')
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

	output1.write('source activate scEnv\n')

	R_input=' \'--args '+snRNA_Data+' '+cellType+' '+DIR+'/ZmatByCells/'+cellType+'.snRNA.zscore.txt\''
	output1.write('R CMD BATCH --no-save --no-restore '+R_input+' 5b_getZmat.byCell.R '+DIR+'/snRNA_zval.'+cellType+'.Rout\n')

	output1.write('\n')
	output1.write('sleep 30\n')
	output1.write('exit 0')
	output1.close()

	os.system('qsub PBS_snRNA_zval.'+cellType+'.pbs')


makeScript()


