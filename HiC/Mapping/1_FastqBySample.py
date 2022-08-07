#!/python2



import os


DIR=os.getcwd()
inputDIR=DIR+'/RawFastq'
ncpu=4
sampleList=[
'X4870NOSN;JW28-JW57-JW230',
'X4689NOSN;JW29-JW58-JW231',
'X5628NOSN;JW79-JW176-JW232',
'X5130NOSN;JW80-JW177-JW233-JW254',
'X5006NOSN;JW81-JW448',
'X4996NOSN;JW82-JW450',
'X5244PDSN;JW40-JW59',
'X4831PDSN;JW41-JW60-JW250',
'X5215PDSN;JW161-JW86-JW191-JW251',
'X5649PDSN;JW157-JW162-JW252',
'X5591PDSN;JW158-JW163-JW192-JW253'
]

os.system('mkdir '+DIR+'/SampleFastq')



def makeScript():

    for i in sampleList:
	sampleID=i.split(';')[0]
	rawFileList=i.split(';')[1].split('-')

	for orient in ['_1','_2']:
	    output1=open('PBS_hic_fastq_'+sampleID+orient+'.pbs','w')
	    output1.write('#PBS -N hic_fastq_'+sampleID+orient+'\n')
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

	    line_elements=['zcat']
	    for rawFile in rawFileList:
		full_id=inputDIR+'/'+rawFile+orient+'.fq.gz'
		line_elements.append(full_id)
	    line_elements.append('| gzip -c >')
	    line_elements.append(DIR+'/SampleFastq/'+sampleID+orient+'.fq.gz')
	    runLine=' '.join(line_elements)

	    output1.write(runLine+'\n')

	    output1.write('\n')
	    output1.write('sleep 30\n')
	    output1.write('exit 0')
	    output1.close()

	    os.system('qsub PBS_hic_fastq_'+sampleID+orient+'.pbs')


makeScript()


