#!python3

import scrublet 
import scipy.io
import numpy
import os
import gzip


DIR=os.getcwd()
os.system('mkdir '+DIR+'/DoubletStat')
sampleList=[
'X5628NOSN',
'X4689NOSN',
'X4870NOSN',
'X4996NOSN',
'X5006NOSN',
'X5130NOSN',

'Public1NOSN',
'Public2_1NOSN',
'Public2_2NOSN',
'Public3_1NOSN',
'Public3_2NOSN',
'Public4NOSN',
'Public5NOSN',

'X5591PDSN',
'X4831PDSN',
'X5215PDSN',
'X5244PDSN',
'X5742PDSN',
'X5778PDSN'
]

metaDict={   # Gender / Age of death / PM interval
'X5628NOSN':['Female','80','12'],
'X4689NOSN':['Female','79','6'],
'X4870NOSN':['Female','63','8'],
'X4996NOSN':['Male','91','6'],
'X5006NOSN':['Male','69','24'],
'X5130NOSN':['Male','71','2'],

'Public1NOSN':['Male','59','NA'],
'Public2_1NOSN':['Male','70','NA'],
'Public2_2NOSN':['Male','70','NA'],
'Public3_1NOSN':['Female','56','NA'],
'Public3_2NOSN':['Female','56','NA'],
'Public4NOSN':['Male','70','NA'],
'Public5NOSN':['Male','55','NA'],

'X5591PDSN':['Male','76','8'],
'X4831PDSN':['Female','75','3'],
'X5215PDSN':['Female','89','NA'],
'X5244PDSN':['Female','82','5'],
'X5742PDSN':['Male','77','8'],
'X5778PDSN':['Male','71','8']
}


def attachStat():

	output1=open(DIR+'/DoubletStat/DoubletStatMerged.txt','w')
	output1.write('sampleID\tbarcodeID\tdoubletScore\tdoubletStat\tsex\tage\tPMI\n')
	for sampleID in sampleList:
		print(sampleID)

		inputDIR=DIR+'/'+sampleID+'/outs/filtered_feature_bc_matrix'

		os.system('zcat '+inputDIR+'/barcodes.tsv.gz > '+inputDIR+'/barcodes.tsv')
		os.system('zcat '+inputDIR+'/features.tsv.gz > '+inputDIR+'/features.tsv')
		os.system('zcat '+inputDIR+'/matrix.mtx.gz > '+inputDIR+'/matrix.mtx')

		## Load barcodes
		barcodeList=[]
		input1=open(inputDIR+'/barcodes.tsv', 'r')
		all_input1=input1.readlines()
		sampleIdx=str(sampleList.index(sampleID)+1)
		for line in all_input1:
			each=line.strip().split('\t')
			barcodeID=each[0].split('-')[0]+'-'+sampleIdx
			barcodeList.append(barcodeID)
		input1.close()

		## run scrublet
		counts_matrix = scipy.io.mmread(inputDIR +'/matrix.mtx').T.tocsc()
		genes = numpy.array(scrublet.load_genes(inputDIR +'/features.tsv', delimiter='\t', column=1))


		scrub = scrublet.Scrublet(counts_matrix, expected_doublet_rate=0.1)

		doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
									  min_cells=3, 
									  min_gene_variability_pctl=85, 
									  n_prin_comps=30)

		## make data
		output2=open(DIR+'/DoubletStat/'+sampleID+'.doublet.txt','w')
		output2.write('sampleID\tbarcodeID\tdoubletScore\tdoubletStat\tsex\tage\tPMI\n')
		for idx in range(0, len(doublet_scores)):
			newLine=[sampleID, barcodeList[idx], str(doublet_scores[idx]), str(predicted_doublets[idx]), metaDict[sampleID][0], metaDict[sampleID][1], metaDict[sampleID][2]]
			output1.write('\t'.join(newLine)+'\n')
			output2.write('\t'.join(newLine)+'\n')
		output2.close()
	output1.close()



attachStat()



