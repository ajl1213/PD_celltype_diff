#!python3

import scrublet 
import scipy.io
import numpy
import os


DIR=os.getcwd()
os.system('mkdir '+DIR+'/DoubletStat')
sampleList=[
'X4870NOSN',
'X4996NOSN',
'X5006NOSN',
'X5130NOSN',
'X5628NOSN',

'Public1NOSN',
'Public2NOSN',

'X4831PDSN',
'X5215PDSN',
'X5244PDSN',
'X5532PDSN',
'X5591PDSN',
'X5627PDSN',
'X5742PDSN',
'X5778PDSN'
]

metaDict={
'X4870NOSN':['Female','63','8'],
'X4996NOSN':['Male','91','6'],
'X5006NOSN':['Male','69','24'],
'X5130NOSN':['Male','71','2'],
'X5628NOSN':['Female','80','12'],

'Public1NOSN':['Male','86','NA'], 
'Public2NOSN':['Female','83','NA'],

'X4831PDSN':['Female','75','3'],
'X5215PDSN':['Female','89','NA'],
'X5244PDSN':['Female','82','5'],
'X5532PDSN':['Male','75','6'],
'X5591PDSN':['Male','76','8'],
'X5627PDSN':['Male','68','7'],
'X5742PDSN':['Male','77','8'],
'X5778PDSN':['Male','71','8']
}


def attachStat():

	outFinal=open(DIR+'/DoubletStat/DoubletStatMerged.txt','w')
	outFinal.write('barcodeID\torig.ident\tdoubletScore\tdoubletStat\tsex\tage\tPMI\n')
	for sampleID in sampleList:
		print(sampleID)

		inputDIR=DIR+'/'+sampleID+'/outs/filtered_peak_bc_matrix'

		## Make feature input
		input1=open(inputDIR+'/peaks.bed','r')
		output1=open(inputDIR+'/peaks.txt','w')
		all_input1=input1.readlines()
		for line in all_input1:
			each=line.strip().split('\t')
			output1.write('-'.join(each)+'\n')
		input1.close()
		output1.close()

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
		genes = numpy.array(scrublet.load_genes(inputDIR +'/peaks.txt', delimiter='\t', column=0))

		scrub = scrublet.Scrublet(counts_matrix, expected_doublet_rate=0.1)

		doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
									  min_cells=3, 
									  min_gene_variability_pctl=85, 
									  n_prin_comps=30)
		## make data
		output2=open(DIR+'/DoubletStat/'+sampleID+'.doublet.txt','w')
		output2.write('barcodeID\torig.ident\tdoubletScore\tdoubletStat\tsex\tage\tPMI\n')
		for idx in range(0, len(doublet_scores)):
			newLine=[barcodeList[idx], sampleID, str(doublet_scores[idx]), str(predicted_doublets[idx]), metaDict[sampleID][0], metaDict[sampleID][1], metaDict[sampleID][2]]
			outFinal.write('\t'.join(newLine)+'\n')
			output2.write('\t'.join(newLine)+'\n')
		output2.close()
	outFinal.close()


attachStat()



