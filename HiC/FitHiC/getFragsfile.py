#!/usr/bin/python


import os
import sys
import gzip


DIR=os.getcwd()
inputDIR=sys.argv[1]
sampleID=sys.argv[2]
chrID=sys.argv[3]
resolution=sys.argv[4]


input1=open(inputDIR+'/'+sampleID+'.'+chrID+'.'+resolution+'.coverage','r')
all_input1=input1.readlines()
output1=gzip.open(DIR+'/Fragsfile/fragsfile.'+sampleID+'.'+chrID+'.'+resolution+'.gz','w')
for i in all_input1:
    each=i.split()

    coordID=each[0]+'.'+each[1]+'.'+each[2]
    mid=(int(each[1])+int(each[2]))/2
    coverage=each[3]

    output1.write(chrID+'\t'+coordID+'\t'+str(mid)+'\t'+coverage+'\t.\n')  

input1.close()
output1.close()









