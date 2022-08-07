#!/usr/bin/python

import fileinput
import sys



maxDist=int(sys.argv[1])
resolution=sys.argv[2]

interDict={}

if resolution=='5kb':
    binSize=5000
if resolution=='10kb':
    binSize=10000
if resolution=='20kb':
    binSize=20000
if resolution=='40kb':
    binSize=40000


for line in sys.stdin:
        # HWI-ST1113:549:HKWJFADXX:2:1115:20976:94471   163     chr22   16052163        41      64M36S  chrY    16541071        0       TCAAAAAACAAACAAACAAACAAACAAACAAACAAACTGTCAAAATCTGTACAGTATGTGAAGAAGCTTAGTTTTTTCGCTCTTTGCAATAAATCTTGCT    CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJIGIJDGJJJJIJJJJJJIIIIJGHIIJJHHGHHFFCEFFDDDDDDDDDDDDEDDDDDEDDDDDC    NM:i:0  MD:Z:64 AS:i:64 XS:i:51 SA:Z:chrY,16540953,+,62S38M,27,0;
        [id1,flag1,chr_from1,loc_from1,mapq1,cigar1,chr_from2, loc_from2, dist, read1, read_qual1]=line.split('\t')[0:11]

	if abs(int(dist))<=maxDist:

	       	frag1=int(loc_from1)/binSize*binSize
	       	frag2=int(loc_from2)/binSize*binSize
		 
		if frag1<frag2:
			id=chr_from1+'.'+str(frag1)+'.'+str(frag2)
		else:
			id=chr_from1+'.'+str(frag2)+'.'+str(frag1)
		if interDict.has_key(id):
			interDict[id]+=1
		else:
			interDict[id]=1

for i in interDict:
	info=i.split('.')
	print info[0]+'\t'+str(int(info[1])+(binSize/2))+'\t'+info[0]+'\t'+str(int(info[2])+(binSize/2))+'\t'+str(interDict[i]/2)
	
		
			

