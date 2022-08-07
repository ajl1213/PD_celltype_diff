#!/usr/bin/python

import sys
import fileinput


#file_name=str(sys.argv[1])
output_name=str(sys.argv[1])

def identify_usable_cis():
    trans=0
    self=0
    cis=0
    check_id={}

    for line in sys.stdin:
		if line[0]!='@':
	            [id1,flag1,chr_from1,loc_from1,mapq1,cigar1,chr_from2, loc_from2, dist, read1, read_qual1]=line.split('\t')[0:11]
		    pos1_chr=chr_from1
		    pos2_chr=chr_from2
		    pos1=int(loc_from1)
		    pos2=int(loc_from2)
		    dist=abs(int(dist))
		    if check_id.has_key(id1):
			    del check_id[id1]
		    else:
			    check_id[id1]=1
			    if chr_from2=='=':
				    if dist>15000:
					    cis+=1
				    else:
					    self+=1
			    else:
				    trans+=1


    total=float(trans+cis+self)
    output1=open(output_name,'w')
    output1.write('trans\tself\tcis\n')
    output1.write(str(trans)+'\t'+str(self)+'\t'+str(cis)+'\n')
    output1.write(str(trans/total)+'\t'+str(self/total)+'\t'+str(cis/total)+'\n')
    output1.close()


identify_usable_cis()


