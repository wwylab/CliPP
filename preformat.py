'''
Debug use:
sys.argv = ['','Sample_data/Tumour2.mutect.vcf', 'Sample_data/Tumour2.battenberg.txt']

'''
import numpy as np
import sys
import os
os.environ["OMP_NUM_THREADS"] = "30"
#### First process the SNV file in VCF format ####
fin = open(sys.argv[1],'r')
lines = fin.readlines()
fin.close()
num_lines = len(lines)

fout = open('/intermediate/snv_input.txt','w')
for i in range(num_lines):
	if lines[i].startswith('##'):
		continue
	elif lines[i].startswith('#C'):
		tmp = np.array(lines[i].strip().split('\t'))
		tumorIndex = np.where(tmp=='tumor')[0][0]
	else:
		tmp = lines[i].strip().split('\t')
		output = tmp[:2]
		tumorInfo = tmp[tumorIndex]
		tmpReads = tumorInfo.strip().split(':')[1].split(',')
		output.append(tmpReads[1])
		output.append(tmpReads[0])
		_ = fout.write('\t'.join(output) + '\n')
fout.close()

#### Taking care of the CNV file ####
fin = open(sys.argv[2],'r')
lines = fin.readlines()
fin.close()
num_lines = len(lines)

fout = open('/intermediate/cnv_input.txt','w')
for i in range(1,num_lines):
	tmp = lines[i].strip().split('\t')
	output = tmp[:3]
	output += tmp[7:9]
	_ = fout.write('\t'.join(output) + '\n')
fout.close()

#### Taking care of the purity file ####
fin = open(sys.argv[3], 'r')
lines = fin.readlines()
fin.close()
fout = open('/intermediate/purity.txt','w')
_ = fout.write(lines[1].strip().split('\t')[0])
fout.close()





