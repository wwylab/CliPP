'''----------------------------------------------------------------------
This script takes care of the running of CliP
Usually you will need to run CliP on HPC; 
the preprocess script makes the input anonymous, so that you can up load them to HPC
Initialized by Kaixian Yu
Date: 05/04/2018
Email: kaixiany@163.com
----------------------------------------------------------------------
This script takes the following argument: path_to_input path_to_output path_to_clip lam
-----------------------------------------------------------------------
Debug use
sys.argv = ['/Users/kaixiany/Working/CliP/Sample_data/intermediate/', '/Users/kaixiany/Working/CliP/Sample_data/results/', '/Users/kaixiany/Working/CliP/', '1.5']
'''
import os
import sys
import numpy as np
sys.path.insert(0,sys.argv[3])
from CliP import *
from numpy import genfromtxt
prefix        = sys.argv[1]
if not os.path.exists(sys.argv[2]):
    os.makedirs(sys.argv[2])
Lambda        = float(sys.argv[4])
r             = genfromtxt(prefix+"r.txt", delimiter="\t")
n             = genfromtxt(prefix+"n.txt", delimiter="\t")
minor         = genfromtxt(prefix+"minor.txt", delimiter="\t")
total         = genfromtxt(prefix+"total.txt", delimiter="\t")
purity        = genfromtxt(prefix+"purity_ploidy.txt", delimiter="\t")
coef          = genfromtxt(prefix+"coef.txt", delimiter="\t")
phicut        = genfromtxt(prefix+"cutbeta.txt", delimiter=" ")
No_mutation   = len(r);
ploidy        = 2
alpha         = 0.8;
gamma         = 3.7;
rho           = 1.02;
precision     = 0.01;
Run_limit     = 1e4;
control_large = 5;
post_th       = 0.05;
least_diff    = 0.01;
least_mut     = np.ceil(0.05 * No_mutation);

wcut=phicut

res = CliP(r, n, minor, total, ploidy, Lambda, alpha, rho, gamma, Run_limit, precision, control_large, least_mut, post_th, least_diff, coef, wcut, purity)

labl = np.unique(res['label'])
summary = np.zeros([len(labl),3])

for i in range(len(labl)):
	summary[i,0] = labl[i]
	summary[i,2] = np.round(np.unique( res['phi'][np.where(res['label']==labl[i] )[0]])[0],3)
	summary[i,1] = len(np.where(res['label']==labl[i] )[0])	

np.savetxt('%s/lam%s_phi.txt'%(sys.argv[2],str(Lambda)), res['phi'],fmt='%.3f', delimiter = ',')
np.savetxt('%s/lam%s_label.txt'%(sys.argv[2], str(Lambda)), res['label'],fmt='%d', delimiter = ',')
np.savetxt('%s/lam%s_summary_table.txt'%(sys.argv[2], str(Lambda)), summary ,fmt='%d\t%d\t%.3f')
