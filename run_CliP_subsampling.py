'''----------------------------------------------------------------------
This script takes care of the running of CliP in case subsampling is needed
Usually you will need to run CliP on HPC; 
the preprocess script makes the input anonymous, so that you can up load them to HPC
Initialized by Kaixian Yu
Date: 05/04/2018
Email: kaixiany@163.com
----------------------------------------------------------------------
This script takes the following argument: path_to_input path_to_output path_to_clip lam No_subsampling Rep_id window_size overlap_size
-----------------------------------------------------------------------
Debug use
sys.argv = ['/Users/kaixiany/Working/CliP/Sample_data/intermediate/', '/Users/kaixiany/Working/CliP/Sample_data/results/', '/Users/kaixiany/Working/CliP/', '1.5', '1500', '1', '0.05', '0']
'''
import os
import sys
import numpy as np
sys.path.insert(0,sys.argv[3])
from CliP import *
from numpy import genfromtxt
prefix = sys.argv[1]
if not os.path.exists(sys.argv[2]):
    os.makedirs(sys.argv[2])
Lambda    = float(sys.argv[4])
No_subsampling = int(sys.argv[5])
rep = int(sys.argv[6])
window_size = float(sys.argv[7])
overlap = float(sys.argv[8])
r_all               = genfromtxt(prefix+"r.txt", delimiter="\t")
n_all               = genfromtxt(prefix+"n.txt", delimiter="\t")
minor_all           = genfromtxt(prefix+"minor.txt", delimiter="\t")
total_all           = genfromtxt(prefix+"total.txt", delimiter="\t")
purity              = genfromtxt(prefix+"purity_ploidy.txt", delimiter="\t")
coef_all            = genfromtxt(prefix+"coef.txt", delimiter="\t")
phicut_all          = genfromtxt(prefix+"cutbeta.txt", delimiter=" ")
No_mutation_all     = len(r_all)
ploidy              = 2
theta_hat_all 	    = r_all / n_all
phi_hat_all         = ploidy / (minor_all / theta_hat_all - total_all + ploidy)
scale_parameter_all = np.max( [1, np.max(phi_hat_all)] )
phi_all             = phi_hat_all / scale_parameter_all

ending = np.array([window_size])
current = window_size
flag = True
while flag:
	flag = False
	current += window_size - overlap
	ending = np.append(ending, np.min([1,current]))
	if current < 1:
		flag = True
starting = ending - window_size

index = []
SNVcount = np.zeros(len(ending))
index = [np.where( (phi_all > starting[i]) * (phi_all <= ending[i])  )[0] for i in range(len(ending))]
SNVcount = [len(index[i]) for i in range(len(index))]
sampling_proportion = No_subsampling / No_mutation_all
SNVtoSample = [np.round(SNVcount[i]*sampling_proportion) for i in range(len(SNVcount))]
np.random.seed(rep)
sample_index = [ np.random.choice(index[i],int(SNVtoSample[i]), False) for i in np.where(np.array(SNVtoSample,dtype='int') > 0)[0] ]
sample_index = np.concatenate(sample_index).ravel()

r = r_all[sample_index]
n = n_all[sample_index]
minor = minor_all[sample_index]
total = total_all[sample_index]
coef = coef_all[sample_index,:]
phicut = phicut_all[sample_index,:]
No_mutation = len(r)

alpha             = 0.8;
gamma             = 3.7;
rho               = 1.02;
precision         = 0.01;
Run_limit         = 1e4;
control_large     = 5;
post_th           = 0.05;
least_diff        = 0.01;
least_mut         = np.ceil(0.05 * No_mutation);
wcut=phicut
res = CliP(r, n, minor, total, ploidy, Lambda*0.01, alpha, rho, gamma, Run_limit, precision,
         control_large, least_mut, post_th, least_diff, coef, wcut)
labl = np.unique(res['label'])
summary = np.zeros([len(labl),3])

for i in range(len(labl)):
	summary[i,0] = labl[i]
	summary[i,2] = np.round(np.unique( res['phi'][np.where(res['label']==labl[i] )[0]])[0],3)
	summary[i,1] = len(np.where(res['label']==labl[i] )[0])	 
		 
np.savetxt('%s/lam%s_rep%s.txt'%(sys.argv[2], str(Lambda),str(rep)), summary ,fmt='%d\t%d\t%.3f')