'''----------------------------------------------------------------------
This script takes care of the running of CliP in case subsampling is needed
Usually you will need to run CliP on HPC; 
the preprocess script makes the input anonymous, so that you can up load them to HPC
Authors: Kaixian Yu, Yujie Jiang
Date: 04/02/2021
Email: yujiejiang679@gmail.com
----------------------------------------------------------------------
This script takes the following argument: path_to_input path_to_output path_to_clip lam No_subsampling Rep_num window_size overlap_size
-----------------------------------------------------------------------
Debug use
sys.argv = ['/Users/kaixiany/Working/CliP/Sample_data/intermediate/', '/Users/kaixiany/Working/CliP/Sample_data/results/', '/Users/kaixiany/Working/CliP/', '1.5', '1500', '1', '0.05', '0']
'''
import os
import sys
import time
import numpy as np
import ctypes
## For parallel
# mkl_rt = ctypes.CDLL('libmkl_rt.so')
# mkl_get_max_threads = mkl_rt.mkl_get_max_threads
# mkl_rt.mkl_set_num_threads(ctypes.byref(ctypes.c_int(16)))

sys.path.insert(0,sys.argv[3])
from kernel import *
from numpy import genfromtxt
prefix = sys.argv[1]
if not os.path.exists(sys.argv[2]):
    os.makedirs(sys.argv[2])
No_subsampling = int(sys.argv[4])
rep = int(sys.argv[5])
window_size = float(sys.argv[6])
overlap = float(sys.argv[7])
#sample_id = sys.argv[9]

r_all               = genfromtxt(os.path.join(prefix, "r.txt"), delimiter="\t")
n_all               = genfromtxt(os.path.join(prefix, "n.txt"), delimiter="\t")
minor_all           = genfromtxt(os.path.join(prefix, "minor.txt"), delimiter="\t")
total_all           = genfromtxt(os.path.join(prefix, "total.txt"), delimiter="\t")
purity              = genfromtxt(os.path.join(prefix, "purity_ploidy.txt"), delimiter="\t")
coef_all            = genfromtxt(os.path.join(prefix, "coef.txt"), delimiter="\t")
phicut_all          = genfromtxt(os.path.join(prefix, "cutbeta.txt"), delimiter=" ")
No_mutation_all     = len(r_all)
ploidy              = 2
theta_hat_all       = r_all / n_all
phi_hat_all         = theta_hat_all * ((ploidy - purity*ploidy + purity*total_all) / minor_all)
scale_parameter_all = np.max( [1, np.max(phi_hat_all)] )
phi_all             = phi_hat_all / scale_parameter_all
#############################
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


# start = time.time()
for j in range(1,rep+1):
    np.random.seed(j)
    sample_index = [ np.random.choice(index[i],int(SNVtoSample[i]), False) for i in np.where(np.array(SNVtoSample,dtype='int') > 0)[0] ]
    sample_index = np.concatenate(sample_index).ravel()
    ##############################
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
    Lambda_list = [0.01, 0.03, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25]
    wcut=phicut

    if len(sys.argv) < 9:
        for Lambda in Lambda_list:
            res = CliP(r, n, minor, total, ploidy, Lambda*0.01, alpha, rho, gamma, Run_limit, precision, 
            control_large, least_mut, post_th, least_diff, coef, wcut, purity)
            labl = np.unique(res['label'])
            summary = np.zeros([len(labl),3])
		
            for i in range(len(labl)):
                summary[i,0] = labl[i]
                summary[i,2] = np.round(np.unique( res['phi'][np.where(res['label']==labl[i] )[0]])[0],3)
                summary[i,1] = len(np.where(res['label']==labl[i] )[0])
		
            np.savetxt('%s/lam%s_rep%s.txt'%(sys.argv[2], str(Lambda),str(j)), summary ,fmt='%d\t%d\t%.3f')
    
    else:
        Lambda    = float(sys.argv[8])
        res = CliP(r, n, minor, total, ploidy, Lambda*0.01, alpha, rho, gamma, Run_limit, precision, 
        control_large, least_mut, post_th, least_diff, coef, wcut, purity)
        labl = np.unique(res['label'])
        summary = np.zeros([len(labl),3])
	
        for i in range(len(labl)):
            summary[i,0] = labl[i]
            summary[i,2] = np.round(np.unique( res['phi'][np.where(res['label']==labl[i] )[0]])[0],3)
            summary[i,1] = len(np.where(res['label']==labl[i] )[0])
	
        np.savetxt('%s/lam%s_rep%s.txt'%(sys.argv[2], str(Lambda),str(j)), summary ,fmt='%d\t%d\t%.3f')
    
# end = time.time()
# print(" Time elapsed: ", end - start, "seconds")
