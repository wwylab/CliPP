'''----------------------------------------------------------------------
This script takes care of the running of CliP in case subsampling is needed
Usually you will need to run CliP on HPC; 
the preprocess script makes the input anonymous, so that you can up load them to HPC
Authors: Kaixian Yu, Yujie Jiang
Date: 04/02/2021
Email: yujiejiang679@gmail.com
----------------------------------------------------------------------
This script takes the following argument: path_to_input path_to_output path_to_clip No_subsampling Rep_num window_size overlap_size lam
-----------------------------------------------------------------------
Debug use
sys.argv = ['/Users/kaixiany/Working/CliP/Sample_data/intermediate/', '/Users/kaixiany/Working/CliP/Sample_data/results/', '/Users/kaixiany/Working/CliP/', '1.5', '1500', '1', '0.05', '0']
'''
import concurrent
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, ALL_COMPLETED
import numpy as np

sys.path.insert(0,"./src/")

from src.kernel import CliP

sys.path.insert(0,sys.argv[3])
from numpy import genfromtxt
prefix = sys.argv[1]

if not os.path.exists(sys.argv[2]):
    os.makedirs(sys.argv[2])

'''----------------------------------------------------------------------
This script takes care of the running of CliP
Usually you will need to run CliP on HPC;
the preprocess script makes the input anonymous, so that you can up load them to HPC
Authors: Kaixian Yu, Yujie Jiang, Yuxin Tang
Date: 04/02/2021
Email: yujiejiang679@gmail.com
----------------------------------------------------------------------
This script takes the following argument: path_to_input path_to_output path_to_clip lam
-----------------------------------------------------------------------
Debug use
sys.argv = ['/Users/kaixiany/Working/CliP/Sample_data/intermediate/', '/Users/kaixiany/Working/CliP/Sample_data/results/', '/Users/kaixiany/Working/CliP/', '1.5']
'''


def clip_kernel_sub(preliminary_result,j, r, n, minor, total, ploidy, Lambda, alpha, rho, gamma, Run_limit, precision, control_large, least_mut, post_th, least_diff, coef, wcut, purity):
    res = CliP(r, n, minor, total, ploidy, Lambda, alpha, rho, gamma, Run_limit, precision, control_large, least_mut, post_th, least_diff, coef, wcut, purity)
    labl = np.unique(res['label'])
    summary = np.zeros([len(labl), 3])

    for i in range(len(labl)):
        summary[i, 0] = labl[i]
        summary[i, 2] = np.round(np.unique(res['phi'][np.where(res['label'] == labl[i])[0]])[0], 3)
        summary[i, 1] = len(np.where(res['label'] == labl[i])[0])

    np.savetxt('%s/lam%s_rep%s.txt' % (preliminary_result, str(Lambda), str(j)), summary, fmt='%d\t%d\t%.3f')
    return 1


def run_clip_sub(prefix, preliminary_result, _lambda_list, No_subsampling, rep, window_size, overlap):
    if not os.path.exists(preliminary_result):
        os.makedirs(preliminary_result)

    r_all = genfromtxt(os.path.join(prefix, "r.txt"), delimiter="\t")
    n_all = genfromtxt(os.path.join(prefix, "n.txt"), delimiter="\t")
    minor_all = genfromtxt(os.path.join(prefix, "minor.txt"), delimiter="\t")
    total_all = genfromtxt(os.path.join(prefix, "total.txt"), delimiter="\t")
    purity = genfromtxt(os.path.join(prefix, "purity_ploidy.txt"), delimiter="\t")
    coef_all = genfromtxt(os.path.join(prefix, "coef.txt"), delimiter="\t")
    phicut_all = genfromtxt(os.path.join(prefix, "cutbeta.txt"), delimiter=" ")
    No_mutation_all = len(r_all)
    ploidy = 2
    theta_hat_all = r_all / n_all
    phi_hat_all = theta_hat_all * ((ploidy - purity * ploidy + purity * total_all) / minor_all)
    scale_parameter_all = np.max([1, np.max(phi_hat_all)])
    phi_all = phi_hat_all / scale_parameter_all
    #############################
    ending = np.array([window_size])
    current = window_size
    flag = True
    while flag:
        flag = False
        current += window_size - overlap
        ending = np.append(ending, np.min([1, current]))
        if current < 1:
            flag = True
    starting = ending - window_size

    index = []
    SNVcount = np.zeros(len(ending))
    index = [np.where((phi_all > starting[i]) * (phi_all <= ending[i]))[0] for i in range(len(ending))]
    SNVcount = [len(index[i]) for i in range(len(index))]
    sampling_proportion = No_subsampling / No_mutation_all
    SNVtoSample = [np.round(SNVcount[i] * sampling_proportion) for i in range(len(SNVcount))]

    for j in range(1, rep + 1):
        np.random.seed(j)
        sample_index = [np.random.choice(index[i], int(SNVtoSample[i]), False) for i in
                        np.where(np.array(SNVtoSample, dtype='int') > 0)[0]]
        sample_index = np.concatenate(sample_index).ravel()
        ##############################
        r = r_all[sample_index]
        n = n_all[sample_index]
        minor = minor_all[sample_index]
        total = total_all[sample_index]
        coef = coef_all[sample_index, :]
        phicut = phicut_all[sample_index, :]
        No_mutation = len(r)

        alpha = 0.8
        gamma = 3.7
        rho = 1.02
        precision = 0.01
        Run_limit = 1e4
        control_large = 5
        post_th = 0.05
        least_diff = 0.01
        least_mut = np.ceil(0.05 * No_mutation)
        Lambda_list = [0.01, 0.03, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25]
        wcut = phicut

        # It is better to set # of cores = max_workers
        pool = ThreadPoolExecutor(max_workers=11)
        future_to_lambdas = []
        for Lambda in Lambda_list:
            future_to_lambdas.append(pool.submit(clip_kernel_sub, preliminary_result,j, r, n, minor, total, ploidy, Lambda, alpha, rho, gamma, Run_limit, precision,
                                             control_large, least_mut, post_th, least_diff, coef, wcut, purity))

        concurrent.futures.wait(future_to_lambdas,timeout=None,return_when=ALL_COMPLETED)
