
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

import concurrent
import os
import sys
from concurrent.futures import ThreadPoolExecutor, ALL_COMPLETED
import numpy as np
from numpy import genfromtxt

sys.path.insert(0,"./src/")

from src.kernel import CliP

def clip_kernel_nosub(preliminary_result, r, n, minor, total, ploidy, Lambda, alpha, rho, gamma, Run_limit, precision, control_large, least_mut, post_th, least_diff, coef, wcut, purity):
    res = CliP(r, n, minor, total, ploidy, Lambda, alpha, rho, gamma, Run_limit, precision, control_large, least_mut, post_th, least_diff, coef, wcut, purity)
    labl = np.unique(res['label'])
    summary = np.zeros([len(labl), 3])

    for i in range(len(labl)):
        summary[i, 0] = labl[i]
        summary[i, 2] = np.round(np.unique(res['phi'][np.where(res['label'] == labl[i])[0]])[0], 3)
        summary[i, 1] = len(np.where(res['label'] == labl[i])[0])

    np.savetxt('%s/lam%s_phi.txt' % (preliminary_result, str(Lambda)), res['phi'], fmt='%.3f', delimiter=',')
    np.savetxt('%s/lam%s_label.txt' % (preliminary_result, str(Lambda)), res['label'], fmt='%d', delimiter=',')
    np.savetxt('%s/lam%s_summary_table.txt' % (preliminary_result, str(Lambda)), summary, fmt='%d\t%d\t%.3f')
    return 1

def run_clip_nosub(prefix, preliminary_result):
    if not os.path.exists(preliminary_result):
        os.makedirs(preliminary_result)

    r = genfromtxt(prefix+"r.txt", delimiter="\t")
    n = genfromtxt(prefix+"n.txt", delimiter="\t")
    minor = genfromtxt(prefix+"minor.txt", delimiter="\t")
    total = genfromtxt(prefix+"total.txt", delimiter="\t")
    purity = genfromtxt(prefix+"purity_ploidy.txt", delimiter="\t")
    coef = genfromtxt(prefix+"coef.txt", delimiter="\t")
    phicut = genfromtxt(prefix+"cutbeta.txt", delimiter=" ")
    wcut = phicut
    No_mutation = len(r)
    least_mut = np.ceil(0.05 * No_mutation)

    ploidy = 2
    alpha = 0.8
    gamma = 3.7
    rho = 1.02
    precision = 0.01
    Run_limit = 1e4
    control_large = 5
    post_th = 0.05
    least_diff = 0.01
    Lambda_list = [0.01, 0.03, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25]

    # It is better to set # of cores = max_workers
    pool = ThreadPoolExecutor(max_workers=11)
    future_to_lambdas = []
    for Lambda in Lambda_list:
        future_to_lambdas.append(pool.submit(clip_kernel_nosub, preliminary_result, r, n, minor, total, ploidy, Lambda, alpha, rho, gamma, Run_limit, precision,
                                             control_large, least_mut, post_th, least_diff, coef, wcut, purity))

    concurrent.futures.wait(future_to_lambdas,timeout=None,return_when=ALL_COMPLETED)
