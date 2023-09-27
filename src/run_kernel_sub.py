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
import os
import sys
import numpy as np
from numpy import genfromtxt
import ctypes
import glob

#prefix = sys.argv[1]

#if not os.path.exists(sys.argv[2]):
#    os.makedirs(sys.argv[2])

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


def run_clip_sub(prefix, preliminary_result, Lambda_list, No_subsampling, rep, window_size, overlap):
    if not os.path.exists(preliminary_result):
        os.makedirs(preliminary_result)

    current_folder = os.path.dirname(os.path.abspath(__file__))
    clip_lib_path = glob.glob(os.path.join(current_folder, "../build/*/CliP*%s*.so" %(sys.platform)))
    if not clip_lib_path:
        sys.stderr.write("Cannot find shared lib. Make sure run run python setup.py build first.\n")
        sys.exit(0)

    clip_lib_path = clip_lib_path[0]
    if not os.path.isfile(clip_lib_path):
        sys.stderr.write("Cannot find shared lib. Make sure run run python setup.py build first.\n")
        sys.exit(0)
        
    clip_lib = ctypes.CDLL(clip_lib_path)

    Lambda_num = len(Lambda_list)
    Lambda_list = np.array(Lambda_list).astype(np.float64)
    
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

    alpha = 0.8
    gamma = 3.7
    rho = 1.02
    precision = 0.01
    Run_limit = 1e4
    control_large = 5
    post_th = 0.05
    least_diff = 0.01

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

        least_mut = np.ceil(0.05 * No_mutation)
        wcut = phicut
        
        r = r.astype(np.int32)
        n = n.astype(np.int32)
        minor = minor.astype(np.int32)
        total = total.astype(np.int32)
        coef = coef.flatten()
        wcut = wcut.flatten()
        
        #ploidy = double(ploidy)
        Run_limit = int(Run_limit)
        least_mut = int(least_mut)
        
        coef = coef.astype(np.float64)
        wcut = wcut.astype(np.float64)
        
        
        clip_lib.CliP.restype = None
        clip_lib.CliP.argtypes = [ctypes.c_int, 
                                  np.ctypeslib.ndpointer(dtype=np.int32),
                                  np.ctypeslib.ndpointer(dtype=np.int32),
                                  np.ctypeslib.ndpointer(dtype=np.int32), 
                                  np.ctypeslib.ndpointer(dtype=np.int32), 
                                  ctypes.c_double,
                                  np.ctypeslib.ndpointer(dtype=np.float64), 
                                  ctypes.c_int, 
                                  ctypes.c_double, 
                                  ctypes.c_double, 
                                  ctypes.c_double, 
                                  ctypes.c_int, 
                                  ctypes.c_double, 
                                  ctypes.c_int,
                                  ctypes.c_int, 
                                  ctypes.c_double, 
                                  ctypes.c_double, 
                                  np.ctypeslib.ndpointer(dtype=np.float64), 
                                  np.ctypeslib.ndpointer(dtype=np.float64), 
                                  ctypes.c_double, 
                                  ctypes.c_char_p]
        
        clip_lib.CliP(No_mutation, r, n, minor, total, ploidy,
                Lambda_list, Lambda_num, alpha, rho, gamma, Run_limit, precision,
                control_large, least_mut, post_th, least_diff,
                coef, wcut, purity, preliminary_result.encode('utf-8'))

        for _lambda in Lambda_list:
            _phi_file_path = os.path.join(preliminary_result, "lam%s_phi.txt" %(_lambda))
            _phi_file_path_new = os.path.join(preliminary_result, "lam%s_phi_rep%s.txt" %(_lambda, j))

            _label_file_path = os.path.join(preliminary_result, "lam%s_label.txt" %(_lambda))
            _label_file_path_new = os.path.join(preliminary_result, "lam%s_label_rep%s.txt" %(_lambda, j))

            _sum_old_file_path = os.path.join(preliminary_result, "lam%s_summary_table.txt" %(_lambda))
            _sum_new_file_path = os.path.join(preliminary_result, "lam%s_rep%s.txt" %(_lambda, j))

            try:
                # os.remove(_phi_file_path)
                # os.remove(_label_file_path)
                os.rename(_sum_old_file_path, _sum_new_file_path)
                os.rename(_phi_file_path, _phi_file_path_new)
                os.rename(_label_file_path, _label_file_path_new)
            except Exception as err:
                sys.stderr.write(err)
