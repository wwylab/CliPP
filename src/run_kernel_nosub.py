
'''----------------------------------------------------------------------
This script takes care of the running of CliP
Usually you will need to run CliP on HPC;
the preprocess script makes the input anonymous, so that you can up load them to HPC
Authors: Kaixian Yu, Yujie Jiang, Shuangxi Ji. Yuxin Tang
Date: 04/02/2021
Email: yujiejiang679@gmail.com
----------------------------------------------------------------------
This script takes the following argument: path_to_input path_to_output path_to_clip lam
-----------------------------------------------------------------------
Debug use
sys.argv = ['/Users/kaixiany/Working/CliP/Sample_data/intermediate/', '/Users/kaixiany/Working/CliP/Sample_data/results/', '/Users/kaixiany/Working/CliP/', '1.5']
'''

import os
import sys
import numpy as np
from numpy import genfromtxt
import ctypes
import glob


def run_clip_nosub(prefix, preliminary_result, lambda_list):
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
    Lambda_list = lambda_list
    

    r = r.astype(np.int32)
    n = n.astype(np.int32)
    minor = minor.astype(np.int32)
    total = total.astype(np.int32)
    coef = coef.flatten()
    wcut = wcut.flatten()
        
    Lambda_num = len(Lambda_list)
    Lambda_list = np.array(Lambda_list).astype(np.float64)
        
    #ploidy = double(ploidy)
    Run_limit = int(Run_limit)
    least_mut = int(least_mut)
        
    coef = coef.astype(np.float64)
    wcut = wcut.astype(np.float64)

    clip_lib.CliP.restype = None

    clip_lib.CliP.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(dtype=np.int32),np.ctypeslib.ndpointer(dtype=np.int32),np.ctypeslib.ndpointer(dtype=np.int32), np.ctypeslib.ndpointer(dtype=np.int32), ctypes.c_double,np.ctypeslib.ndpointer(dtype=np.float64), ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_double, ctypes.c_int,ctypes.c_int, ctypes.c_double, ctypes.c_double, np.ctypeslib.ndpointer(dtype=np.float64), np.ctypeslib.ndpointer(dtype=np.float64), ctypes.c_double, ctypes.c_char_p]

    clip_lib.CliP(No_mutation, r, n, minor, total, ploidy,
                  Lambda_list, Lambda_num, alpha, rho, gamma, Run_limit, precision,
                  control_large, least_mut, post_th, least_diff,
                  coef, wcut, purity, preliminary_result.encode('utf-8'))
