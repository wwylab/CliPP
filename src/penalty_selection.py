import os
import os.path
from os import path
import sys
import pandas as pd
import numpy as np
import shutil
from shutil import copyfile
from itertools import repeat
import operator 

Lambda_list = [0.01,0.03,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25]

purity_file = sys.argv[1]
final_result = sys.argv[2]

for i in range(11):
    datafile = pd.read_csv(os.path.join(final_result, "subclonal_structure_lam%s.txt" % (Lambda_list[i])), sep = "\t")
    max_cp_value = max(datafile["cellular_prevalence"])
    purity_value = pd.read_csv(purity_file, sep="\t", header=None)[0][0]
    #print(purity_value)
    source_1 = os.path.join(final_result, "subclonal_structure_lam%s.txt" % (Lambda_list[i]))
    source_2 = os.path.join(final_result, "mutation_assignments_lam%s.txt" % (Lambda_list[i]))
    destination = os.path.join(final_result, "Best_lambda/")
    if not os.path.exists(destination):
        os.makedirs(destination)
    if abs(max_cp_value - purity_value)/purity_value < 0.05:
        shutil.copy(source_1, destination) 
        shutil.copy(source_2, destination) 
        
existing_fnames = []
for i in range(1, 12):
    fname = os.path.join(destination, "subclonal_structure_lam%s.txt" % (Lambda_list[i-1]))
    try:
        f = open(fname, "r")
        f.close()
        existing_fnames.append(fname)
    except FileNotFoundError:
        #print("not found: " + fname)
        pass

# remove files
for i in range(len(existing_fnames) - 1):
    os.remove(existing_fnames[i])

if len(existing_fnames) == 0:
    lam_rev = [0.25,0.225,0.2,0.175,0.15,0.125,0.1,0.075,0.05,0.03,0.01]
    
    lst_value = []
    lst_index = []
    for i in range(11):
        datafile_1 = pd.read_csv(os.path.join(final_result, "subclonal_structure_lam%s.txt" % (lam_rev[i])), sep = "\t")
        max_cp_value = max(datafile_1["cellular_prevalence"])
        purity_value = pd.read_csv(purity_file, sep="\t", header=None)[0][0]
        #print(purity_value)
        value = abs(max_cp_value - purity_value)/purity_value
        lst_value.append(value)
    lst_pos = lst_value.index(min(lst_value))
    lam_pos = lam_rev[lst_pos]   
    source_1 = os.path.join(final_result, "subclonal_structure_lam%s.txt" % (lam_pos))
    source_2 = os.path.join(final_result, "mutation_assignments_lam%s.txt" % (lam_pos))
    destination = os.path.join(final_result, "Best_lambda/")
    shutil.copy(source_1, destination) 
    shutil.copy(source_2, destination) 
