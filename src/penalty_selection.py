#------------------------------------------------------------#
# Given the fact that there is a hyper-parameter lambda in our model, its selection is important for
# finding a balance between over- and under-fitting. The size of lambda controls the number of clusters. In
# general, higher lambda values return fewer clusters. We introduce a selection approach, which focuses on
# the match of clonal mutations, rather than a perfect match of the overall structure. 

# This script is a rule of thumb method to select one lambda when you run CliPP with multiple lambdas.
#------------------------------------------------------------#

#------------------------------------------------------------#
# The automated lambda selection pipeline is as follows: 
# 1) Run CliPP with 11 different lambdas: [0.01, 0.03, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25]. 
# 2) For each sample, compute a score: A = |(max(CP)−purity)/purity|. 
# If there are one or more results that satisfy A < 0.01, we select the largest lambda associated with those results. 
# (We observed A<0.01 in more than 70% of the PCAWG samples)
# If all scores A are greater than 0.01, we choose the lambda associated with the smallest A.
#------------------------------------------------------------#

#------------------------------------------------------------#
# Notation:
# CP: cellular prevalence
# CCF: cancer cell fraction, which is the proportion of tumor cells bearing the mutation
# Clonal mutation: mutations belonging to the initiating tumor cell and are expected to occur in every cell in the tumor
# Subclonal mutation: mutations which arose in descendant subpopulations
# Clonal fraction: number of clonal mutation divided by the number of total mutation
# Subclonal fraction: number of subclonal mutation divided by the number of total mutation
#------------------------------------------------------------#

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

# The 11 lambda values in the default lambda list
#Lambda_list = [0.01,0.03,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25]

def run_lambda_selection(purity_file, final_result):

    subclonal_structure_lam_map = dict()
    for f in os.listdir(final_result):
        if f.startswith('subclonal_structure_lam') and f.endswith('.txt'):
            _lambda = f.replace('subclonal_structure_lam', '').replace('.txt', '')
            _lambda = float(_lambda)

            subclonal_structure_lam_map[_lambda] = f

    if not subclonal_structure_lam_map:
        return

    purity_value = pd.read_csv(purity_file, sep="\t", header=None)[0][0]
    passed_lambda = []

    for _lambda in subclonal_structure_lam_map:
        f = subclonal_structure_lam_map[_lambda]
        datafile = pd.read_csv(os.path.join(final_result, f), sep = "\t")
        max_cp_value = max(datafile["cellular_prevalence"])
        
        # Filter 1: Drop the output with corresponding lambda value if the percent difference between its max CP value and input purity is >1% 
        if abs(max_cp_value - purity_value)/purity_value < 0.01:
            passed_lambda.append(_lambda)
            
    # If not all lambdas are dropped
    if passed_lambda:
        # Select the output given by the largest lambda in the current list
        selected_lambda = max(passed_lambda) 
        source_1 = os.path.join(final_result, "subclonal_structure_lam%s.txt" % (selected_lambda))
        source_2 = os.path.join(final_result, "mutation_assignments_lam%s.txt" % (selected_lambda))
        destination = os.path.join(final_result, "Best_lambda/")
        if not os.path.exists(destination):
            os.makedirs(destination)
        shutil.copy(source_1, destination) 
        shutil.copy(source_2, destination) 
    else:
        # If none of the lambdas pass filter 1
        lam_rev = subclonal_structure_lam_map.keys()
        lam_rev = list(lam_rev)
        lam_rev.sort(reverse=True)

        lst_value = []
        lst_index = []

        for _lambda in lam_rev:
            f = subclonal_structure_lam_map[_lambda]
            datafile_1 = pd.read_csv(os.path.join(final_result, f), sep = "\t")
            max_cp_value = max(datafile_1["cellular_prevalence"])
            value = abs(max_cp_value - purity_value)/purity_value
            lst_value.append(value)
            
        # Select the output given by the lambda that returns the smallest percent difference between its max CP value and input purity
        lst_pos = lst_value.index(min(lst_value))
        lam_pos = lam_rev[lst_pos]   
        source_1 = os.path.join(final_result, "subclonal_structure_lam%s.txt" % (lam_pos))
        source_2 = os.path.join(final_result, "mutation_assignments_lam%s.txt" % (lam_pos))
        destination = os.path.join(final_result, "Best_lambda/")
        if not os.path.exists(destination):
            os.makedirs(destination)
        shutil.copy(source_1, destination) 
        shutil.copy(source_2, destination) 
