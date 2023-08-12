import sys
import os
from os import listdir

#if not (sys.version_info[0] == 3 and sys.version_info[1] >= 5 and sys.version_info[2] >= 1):
version_morph = sys.version_info[0]*10000+sys.version_info[1]*100+sys.version_info[2]
version_base = 30501
if not (version_morph >= version_base):
    sys.stderr.write("Error message: CliP can only run with python >=3.5.1\n")
    sys.exit(-1)

import argparse
import subprocess
import shutil
import time

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(current_dir, "src"))

from run_kernel_nosub import run_clip_nosub
from run_kernel_sub import run_clip_sub
from penalty_selection import run_lambda_selection

parser = argparse.ArgumentParser()

parser.add_argument("snv_input", type=str, help="Path of the snv input.")
parser.add_argument("cn_input", type=str, help="Path of the copy number input.")
parser.add_argument("purity_input", type=str, help="Path of the purity input.")
parser.add_argument("-i", "--sample_id", type=str, default="sample_id", help="Name of the sample being processed. Default is 'sample'.")
parser.add_argument("-b", "--subsampling", action='store_true', help="Whether doing subsampling or not. Default is not doing the subsampling, and a flag -b is needed when you want to do subsampling.")
parser.add_argument("-l", "--lam", type=float, help="The penalty parameter (lambda), which usually takes values from 0.01-0.25. If skipping this parameter, it will return a list of results that take value of [0.01, 0.03, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25] by default.")
parser.add_argument("-p", "--preprocess", type=str, default="preprocess_result/", help="Directory that stores the preprocess results. Default name is 'preprocess_result/'.")
parser.add_argument("-f", "--final", type=str, default="final_result/", help="Directory that stores the final results after postprocessing. Default name is 'final_result/'.")

parser.add_argument("-s", "--subsample_size", type=int, help="(Required if doing subsampling) The number of SNVs you want to include in each subsamples.")
parser.add_argument("-n", "--rep_num", type=int, help="(Required if doing subsampling) The number of random subsamples needed.")
parser.add_argument("-w", "--window_size", type=float, default=0.05, help="Controls the length of the window. Takes value between 0 and 1. Default is 0.05.")
parser.add_argument("-o", "--overlap_size", type=float, default=0.0, help="Controls the overlapped length of two consecutive windows. Takes value between 0 and 1. Default is 0.")

args = parser.parse_args()

run_preprocess = os.path.join(current_dir, "src/preprocess.R")
sample_id = args.sample_id
sample_id = sample_id.strip()
final_result = args.final.strip()

if sample_id == "":
	sys.stderr.write("User specified sample id is empty. Use default sample_id instead.\n")
	sample_id = "sample_id"

if final_result == "":
    sys.stderr.write("User specified final is empty. Use default final_result instead.\n")
    final_result = "final_result"

if args.subsampling:
    if args.subsample_size is None:
        sys.stderr.write("Please specify subsample_size\n")
        sys.exit(-1)
        
    if args.rep_num is None:
        sys.stderr.write("Please specify rep_num\n")
        sys.exit(-1)
        

lambda_list = [0.01, 0.03, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25]
if args.lam is not None:
    try:  
        lambda_list = [float(args.lam)]
    except Exception as err:
        sys.stderr.write(err)
        sys.exit(-1)

result_dir = os.path.join(current_dir, sample_id)
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

path_for_preprocess = os.path.join(result_dir, args.preprocess)
path_for_preliminary = os.path.join(result_dir, "preliminary_result")
path_for_final = os.path.join(result_dir, final_result)

# Run preprocessing
print("Running preprocessing...")
p_preprocess = subprocess.Popen(["Rscript", 
                                 run_preprocess, 
                                 args.snv_input, 
                                 args.cn_input, 
                                 args.purity_input, 
                                 args.sample_id, 
                                 path_for_preprocess], 
                                 stdout=subprocess.PIPE, 
                                 stderr=subprocess.PIPE)

_stdout, _stderr = p_preprocess.communicate()

if "error" in _stderr.decode().strip().lower():
    print(_stderr.decode().strip())
    sys.exit(-1)
print("Preprocessing finished.")

run_postprocess = os.path.join(current_dir, "src/postprocess.R")
# run_lambda_selection = os.path.join(current_dir, "src/penalty_selection.py")

# Run the main CliP function (without subsampling)
print("Running the main CliP function...")
if not args.subsampling:
    start = time.time()
    run_clip_nosub(path_for_preprocess, path_for_preliminary, lambda_list)
    end = time.time()
    elapsed_time = end - start
    print("\nElapsed time: %.2fsec" % elapsed_time + "\n")
	
    # Run postprocessing

    p_postprocess = subprocess.Popen(["Rscript", 
                                      run_postprocess, 
                                      path_for_preliminary, 
                                      path_for_preprocess, 
                                      path_for_final, 
                                      str(1)], 
                                     stdout=subprocess.PIPE, 
                                     stderr=subprocess.PIPE)
 
    _stdout, _stderr = p_postprocess.communicate()
        
    if "error" in _stderr.decode().strip().lower():
        print(_stderr.decode().strip())
        sys.exit(-1)
		
    # The lambda selection methods:
    if args.lam is None:
        run_lambda_selection(args.purity_input, path_for_final)

# Run the main CliP function (with subsampling)
else:    
    start = time.time()
    run_clip_sub(path_for_preprocess, 
                 path_for_preliminary, 
                 lambda_list,
                 args.subsample_size,
                 args.rep_num, 
                 args.window_size, 
                 args.overlap_size)
 
    end = time.time()
    elapsed_time = end - start
    print("\nElapsed time: %.2fsec" % elapsed_time + "\n")
    
    # Run postprocessing
    p_postprocess = subprocess.Popen(["Rscript", 
                                      run_postprocess, 
                                      path_for_preliminary, 
                                      path_for_preprocess, 
                                      path_for_final, str(1)], 
                                      stdout=subprocess.PIPE, 
                                      stderr=subprocess.PIPE)
    _stdout, _stderr = p_postprocess.communicate()
    if "error" in _stderr.decode().strip().lower():
        print(_stderr.decode().strip())
        sys.exit(-1)
        
    # The lambda selection methods:
    if args.lam is None:
        run_lambda_selection(args.purity_input, path_for_final)

## Annotate problematic SNVs. 

warning_tag = dict()
warning_tag[0] = []
warning_tag[1] = []

try:
    for f in listdir(path_for_preliminary):
        if f.endswith('_problematic_snvs.txt'):

            problematic_snv_indices = []
            with open(os.path.join(path_for_preliminary, f), 'r') as handle:
                for line in handle:
                    line = line.strip()
                    if len(line) > 0:
                        line = int(line)
                        problematic_snv_indices.append(line)

            lambda_num = f.replace('_problematic_snvs.txt', '')
            warning_tag[0].append(lambda_num.replace('lam', ''))
            
            mutation_assignment_f = 'mutation_assignments_%s.txt' % (lambda_num)
            mutation_assignments = []
            with open(os.path.join(path_for_final, mutation_assignment_f), 'r') as handle:
                for line in handle:
                    line = line.strip()
                    if len(line) > 0:
                        mutation_assignments.append(line)

            output_handle = open(os.path.join(path_for_final, mutation_assignment_f), 'w')
            count = 0

            for line in mutation_assignments:
                if line.startswith('chromosome_index'):
                    line_out = line + '\t' + 'Warning'
                    output_handle.write(line_out + '\n')
                else:
                    if count in problematic_snv_indices:
                        line_out = line + '\t' + '1'
                    else:
                        line_out = line + '\t' + '0'
                    output_handle.write(line_out + '\n')
                    count += 1
                    
            
    for f in listdir(os.path.join(path_for_final, 'Best_lambda')):
        if f.startswith('mutation_assignments_'):
            shutil.copy(os.path.join(path_for_final, f), os.path.join(path_for_final, 'Best_lambda', f))
        

    shutil.rmtree(path_for_preliminary)
except:
    pass

## Write warning file to final_result/Best_lambda folder if applicable

lambda_list_with_results = []
lambda_list_without_results = []
for f in listdir(path_for_final):
    if f.startswith('mutation_assignments_lam') and f.endswith('.txt'):
        lambda_num = f.replace('mutation_assignments_lam', '').replace('.txt', '')
        lambda_list_with_results.append(float(lambda_num))

for lambda_num in lambda_list:
    if lambda_num not in lambda_list_with_results:
        lambda_list_without_results.append(lambda_num)

warning_tag[1] = lambda_list_without_results

if warning_tag[0] or warning_tag[1]:
    if not os.path.exists(os.path.join(path_for_final, 'Best_lambda')):
        os.mkdir(os.path.join(path_for_final, 'Best_lambda'))

    output_handle = open(os.path.join(path_for_final, 'Best_lambda/WARNING.txt'), 'w')
    output_handle.write('This sample is problematic due to the reason(s) below. Please take caution when use the CliP result.\n')
    
    if warning_tag[0]:
        best_lambda_tag = 0
        for f in listdir(os.path.join(path_for_final, 'Best_lambda')):
            if f.startswith('mutation_assignments_lam') and f.endswith('.txt'):
                with open(os.path.join(path_for_final, 'Best_lambda', f)) as handle:
                    for line in handle:
                        line = line.strip().split()
                        if len(line) == 4:
                            if line[-1] == '1':
                                best_lambda_tag = 1
                                break

        
        _lambdas = ', '.join(warning_tag[0])
        _mutation_assignment = ['mutation_assignments_lam%s.txt' % (_lambda) for _lambda in warning_tag[0]]
        _mutation_assignment = ', '.join(_mutation_assignment)
            
        output_string = '1. The clustering assignment for some of the SNVs may be not correct for these lambdas: %s. Please check the Warning column in the %s file(s) in the final_result folder' % (_lambdas, _mutation_assignment)

        if best_lambda_tag:
            output_string = output_string + '. The mutation_assignments file in the Best_lambda folder also has such SNVs. Please check the Warning column as well.'
        
        output_handle.write(output_string + '\n')

    if warning_tag[1]:
        if warning_tag[0]:
            output_handle.write('2. ')
        else:
            output_handle.write('1. ')

        _lambdas = [str(_lambda) for _lambda in  warning_tag[1]]
        _lambdas = ', '.join(_lambdas)
        output_handle.write('These lambdas do not have mutation clustering results: %s. The final CliP result in the Best_lambda folder is selected based on the lambdas whose results are available in the final_result folder.\n' % (_lambdas))

        
print("Main CliP function finished.")
