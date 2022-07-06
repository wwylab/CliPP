import sys
import os

if not (sys.version_info[0] == 3 and sys.version_info[1] >= 5 and sys.version_info[2] >= 1):
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

print(_stderr)
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
  
try:
    shutil.rmtree(path_for_preliminary)
except:
    pass
    
print("Main CliP function finished.")
