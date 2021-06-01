import argparse
import sys
import os
from os import listdir
import json
import subprocess
import shutil
import time
import threading

from src.run_kernel_nosub import run_clip_nosub
from src.run_kernel_sub import run_clip_sub

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

current_dir = os.path.dirname(os.path.abspath(__file__))
run_preprocess = os.path.join(current_dir, "src/preprocess.R")

result_dir = os.path.join(current_dir, args.sample_id)
if not os.path.exists(result_dir):
	os.makedirs(result_dir)

path_for_preprocess = os.path.join(result_dir, args.preprocess)
path_for_preliminary = os.path.join(result_dir, "preliminary_result")
path_for_final = os.path.join(result_dir, args.final)

# Run preprocess
print("Running preprocessing...")
p_preprocess = subprocess.Popen(["Rscript", run_preprocess, args.snv_input, args.cn_input, args.purity_input, args.sample_id, path_for_preprocess], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
_stdout, _stderr = p_preprocess.communicate()
if _stderr:
	print(_stderr.decode().strip())
	sys.exit()
print("Preprocessing finished.")

run_CliP = os.path.join(current_dir, "src/run_kernel_nosub.py")
python_clip = os.path.join(current_dir, "src/kernel.py")
run_postprocess = os.path.join(current_dir, "src/postprocess.R")
run_lambda_selection = os.path.join(current_dir, "src/penalty_selection.py")

# Run the main CliP function (without subsampling)
print("Running the main CliP function...")
if args.subsampling == False:
	if args.lam == None:
		start = time.time()
		t = threading.Thread(name="Running the main CliP function", target=run_clip_nosub, args=(path_for_preprocess, path_for_preliminary))
		t.start()
		t.join()
		end = time.time()
		elapsed_time = end - start
		print("\nElapsed time: %.2fsec" % elapsed_time + "\n")
		
		# Run postprocess
		p_postprocess = subprocess.Popen(["Rscript", run_postprocess, path_for_preliminary, path_for_preprocess, path_for_final, str(1)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		_stdout, _stderr = p_postprocess.communicate()
		if _stderr:
			print(_stderr.decode().strip())
			sys.exit()
		
		# The lambda selection methods:
		p_lambda_selection = subprocess.Popen(["python", run_lambda_selection, args.purity_input, path_for_final], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		_stdout, _stderr = p_lambda_selection.communicate()
		if _stderr:
			print(_stderr.decode().strip())
			sys.exit()
		
	else:
		p_run_CliP = subprocess.Popen(["python", run_CliP, path_for_preprocess, path_for_preliminary, python_clip, str(args.lam)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		_stdout, _stderr = p_run_CliP.communicate()
		if _stderr:
			print(_stderr.decode().strip())
			sys.exit()
		
		# Run postprocess
		p_postprocess = subprocess.Popen(["Rscript", run_postprocess, path_for_preliminary, path_for_preprocess, path_for_final, str(1), str(args.lam)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		_stdout, _stderr = p_postprocess.communicate()
		if _stderr:
			print(_stderr.decode().strip())
			sys.exit()

			
# Run the main CliP function (with subsampling)
else:
	subsampling_clip = os.path.join(current_dir, "src/run_kernel_sub.py")
	if args.subsample_size == None:
		sys.exit("Need an input for subsample_size")
		
	if args.rep_num == None:
		sys.exit("Need an input for rep_num")
	
	if args.lam == None:

		start = time.time()
		t = threading.Thread(name="Running the main CliP function", target=run_clip_sub, args=(path_for_preprocess,path_for_preliminary,python_clip, args.subsample_size,
			args.rep_num, args.window_size, args.overlap_size))
		
		t.start()
		t.join()
		end = time.time()
		elapsed_time = end - start
		print("\nElapsed time: %.2fsec" % elapsed_time + "\n")
		
		# Run postprocess
		p_postprocess = subprocess.Popen(["Rscript", run_postprocess, path_for_preliminary, path_for_preprocess, path_for_final, str(1)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		_stdout, _stderr = p_postprocess.communicate()
		if _stderr:
			print(_stderr.decode().strip())
			sys.exit()
		
		# The lambda selection methods:
		p_lambda_selection = subprocess.Popen(["python", run_lambda_selection, args.purity_input, path_for_final], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		_stdout, _stderr = p_lambda_selection.communicate()
		if _stderr:
			print(_stderr.decode().strip())
			sys.exit()
		
	else:
		p_run_subsampling = subprocess.Popen(["python", subsampling_clip, path_for_preprocess, path_for_preliminary, python_clip, str(args.subsample_size), str(args.rep_num), str(args.window_size), str(args.overlap_size), str(args.lam)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		_stdout, _stderr = p_run_subsampling.communicate()
		if _stderr:
			print(_stderr.decode().strip())
			sys.exit()
		
		# Run postprocess
		p_postprocess = subprocess.Popen(["Rscript", run_postprocess, path_for_preliminary, path_for_preprocess, path_for_final, str(1), str(args.lam)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		_stdout, _stderr = p_postprocess.communicate()
		if _stderr:
			print(_stderr.decode().strip())
			sys.exit()

shutil.rmtree(path_for_preliminary)

print("Main CliP function finished.")


