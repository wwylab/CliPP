import argparse
import sys
import os
from os import listdir
import json
import subprocess
import time

parser = argparse.ArgumentParser()
#parser.add_argument("-s", "--snv", type=str, required=True, help="Root path of the snv input")
#parser.add_argument("-c", "--cn", type=str, required=True, help="Root path of the cna input")
#parser.add_argument("-p", "--purity", type=str, required=True, help="Root path of the purity input")
parser.add_argument("snv_input", type=str, help="Path of the snv input.")
parser.add_argument("cn_input", type=str, help="Path of the copy number input.")
parser.add_argument("purity_input", type=str, help="Path of the purity input.")
parser.add_argument("-i", "--sample_id", type=str, default="sample", help="Name of the sample being processed. Default is 'sample'.")
parser.add_argument("-p", "--preprocess", type=str, default="intermediate/", help="Directory that stores the preprocess results. Default name is 'intermediate/'.")
parser.add_argument("-b", "--subsampling", action='store_true', help="Whether doing subsampling or not. Default is not doing the subsampling, and a flag -b is needed when you want to do subsampling.")
parser.add_argument("-r", "--preliminary", type=str, default="preliminary_result/", help="Directory that stores the output of the kernel function, which is considered as the preliminary results. Default name is 'preliminary_result/'.")

<<<<<<< HEAD
parser.add_argument("-f", "--final", type=str, default="final_result/", help="Directory that stores the final results after postprocessing. Default name is 'final_result/'.")
parser.add_argument("-nf", "--no_filtering", action='store_false', help="If filtering is not wanted. Default is doing the filtering, and a flag -nf is needed when you don't want to do the filtering.")
parser.add_argument("-l", "--Lambda", type=float, help="The penalty parameter, which usually takes values from 0.01-0.25. If skipping this parameter, it will return a list of results that take value of [0.01, 0.03, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25] by default.")
parser.add_argument("-s", "--subsample_size", type=int, help="(Required if doing subsampling) The number of SNVs you want to include in each subsamples.")
parser.add_argument("-n", "--rep_num", type=int, help="(Required if doing subsampling) The number of random subsamples needed.")
parser.add_argument("-w", "--window_size", type=float, default=0.05, help="Controls the length of the window. Takes value between 0 and 1. Default is 0.05.")
parser.add_argument("-o", "--overlap_size", type=float, default=0.0, help="Controls the overlapped length of two consecutive windows. Takes value between 0 and 1. Default is 0.")
=======
parser.add_argument("-f", "--final_result", type=str, default="final_result/", help="directory that stores the final results")
parser.add_argument("-g", "--filtering_flag", type=int, default = 1, help="whether filtering is needed during the postprocessing")
parser.add_argument("-l", "--Lambda", type=float, help="The penalty parameter")
parser.add_argument("-n", "--No_subsampling", type=int, help="The number of SNVs you want to include in each subsample.")
parser.add_argument("-m", "--Rep_num", type=int, help="The number of random subsamples needed.")
parser.add_argument("-w", "--window_size", type=float, default=0.05, help="controls the length of the interval between each cellular prevalence")
parser.add_argument("-o", "--overlap_size", type=float, default=0.0, help="controls how large two consecutive windows overlaps")
>>>>>>> 142e174647f1677825fbe9b1849d4aaab8b9468e


args_1 = parser.parse_args()

current_dir = os.path.dirname(os.path.abspath(__file__))
run_preprocess = current_dir + "/src/preprocess.R"

print(int(args_1.no_filtering))

start = time.time()
# Run preprocess
p_preprocess = subprocess.Popen(["Rscript", run_preprocess, args_1.snv_input, args_1.cn_input, args_1.purity_input, args_1.sample_id, args_1.preprocess], stdout=subprocess.PIPE, stderr=sys.stdout.buffer)
p_preprocess.communicate()


run_CliP = current_dir + "/src/run_kernel_nosub.py"
python_clip = current_dir + "/src/kernel.py"
path_to_input_with_prefix = args_1.preprocess + args_1.sample_id
run_postprocess = current_dir + "/src/postprocess.R"
run_lambda_selection = current_dir + "/src/penalty_selection.py"

print(args_1.subsampling)

# Run the main CliP function (with subsampling)
if args_1.subsampling == False:
	if args_1.Lambda == None:
		p_run_CliP = subprocess.Popen(["python", run_CliP, path_to_input_with_prefix, args_1.preliminary, python_clip], stdout=subprocess.PIPE, stderr=sys.stdout.buffer)
		p_run_CliP.communicate()
		
		# Run postprocess
		p_postprocess = subprocess.Popen(["Rscript", run_postprocess, args_1.preliminary, path_to_input_with_prefix, args_1.final, str(int(args_1.no_filtering))], stdout=subprocess.PIPE, stderr=sys.stdout.buffer)
		p_postprocess.communicate()
		
		# The lambda selection methods:
		p_lambda_selection = subprocess.Popen(["python", run_lambda_selection, args_1.purity_input, args_1.final], stdout=subprocess.PIPE, stderr=sys.stdout.buffer)
		p_lambda_selection.communicate()
		
	else:
		p_run_CliP = subprocess.Popen(["python", run_CliP, path_to_input_with_prefix, args_1.preliminary, python_clip, str(Lambda)], stdout=subprocess.PIPE, stderr=sys.stdout.buffer)
		p_run_CliP.communicate()
		
		# Run postprocess
		p_postprocess = subprocess.Popen(["Rscript", run_postprocess, args_1.preliminary, path_to_input_with_prefix, args_1.final, str(int(args_1.no_filtering)), str(args_1.Lambda)], stdout=subprocess.PIPE, stderr=sys.stdout.buffer)
		p_postprocess.communicate()

			
# Run the main CliP function (without subsampling)
else:
	subsampling_clip = current_dir + "/src/run_kernel_sub.py"
	if args_1.subsample_size == None:
		sys.exit("Need an input for subsample_size")
		
	
	if args_1.rep_num == None:
		sys.exit("Need an input for rep_num")
	
	if args_1.Lambda == None:
		p_run_subsampling = subprocess.Popen(["python", subsampling_clip, path_to_input_with_prefix, args_1.preliminary, python_clip, str(args_1.subsample_size), str(args_1.rep_num), str(args_1.window_size), str(args_1.overlap_size)], stdout=subprocess.PIPE, stderr=sys.stdout.buffer)
		p_run_subsampling.communicate()
		
		# Run postprocess
		p_postprocess = subprocess.Popen(["Rscript", run_postprocess, args_1.preliminary, path_to_input_with_prefix, args_1.final, str(int(args_1.no_filtering))], stdout=subprocess.PIPE, stderr=sys.stdout.buffer)
		p_postprocess.communicate()
		
		# The lambda selection methods:
		p_lambda_selection = subprocess.Popen(["python", run_lambda_selection, args_1.purity_input, args_1.final], stdout=subprocess.PIPE, stderr=sys.stdout.buffer)
		p_lambda_selection.communicate()
		
	else:
		p_run_subsampling = subprocess.Popen(["python", subsampling_clip, path_to_input_with_prefix, args_1.preliminary, python_clip, str(args_1.subsample_size), str(args_1.rep_num), str(args_1.window_size), str(args_1.overlap_size), str(args_1.Lambda)], stdout=subprocess.PIPE, stderr=sys.stdout.buffer)
		p_run_subsampling.communicate()
		
		# Run postprocess
		p_postprocess = subprocess.Popen(["Rscript", run_postprocess, args_1.preliminary, path_to_input_with_prefix, args_1.final, str(int(args_1.no_filtering)), str(args_1.Lambda)], stdout=subprocess.PIPE, stderr=sys.stdout.buffer)
		p_postprocess.communicate()
				
end = time.time()
print(" Time elapsed: ", end - start, "seconds")


