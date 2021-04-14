import argparse
import sys
import os
from os import listdir
import json
import subprocess
import time

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--snv", type=str, required=True, help="Root path of the snv input")
parser.add_argument("-c", "--cn", type=str, required=True, help="Root path of the cna input")
parser.add_argument("-p", "--purity", type=str, required=True, help="Root path of the purity input")
parser.add_argument("-i", "--sample_id", type=str, default="sample", help="name of the sample being processed")
parser.add_argument("-e", "--preprocess_result", type=str, default="intermediate/", help="directory that stores the preprocess results")
parser.add_argument("-b", "--If_subsampling", action='store_true', help="whether doing subsampling or not")
parser.add_argument("-r", "--preliminary_result", type=str, default="Results_nosub/", help="directory that stores the preliminary results")

parser.add_argument("-f", "--final_result", type=str, default="final_result/", help="directory that stores the final results")
parser.add_argument("-g", "--filtering_flag", type=int, default = 1, help="whether filtering is needed during the postprocessing")
parser.add_argument("-l", "--Lambda", type=float, help="The penalty parameter")
parser.add_argument("-n", "--No_subsampling", type=int, help="The number of SNVs you want to include in each subsample.")
parser.add_argument("-m", "--Rep_num", type=int, help="The number of random subsamples needed.")
parser.add_argument("-w", "--window_size", type=float, default=0.05, help="controls the length of the interval between each cellular prevalence")
parser.add_argument("-o", "--overlap_size", type=float, default=0.0, help="controls how large two consecutive windows overlaps")


args_1 = parser.parse_args()


current_dir = os.path.dirname(os.path.abspath(__file__))
run_preprocess = current_dir + "/src/preprocess.R"

print(args_1.preliminary_result)

start = time.time()
# Run preprocess
p_preprocess = subprocess.Popen(["Rscript", run_preprocess, args_1.snv, args_1.cn, args_1.purity, args_1.sample_id, args_1.preprocess_result], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
p_preprocess.communicate()


run_CliP = current_dir + "/src/run_kernel_nosub.py"
python_clip = current_dir + "/src/kernel.py"
path_to_input_with_prefix = args_1.preprocess_result + args_1.sample_id
run_postprocess = current_dir + "/src/postprocess.R"
run_lambda_selection = current_dir + "/src/penalty_selection.py"

print(args_1.If_subsampling)

# Run the main CliP function (with subsampling)
if args_1.If_subsampling == False:
	if args_1.Lambda == None:
		p_run_CliP = subprocess.Popen(["python", run_CliP, path_to_input_with_prefix, args_1.preliminary_result, python_clip], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		p_run_CliP.communicate()
		
		# Run postprocess
		p_postprocess = subprocess.Popen(["Rscript", run_postprocess, args_1.preliminary_result, path_to_input_with_prefix, args_1.final_result, str(args_1.filtering_flag)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		p_postprocess.communicate()
		
		# The lambda selection methods:
		p_lambda_selection = subprocess.Popen(["python", run_lambda_selection, args_1.purity, args_1.final_result], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		p_lambda_selection.communicate()
		
	else:
		p_run_CliP = subprocess.Popen(["python", run_CliP, path_to_input_with_prefix, args_1.preliminary_result, python_clip, str(Lambda)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		p_run_CliP.communicate()
		
		# Run postprocess
		p_postprocess = subprocess.Popen(["Rscript", run_postprocess, args_1.preliminary_result, path_to_input_with_prefix, args_1.final_result, str(args_1.filtering_flag), str(args_1.Lambda)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		p_postprocess.communicate()

			
# Run the main CliP function (without subsampling)
else:
	subsampling_clip = current_dir + "/src/run_kernel_sub.py"
	if args_1.No_subsampling == None:
		sys.exit("Need an input for No_subsampling")
		
	
	if args_1.Rep_num == None:
		sys.exit("Need an input for Rep_num")
	
	if args_1.Lambda == None:
		p_run_subsampling = subprocess.Popen(["python", subsampling_clip, path_to_input_with_prefix, args_1.preliminary_result, python_clip, str(args_1.No_subsampling), str(args_1.Rep_num), str(args_1.window_size), str(args_1.overlap_size)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		p_run_subsampling.communicate()
		
		# Run postprocess
		p_postprocess = subprocess.Popen(["Rscript", run_postprocess, args_1.preliminary_result, path_to_input_with_prefix, args_1.final_result, str(args_1.filtering_flag)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		p_postprocess.communicate()
		
		# The lambda selection methods:
		p_lambda_selection = subprocess.Popen(["python", run_lambda_selection, args_1.purity, args_1.final_result], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		p_lambda_selection.communicate()
		
	else:
		p_run_subsampling = subprocess.Popen(["python", subsampling_clip, path_to_input_with_prefix, args_1.preliminary_result, python_clip, str(args_1.No_subsampling), str(args_1.Rep_num), str(args_1.window_size), str(args_1.overlap_size), str(args_1.Lambda)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		p_run_subsampling.communicate()
		
		# Run postprocess
		p_postprocess = subprocess.Popen(["Rscript", run_postprocess, args_1.preliminary_result, path_to_input_with_prefix, args_1.final_result, str(args_1.filtering_flag), str(args_1.Lambda)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		p_postprocess.communicate()
				
end = time.time()
print(" Time elapsed: ", end - start, "seconds")


