# CliP
Clonal structure identification through penalizing pairwise differences

## Introduction
CliP is a subclonal identification tool designed for next-generation sequecing of bulk tumor samples. It is one of the 11  participating methods in the Pan-Cancer Analysis of Whole Genome (PCAWG) working group, Heterogeneity and Evolution, of International Cancer Genome Consortium (ICGC). And the method is described in this manuscript (Cite CliP paper).

## Prerequisitions
Need R(>3.3.1) and python(>3.5.1), the script does not support python2.

## Installing
There is no need to install CliP.
## Input data format
CliP need 3 input files (please see Sample data for a more visualized example of input data):

1. SNV file: a tab separated file containing 4 columns, the first column denotes the chromosome, the second one is the position of the SNV, the third column records alt read, and the last column denotes the ref read.

2. CNV file: a tab separated file containing 5 or 4 columns, the first column denotes the chromosome, the second one is the start position of the CNV segment, the third column records the end position of the CNV segment. For the actual copy number of each segment, CliP accepts both total-only or allele specific copy number.

3. Purity file: A file containing the purity either an estimation from CNA or a guess that will be corrected hopefully by CliP. For details on purity issue, please check our manuscript. 

## Start Your First Example

### Preprocess 
We supply a preprocess script preprocess.R to prepare the actual input files for CliP.
```
Rscript preprocess.R input_SNV input_CNV purity_file Output_dir meta_file
```
The meta file should end with '.Rdata'

In the example we have you may run
```
Rscript preprocess.R Sample_data/sample_vcf.txt Sample_data/sample_cnv1.txt Sample_data/sample_purity.txt Sample_data/intermediate/meta.Rdata
```

### Running CliP
After the preprocess was done, you can run CliP as:
```
python3 run_CliP.py path_to_input path_to_output path_to_clip lam
```
where, path_to_input is the path to the directory stores the preprocessed results, path_to_output denotes the directory where the CliP reaults should be put, path_to_clip is the path to CliP.py script, and lam is still the lambda controls the penalization level, and usually takes values from 0.1-1.5.

In our data, you can now run
```
python3 run_CliP.py Sample_data/intermediate/meta.Rdata/ Sample_data/Results_nosub/ CliP.py 0.2
```

CliP is limited by memory (a good estimates , if you do not have enough memory to process all SNVs, we do supply a downsampling strategy, one may run 
```
python3 run_CliP_subsampling.py path_to_input path_to_output path_to_clip lam No_subsampling Rep_id window_size overlap_size
```
where, path_to_input is the path to the directory stores the preprocessed results, path_to_output denotes the directory where the CliP reaults should be put, path_to_clip is the path to CliP.py script, and lam is still the lambda controls the penalization level, No_subsampling quantifies the number of SNVs you want to include, Rep_id serves as a random seed, and typically we do 5-10 runs if a down sampling is needed. 

The sampling is done for each interval of cellular prevalence, the sampled SNV is proportional to number of total SNVs belong to each interval. The window_size argument controls the length of the interval, overlap_size controls how large two consecutive windows overlaps. 

In the example we have, you may run
```
python3 run_CliP_subsampling.py Sample_data/intermediate/ Sample_data/Results_sub/ CliP.py 0.2 150 1 0.05 0
python3 run_CliP_subsampling.py Sample_data/intermediate/ Sample_data/Results_sub/ CliP.py 0.2 150 2 0.05 0
python3 run_CliP_subsampling.py Sample_data/intermediate/ Sample_data/Results_sub/ CliP.py 0.2 150 3 0.05 0
```
to create 3 subsamples.
### Postprocess
If you run CliP directly, typically you do not need to run postprocess. When the average read depth is low (e.g. ~30X) you may want to run the postprocess with filtering switched on, which in this serves as a denoise step. If you ran CliP with downsampling, you have to run the postprocess script to obtain the final results. You can run post process as
```
Rscript post_analysis_run.R CliP_results_dir preprocessed_output_dir Output_dir lambda filtering_flag
```
Here CliP_results_dir is the directory stores the direct output from CliP, preprocessed_output_dir is the directory stores the preprocessed data, Output_dir records the directory of final output, lambda sepecifies which run you want to postprocess, and 
filtering_flag taking values 0 and 1， meaning no filtering and filtering, respectively. 

In our non-downsampled run, you can run the following:
```
Rscript post_analysis_run.R Sample_data/Results_nosub/ Sample_data/intermediate/ Sample_data/Final_res/ 0.2 1
```
and the downsampled version:
```
Rscript post_analysis_run.R Sample_data/Results_sub/ Sample_data/intermediate/ Sample_data/Final_res/ 0.2 1
```
***IMPORTANT if you put both downsampled and non-downsampled results in the same directory, the non-downsampled is always picked! ***





