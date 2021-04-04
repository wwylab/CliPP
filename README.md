# CliP
Clonal structure identification through penalizing pairwise differences

## Introduction
CliP is a subclonal identification tool designed for next-generation sequecing of bulk tumor samples. It is one of the 11  participating methods in the Pan-Cancer Analysis of Whole Genome (PCAWG) working group, Heterogeneity and Evolution, of International Cancer Genome Consortium (ICGC). And the method is described in this manuscript (Cite CliP paper).

## Prerequisitions
Need R(>3.3.1) and python(>3.5.1), the script does not support python2.

## Installing
There is no need to install CliP.

## Input data sample
Under /sample/. Those are simulated input for CliP.

## Input data format
CliP need 3 input files (please see sample for a more visualized example of input data):

1. SNV file: a tab separated file containing 8 columns. Specifically, it needs to include chromosome information, the position of the SNV, alt read, and ref read information. Please see the sample input for the detailed format. 

2. CNA file: a tab separated file containing 6 columns. The first column denotes the chromosome, the second one is the start position of the CNV segment, the third column records the end position of the CNV segment. The last three columns takes the actual copy number of each segment, which are major & minor allele total and total copy number.

3. Purity file: A file containing the purity either an estimation from CNA or a guess that will be corrected hopefully by CliP. For details on purity issue, please check our manuscript. 

Input files from the same sample need to be saved in the same directory with the same prefix. Specifically, they need to be named as
```
prefix.vcf prefix.cna.txt prefix.purity
```

## Start Your First Example

### Preprocess 
We supply a preprocess script preprocess.R to prepare the actual input files for CliP.
```
Rscript preprocess.R input_SNV input_CNV purity_file Input_prefix Output_dir
```

In the example you may run
```
Rscript preprocess.R path/to/sample/sample.vcf path/to/sample/sample.cna.txt path/to/sample/sample.purity sample path/to/intermediate
```

### Running CliP
After the preprocess was done, you can run CliP as:
```
python run_CliP.py path_to_input_with_prefix path_to_output path_to_CliP lam
```
where, path_to_input_with_prefix is the path to the directory stores the preprocessed results, **along with whatever the prefix chosen above**; path_to_output denotes the directory where the CliP results should be put; path_to_clip is the path to CliP.py script; lam is still the lambda controls the penalization level, and usually takes values from 0.01-0.25.

In our data, you can now run
```
python run_CliP.py path/to/intermediate/sample path/to/Results_nosub/ CliP.py 0.2
```

CliP is limited by memory (a good estimates, if you do not have enough memory to process all SNVs, we do supply a downsampling strategy, one may run 
```
python run_CliP_subsampling.py path_to_input_with_prefix path_to_output path_to_clip lam No_subsampling Rep_num window_size overlap_size
```
where, path_to_input_with_prefix is the path to the directory stores the preprocessed results, **along with whatever the prefix chosen above**; path_to_output denotes the directory where the CliP results should be put; path_to_clip is the path to CliP.py script; and lam is still the lambda controls the penalization level; No_subsampling quantifies the number of SNVs you want to include; Rep_num indicates the number of random subsamples needed.

The sampling is done for each interval of cellular prevalence, the sampled SNV is proportional to number of total SNVs belong to each interval. The window_size argument controls the length of the interval, overlap_size controls how large two consecutive windows overlaps. 

In the example we have, you may run
```
python run_CliP_subsampling.py path/to/intermediate/sample path/to/Results_sub/ CliP.py 0.2 200 5 0.05 0

```
to create 5 subsamples.

### Postprocess
If you run CliP directly, typically you do not need to run postprocess. When the average read depth is low (e.g. ~30X) you may want to run the postprocess with filtering switched on, which in this serves as a denoise step. If you ran CliP with downsampling, you have to run the postprocess script to obtain the final results. You can run post process as
```
Rscript post_analysis_run.R CliP_results_dir preprocessed_output_dir_with_prefix Output_dir lambda filtering_flag
```
Here CliP_results_dir is the directory stores the direct output from CliP; preprocessed_output_dir_with_prefix is the directory stores the preprocessed data, **along with whatever the prefix chosen above**; Output_dir records the directory of final output; lambda sepecifies which run you want to postprocess; 
filtering_flag taking values 0 and 1，meaning no filtering and filtering, respectively. 

In our non-downsampled run, you can run the following:
```
Rscript post_analysis_run.R path/to/Results_nosub/ path/to/intermediate/sample path/to/Final_res/ 0.2 1
```
and the downsampled version:
```
Rscript post_analysis_run.R path/to/Results_sub/ path/to/intermediate/sample path/to/Final_res/ 0.2 1
```
***IMPORTANT if you put both downsampled and non-downsampled results in the same directory, the non-downsampled is always picked! ***





