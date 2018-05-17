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
Rscript preprocess.R Sample_data/sample_vcf.txt Sample_data/sample_cnv1.txt Sample_data/intermediate/ meta.Rdata
```

### Running CliP
After the preprocess was done, you can run CliP as:
```
python3 run_CliP.py path_to_input path_to_output path_to_clip lam
```
where 'lam' usually takes values from 0.5-1.5.

In our data, you can now run
```
python3 run_CliP.py Sample_data/intermediate/ Sample_data/Results/ CliP.py 0.5
```

CliP is limited by memory, if you do not have enough memory to process all SNVs, we do supply a downsampling strategy, one may run 
```
python3 run_CliP_subsampling.py path_to_input path_to_output path_to_clip lam No_subsampling Rep_id window_size overlap_size
```
where, No_subsampling quantifies the number of SNVs you want to include, Rep_id serves as a random seed, and typically we do a couple of runs if a down sampling is needed. Window_size actually 










