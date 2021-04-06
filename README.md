# CliP
Clonal structure identification through penalizing pairwise differences

## Introduction
CliP is a subclonal identification tool designed for next-generation sequecing of bulk tumor samples. It is one of the 11  participating methods in the Pan-Cancer Analysis of Whole Genome (PCAWG) working group, Heterogeneity and Evolution, of International Cancer Genome Consortium (ICGC). And the method is described in this manuscript: https://www.biorxiv.org/content/10.1101/2021.03.31.437383v1.

## Prerequisitions
Need R(>3.3.1) and python(>3.5.1), the script does not support python2.

## Installing
There is no need to install CliP.

## Input data sample
Under /sample/. Those are example input data for CliP.


## Example code

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
After the preprocessing, you can run CliP as:
```
python run_CliP.py path_to_input_with_prefix path_to_output path_to_CliP lam
```
where, path_to_input_with_prefix is the path to the directory stores the preprocessed results, **along with whatever the prefix chosen above**; path_to_output denotes the directory where the CliP results should be put; path_to_clip is the path to CliP.py script; lam is still the lambda controls the penalization level, and usually takes values from 0.01-0.25.

In the example, you can now run
```
python run_CliP.py path/to/intermediate/sample path/to/Results_nosub/ CliP.py 0.2
```





