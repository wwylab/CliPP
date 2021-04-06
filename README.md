# CliP
Clonal structure identification through penalizing pairwise differences

## Introduction
Subpopulations of tumor cells characterized by mutation profiles may confer differential fitness and consequently influence prognosis of cancers. Understanding subclonal architecture has the potential to provide biological insight in tumor evolution and advance precision cancer treatment. Recent methods comprehensively integrate single nucleotide variants (SNVs) and copy number aberrations (CNAs) to reconstruct subclonal architecture using whole-genome or whole-exome sequencing (WGS, WES) data from bulk tumor samples. However, the commonly used Bayesian methods require a large amount of computational resources, a prior knowledge of the number of subclones, and extensive post-processing. Regularized likelihood modeling approach, never explored for subclonal reconstruction, can inherently address these drawbacks. We therefore propose a model-based method, Clonal structure identification through pair-wise Penalization, or CliP, for clustering subclonal mutations without prior knowledge or post-processing. The CliP model is applicable to genomic regions with or without CNAs. CliP demonstrates high accuracy in subclonal reconstruction through extensive simulation studies. A penalized likelihood framework for subclonal reconstruction will help address intrinsic drawbacks of existing methods and expand the scope of computational analysis for cancer evolution in large cancer genomic studies. Also see our preprint: https://www.biorxiv.org/content/10.1101/2021.03.31.437383v1.


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





