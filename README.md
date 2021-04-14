# CliP
Clonal structure identification through penalizing pairwise differences

## Introduction
Subpopulations of tumor cells characterized by mutation profiles may confer differential fitness and consequently influence prognosis of cancers. Understanding subclonal architecture has the potential to provide biological insight in tumor evolution and advance precision cancer treatment. Recent methods comprehensively integrate single nucleotide variants (SNVs) and copy number aberrations (CNAs) to reconstruct subclonal architecture using whole-genome or whole-exome sequencing (WGS, WES) data from bulk tumor samples. However, the commonly used Bayesian methods require a large amount of computational resources, a prior knowledge of the number of subclones, and extensive post-processing. Regularized likelihood modeling approach, never explored for subclonal reconstruction, can inherently address these drawbacks. We therefore propose a model-based method, Clonal structure identification through pair-wise Penalization, or CliP, for clustering subclonal mutations without prior knowledge or post-processing. The CliP model is applicable to genomic regions with or without CNAs. CliP demonstrates high accuracy in subclonal reconstruction through extensive simulation studies. A penalized likelihood framework for subclonal reconstruction will help address intrinsic drawbacks of existing methods and expand the scope of computational analysis for cancer evolution in large cancer genomic studies. Also see our paper: https://www.biorxiv.org/content/10.1101/2021.03.31.437383v1.


## Prerequisitions
Need R(>3.3.1) and python(>3.5.1), the script does not support python2.

## Installing
There is no need to install CliP.

## Structure of CliP implementation
The flow chart below shows the CliP implementation. Raw functions and scripts are under ./scr/. 
![image](https://user-images.githubusercontent.com/14543452/114482762-bf4c1480-9bcc-11eb-8c96-a944611e91d7.png)

CliP can run on samples with up to 50,000 SNVs on a machine with 256GB memory. It requires more memory when there are more SNVs. When that happens, we apply a downsampling strategy, as implemented in all other subclonal reconstruction methods. 

## Input data sample
There are three input files:

1. ```sample.snv.txt```: A tab separated file which stores a data matrix. 
* ```chromosome_index```: The chromosomal location of the SNV.
* ```position```: the single-nucleotide position of the SNV on the corresponding chromosome.
* ```ref_count```: The number of reads covering the locus and containing the reference allele.
* ```alt_count```: The number of reads covering the locus and containing the alternative allele.

2. ```sample.cna.txt```: A tab separated file which stores a data matrix.
* ```chromosome_index```: The chromosome location of the CNA.
* ```start_position```: the start position of the CNA segment on the corresponding chromosome.
* ```end_position```: the end position of the CNA segment on the corresponding chromosome.
* ```major_cn```: The copy number of the major allele in tumor cells. This should be greater than equal to the value in the minor_cn column and greater than 0.
* ```minor_cn```: The copy number of the minor allele. This must be less than equal the value in the major_cn column.
* ```total_cn```: The sum of major_cn and minor_cn.

3. ```sample.purity```: A scalar between 0 and 1.

A simulated sample input data is under ./sample/. 


## Example code (one-step implementation)

The caller function ``clip_main.py'' wraps up the CliP pipeline and ensures users to implement subclonal reconstruction in one-step. A full version of the command is as follows:

```
python run_clip_main.py -s snv -c cn -p purity -i sample_id -e preprocess_result -b -r preliminary_result -f final_result -g filtering_flag -l Lambda -n No_subsampling -m Rep_num -w window_size -o overlap_size
```
Here are the details of the options:

* ```snv```: A mandatory argument. Root path of the snv input.
* ```cn```: A mandatory argument. Root path of the cna input.
* ```purity```: A mandatory argument. Root path of the purity input.
* ```sample_id```: The name of the sample being processed. Default is ```sample```.
* ```preprocess_result```: Directory that stores the preprocess results. Default name is ```intermediate/```.
* ```If_subsampling```: Whether doing subsampling or not. Default is ```False```, and a flag ```-b``` is needed when you want to do subsampling
* ```preliminary_result```: Directory that stores the output of the kernel function, which is considered as the preliminary results. Default name is ```Result_nosub/``` (no_subsampling).
* ```final_result```: Directory that stores the final results after postprocessing. Default name is ```final_result/```.
* ```filtering_flag```: Whether filtering is needed during the postprocessing. Take value of 0 or 1. Default is 1 (need filtering).
* ```Lambda```: The penalty parameter, which usually takes values from 0.01-0.25. If skipping this parameter, it will return a list of results that take value of [0.01, 0.03, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25] by default.

The followings are parameters only needed when doing subsampling. We take partitions from 0-1 (determined by users through windows, default at 0.05), then take the VAF of each SNV as an initial estimation of its CP value and randomly assign SNVs that belong to their corresponding windows. The number of sampled SNV is proportional to number of total SNVs belongs to each window. 
* ```No_subsampling```: Required if doing subsampling. The number of SNVs you want to include in each subsamples.
* ```Rep_num```: Required if doing subsampling. The number of random subsamples needed.
* ```window_size```: Controls the length of the window. Takes value between 0 and 1. Default is ```0.05```.
* ```overlap_size```: Controls the overlapped length of two consecutive windows. Takes value between 0 and 1. Default is ```0```.
