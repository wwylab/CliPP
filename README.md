# CliP
Clonal structure identification through penalizing pairwise differences

## Introduction
Subpopulations of tumor cells characterized by mutation profiles may confer differential fitness and consequently influence prognosis of cancers. Understanding subclonal architecture has the potential to provide biological insight in tumor evolution and advance precision cancer treatment. Recent methods comprehensively integrate single nucleotide variants (SNVs) and copy number aberrations (CNAs) to reconstruct subclonal architecture using whole-genome or whole-exome sequencing (WGS, WES) data from bulk tumor samples. However, the commonly used Bayesian methods require a large amount of computational resources, a prior knowledge of the number of subclones, and extensive post-processing. Regularized likelihood modeling approach, never explored for subclonal reconstruction, can inherently address these drawbacks. We therefore propose a model-based method, Clonal structure identification through pair-wise Penalization, or CliP, for clustering subclonal mutations without prior knowledge or post-processing. The CliP model is applicable to genomic regions with or without CNAs. CliP demonstrates high accuracy in subclonal reconstruction through extensive simulation studies. A penalized likelihood framework for subclonal reconstruction will help address intrinsic drawbacks of existing methods and expand the scope of computational analysis for cancer evolution in large cancer genomic studies. Also see our paper: https://www.biorxiv.org/content/10.1101/2021.03.31.437383v1.


## Prerequisites
- R [>3.3.1]
- python [>3.5.1]
- NumPy
- SciPy
- pandas

## Setting up CliP
```
git clone https://github.com/wwylab/CliP.git
cd CliP
```

## Structure of CliP implementation
The flow chart below shows the CliP implementation. Raw functions and scripts are under `scr/`. 
![image](https://user-images.githubusercontent.com/14543452/114482762-bf4c1480-9bcc-11eb-8c96-a944611e91d7.png)

CliP can run on samples with up to 50,000 SNVs on a machine with 256GB memory. It may require more memory when there are more SNVs. When that happens, we apply a downsampling strategy, as implemented in all other subclonal reconstruction methods. 

Furthermore, CliP is automatically run in parallel when users have multiple cores available.

## Input data sample
There are three required input files:

1. ```sample.snv.txt```: A tab separated file which stores a data matrix with the following named columns:
* ```chromosome_index```: The chromosomal location of the SNV.
* ```position```: the single-nucleotide position of the SNV on the corresponding chromosome.
* ```ref_count```: The number of reads covering the locus and containing the reference allele.
* ```alt_count```: The number of reads covering the locus and containing the alternative allele.

2. ```sample.cna.txt```: A tab separated file which stores a data matrix with the following named columns:
* ```chromosome_index```: The chromosome location of the CNA.
* ```start_position```: the start position of the CNA segment on the corresponding chromosome.
* ```end_position```: the end position of the CNA segment on the corresponding chromosome.
* ```major_cn```: The copy number of the major allele in tumor cells. This should be greater than equal to the value in the minor_cn column and greater than 0.
* ```minor_cn```: The copy number of the minor allele. This must be less than equal the value in the major_cn column.
* ```total_cn```: The sum of major_cn and minor_cn.

3. ```sample.purity.txt```: A file storing a scalar purity value between 0 and 1.

A simulated sample input data is under `sample/`. 


## Running CliP with one-step implementation

The caller function `run_clip_main.py` wraps up the CliP pipeline and enables users to implement subclonal reconstruction in one-step. To try CliP with our sample input, you may run:
```
python3 run_clip_main.py sample/sample.snv.txt sample/sample.cna.txt sample/sample.purity.txt
````

A full manual is as follows:

```
usage: run_clip_main.py [-h] [-i SAMPLE_ID] [-e PREPROCESS] [-b] [-f FINAL] [-l LAM] 
                        [-s SUBSAMPLE_SIZE] [-n REP_NUM] [-w WINDOW_SIZE] [-o OVERLAP_SIZE]
                        snv_input cn_input purity_input

positional arguments:
  snv_input             Path/Filename of the snv input.
  cn_input              Path/Filename of the copy number input.
  purity_input          Path/Filename of the purity input.


optional arguments:
  -h, --help            show this help message and exit
  -i SAMPLE_ID, --sample_id SAMPLE_ID
                        Name of the sample being processed. Default is 'sample_id'.
  -l LAM, --lam LAM
                        The penalty parameter (lambda), which usually takes values from 0.01-0.25. 
                        If skipping this parameter, it will return a list of results that take value 
                        of [0.01, 0.03, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25] 
                        by default, and select a preferable one among them.
  -b, --subsampling     Whether doing subsampling or not. Default is not doing the subsampling, 
                        and a flag -b is needed for subsampling.
  -p PREPROCESS, --preprocess PREPROCESS
                        Directory that stores the preprocess results. Default name is 'preprocess_result/'.
  -f FINAL, --final FINAL
                        Directory that stores the final results after postprocessing. Default name
                        is 'final_result/'.
```

The followings parameters are only needed when doing subsampling. We take partitions from 0 to 1 (determined by the `WINDOW_SIZE` parameter, default at `0.05`), then take the VAF of each SNV as an initial estimation of its cellular prevalence (CP) value and randomly assign SNVs that belong to their corresponding windows. The number of sampled SNV is proportional to number of total SNVs belonging to each window. 
```
  -s SUBSAMPLE_SIZE, --subsample_size SUBSAMPLE_SIZE
                        (Required if doing subsampling) The number of SNVs you want to include in 
                        each subsamples. We use 40,000 for 256GB memory. 
  -n REP_NUM, --rep_num REP_NUM
                        (Required if doing subsampling) The number of random subsamples needed.
  -w WINDOW_SIZE, --window_size WINDOW_SIZE
                        Controls the length of the window. Takes value between 0 and 1. Default 
                        is 0.05.
  -o OVERLAP_SIZE, --overlap_size OVERLAP_SIZE
                        Controls the overlapped length of two consecutive windows. Takes value 
                        between 0 and 1. Default is 0.
```


## The CliP Outputs
By default, all outputs will be stored in the folder named `sample_id`, and this name can be changed with the `-i` or `--sample_id ` option.

In final results, cluster index = 0 indicates clonal mutations, while non-zero cluster indexes indicate subclonal mutations.

The final result for CliP is two-fold:
* The subclonal structure, i.e., clustering results: cluster number, the total number of SNVs in each cluster, and the estimated CP for each cluster.
* The mutation assignment, i.e., cluster id for each mutation. This output can then serve as the basis for inference of phylogenetic trees.

## Citation
If you are using this framework, please cite our paper
```
@article{Jiang2021,
    title = {CliP: subclonal architecture reconstruction of cancer cells in DNA sequencing data using a penalized likelihood model},
    author = {Jiang, Yujie and Yu, Kaixian and Ji, Shuangxi and Shin, Seung Jun and Cao, Shaolong and Montierth, Matthew D and Huang, Licai and Kopetz, Scott and Msaouel, Pavlos and Wang, Jennifer Rui and Kimmel, Marek and Zhu, Hongtu and Wang, Wenyi},
    journal = {bioRxiv},
    year = {2021},
    publisher = {Cold Spring Harbor Laboratory}
}
```
