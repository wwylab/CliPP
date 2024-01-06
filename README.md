# CliPP
Clonal structure identification through penalizing pairwise differences

## Contents

- [Introduction](#introduction)
- [Prerequisites](#prerequisites)
- [Setting up CliPP](#setting-up-clipp)
- [Structure of CliPP implementation](#structure-of-clipp-implementation)
- [Input data sample](#input-data-sample)
- [Running CliPP with one-step implementation](#running-clipp-with-one-step-implementation)
- [The CliPP Outputs](#the-clipp-outputs)
- [Citation](#citation)

## Introduction
Subpopulations of tumor cells characterized by mutation profiles may confer differential fitness and consequently influence prognosis of cancers. Understanding subclonal architecture has the potential to provide biological insight in tumor evolution and advance precision cancer treatment. Recent methods comprehensively integrate single nucleotide variants (SNVs) and copy number aberrations (CNAs) to reconstruct subclonal architecture using whole-genome or whole-exome sequencing (WGS, WES) data from bulk tumor samples. However, the commonly used Bayesian methods require a large amount of computational resources, a prior knowledge of the number of subclones, and extensive post-processing. Regularized likelihood modeling approach, never explored for subclonal reconstruction, can inherently address these drawbacks. We therefore propose a model-based method, Clonal structure identification through Pair-wise Penalization, or CliPP, for clustering subclonal mutations without prior knowledge or post-processing. The CliPP model is applicable to genomic regions with or without CNAs. CliPP demonstrates high accuracy in subclonal reconstruction through extensive simulation studies. A penalized likelihood framework for subclonal reconstruction will help address intrinsic drawbacks of existing methods and expand the scope of computational analysis for cancer evolution in large cancer genomic studies. Also see our paper: https://www.biorxiv.org/content/10.1101/2021.03.31.437383v1.


## Prerequisites
- MacOS or Linux 
  - If MacOS, CliPP only runs with one core due to a lack of support of OpenMP on this OS.
  - If Linux, CliPP requires OpenMP for paralell computing. For instruction on how to install OpenMP, please check the instruction at http://bioinformatics.mdanderson.org/Software/DeMixT/HowtoinstallOpenMP.docx.
- R [>3.3.1]
- python [>3.5.1]
- NumPy
- SciPy
- pandas

## Setting up CliPP

### Manual Install
```
git clone https://github.com/wwylab/CliPP.git
cd CliPP

python setup.py build
```


### Docker container
We also include a `Dockerfile` in the repository. The user can build and run a Docker contrainer of CliPP to avoid any issues caused by package dependencies. The following provides a tutorial about this. 

- Download and install `Docker` (https://docs.docker.com/get-docker/).
- Download the `Dockerfile` from this repository. Alternatively, you can clone the whole repository to your machine.
- Open a terminal on your machine and change directory to the folder of `Dockerfile`. Create a Docker container for CliPP using the command: 

  ```
  docker build -t clipp .
  ```
  
- Use the sample data in the repository as an example to show how the run the docker container: 

  ```
  git clone https://github.com/wwylab/CliPP.git
  cd CliPP
  docker run -v $(pwd):/Sample clipp python3 /CliPP/run_clipp_main.py -i /Sample/test /Sample/sample/sample.snv.txt /Sample/sample/sample.cna.txt /Sample/sample/sample.purity.txt
  ```
  
**Note**: `-v $(pwd):/Sample` is to mount the the current directory in the host (the machine itself) to the container's directory of `/Sample`, so the input data in the current directory `sample/sample.cna.txt`, `sample/sample.cna.txt` and `sample/sample.purity.txt` can be seen by the container through `/Sample/sample/sample.snv.txt`, `/Sample/sample/sample.cna.txt`, `/Sample/sample/sample.purity.txt`. The output folder is located at `/Sample/test` in the container, which can be accessed at `$(pwd)/test` in the host machine. 

<ins>The installation normally takes a few seconds to 1 minute.</ins>

## Structure of CliPP implementation
The flow chart below shows the CliPP implementation. Raw functions and scripts are under `scr/`. 
![image](https://user-images.githubusercontent.com/14543452/114482762-bf4c1480-9bcc-11eb-8c96-a944611e91d7.png)

Approximately, CliPP can run on samples with up to 50,000 SNVs on a machine with 256GB memory. It may require more memory when there are more SNVs. When that happens, we apply a downsampling strategy, as implemented in all other subclonal reconstruction methods. 

Furthermore, CliPP is automatically run in parallel when users have multiple cores available.

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

## Running CliPP with one-step implementation

The caller function `run_clipp_main.py` wraps up the CliPP pipeline and enables users to implement subclonal reconstruction in one-step. To try CliPP with our sample input, you may run:
```
python run_clipp_main.py sample/sample.snv.txt sample/sample.cna.txt sample/sample.purity.txt
````

A full manual is as follows:

```
usage: run_clipp_main.py [-h] [-i SAMPLE_ID] [-e PREPROCESS] [-b] [-f FINAL] [-l LAM] 
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


## The CliPP Outputs
By default, all outputs will be stored in the folder named `sample_id`, and this name can be changed with the `-i` or `--sample_id ` option.

In final results, cluster index = 0 indicates clonal mutations, while non-zero cluster indexes indicate subclonal mutations.

The final result for CliPP is two-fold:
* The subclonal structure, i.e., clustering results: cluster number, the total number of SNVs in each cluster, and the estimated CP for each cluster.
* The mutation assignment, i.e., cluster id for each mutation. This output can then serve as the basis for inference of phylogenetic trees.

Occasionally, a warning file may appear in the output, which can be caused by two factors:
* For SNVs where the read count information that does not match the CNA-based input data, leading to out-of-bound calculations in the CliPP model, we replace the out-of-bound values with 0.01 and these SNVs are flagged as 1. This is an extremely rare case in real data when the CNA data are good, e.g. TCGA and PCAWG data. When many 1's are observed, the users should look into the quality of the CNA input data.
* While rare, some lambda values may not produce an output. In these instances, CliPP will select the best available result and list the failing lambda values in the warning file.

### Demo output
By default, CliPP processes 11 results using 11 pre-selected lambda values (0.01, 0.03, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25), then recommend the result that follows our pre-defined rule (see **Methods** from the paper). Using the simulation data from directory `sample/` with the command above, CliP will output two files to the default directory `sample_id/final_result/Best_lambda`, with filename `subclonal_structure_lam0.25.txt` and `mutation_assignments_lam0.25`, respectively. `0.25` is the selected lambda value in this case.
The first file outputs an overview of subclonal reconstruction result for this simulation sample, which is shown below. An index is assigned to each cluster, ranging from 0 to n, with clusters sorted in descending order based on their Cellular Prevalence (CP), such that the cluster with the highest CP is assigned an index of 0, the second highest an index of 1, and so on.  
<img width="583" alt="subclonal_structure" src="https://github.com/wwylab/CliPP/assets/14543452/d10294a4-b76a-49ae-91bf-afa8860e164a">

The second file outputs a detailed list of mutation assignment result, i.e., for each mutation, which cluster is it assigned to. The figure below shows the top 10 rows of this mutation assignment table.  
<img width="583" alt="mutation_assignments" src="https://github.com/wwylab/CliPP/assets/14543452/f61a0ba7-429f-4adf-911f-a722fdca8b1c">

<ins>The expected run time for this sample on a personal desktop computer or laptop is around 1 second or less.</ins>

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
