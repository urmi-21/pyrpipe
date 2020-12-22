[![Build Status](https://travis-ci.org/urmi-21/pyrpipe.svg?branch=master)](https://travis-ci.org/urmi-21/pyrpipe)
[![Coverage Status](https://coveralls.io/repos/github/urmi-21/pyrpipe/badge.svg?branch=master)](https://coveralls.io/github/urmi-21/pyrpipe?branch=master)
[![Documentation Status](https://readthedocs.org/projects/pyrpipe/badge/?version=latest)](https://pyrpipe.readthedocs.io/en/latest/?badge=latest)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/pyrpipe)
[![install with bioconda](https://anaconda.org/bioconda/plncpro/badges/installer/conda.svg)](https://anaconda.org/bioconda/pyrpipe)
![PyPI - License](https://img.shields.io/pypi/l/pyrpipe)

# pyrpipe: python rna-seq pipeliner



## Introduction
pyrpipe (Pronounced as "pyre-pipe") is a python package to easily develop bioinformatic or any other computational pipeline. pyrpipe provides an easy-to-use framework for importing any UNIX command in python. 
pyrpipe comes with specialized classes and functions to easily code RNA-Seq processing workflows in pure python.
Pipelines in pyrpipe can be created and extended by integrating third-party tools, executable scripts or python libraries in an object oriented manner.
pyrpipe can be used on local computers or on HPC environments to manage analysis of RNA-Seq data.

Users can use the easy-to-use APIs to popular bioinformatic tools provided
with pyrpipe, or can use the methods provided in the pyrpipe_engine module to integrate any other third-party program in their pipeline.

Preprint is available [here](https://www.biorxiv.org/content/10.1101/2020.03.04.925818v3)

Read the docs [here](https://pyrpipe.readthedocs.io/en/latest/?badge=latest)

### Note: pyrpipe version 0.0.5 and above is not compatible with lower versions due to changes in the API design.



## Key Features
* Import any UNIX command in python
* Dry-run feature to check dependencies and commands before execution
* Flexible and robust handling of options and arguments (both Linux and Java style options)
* Auto load command options from .yaml files
* Easily override threads and memory options using global values
* Extensive logging for all the commands
* Automatically verify Integrity of output targets
* Resume feature to restart pipelines/jobs from where interrupted
* Create reports, MultiQC reports for bioinformatic pipelines
* Easily integrated into workflow managers like Snakemake and NextFlow



## What it does
Allows fast and easy development of bioinformatics pipelines in python by providing 
* a high level api to popular RNA-Seq processing tools -- downloading, trimming, alignment, quantificantion and assembly
* optimizes program parameters based on the data
* a general framework to execute any linux command from python
* comprehensive logging features to log all the commands, output and their return status
* report generating features for easy sharing, reproducing, benchmarking and debugging

## What it CAN NOT do by itself
* Schedule jobs
* Scale jobs on HPC/cloud
* Identify parallel steps in pipelines


## Prerequisites
* python 3.6 or higher
* OS: Linux, Mac


## pyrpipe RNA-Seq tools include:

| Tool                                                                                 | Purpose             |
|--------------------------------------------------------------------------------------|---------------------|
| [SRA Tools](https://github.com/ncbi/sra-tools) (v. 2.9.6 or higher)                  | SRA access          |
| [Trimgalore](https://github.com/FelixKrueger/TrimGalore)                             | QC                  |
| [BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/) | QC                  |
| [Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml)                            | Alignment           |
| [STAR](https://github.com/alexdobin/STAR)                                            | Alignment           |
| [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)                     | Alignment           |
| [Kallisto](https://pachterlab.github.io/kallisto/)                                   | Quant               |
| [Salmon](https://combine-lab.github.io/salmon/)                                      | Quant               |
| [Stringtie](https://github.com/gpertea/stringtie)                                    | Transcript Assembly |
| [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/)                           | Transcript Assembly |
| [Samtools](https://github.com/samtools/samtools)                                     | Tools               |



### Read the documentation [here](https://pyrpipe.readthedocs.io/en/latest/?badge=latest)

## Installation
Please follow these instruction: 

### To create a new Conda environment (optional):
1. Download and install Conda
2. `conda create -n pyrpipe python=3.7`
3. `conda activate pyrpipe`

**NOTE: You need to install the third-party tools to work with pyrpipe. We recomend installing these through [bioconda](https://bioconda.github.io/) where possible. 
A recommended environment file, using python 3.7, is given [here](https://github.com/urmi-21/pyrpipe/blob/master/tests/test_environment.yml).
It is best to [share your conda environment files](https://stackoverflow.com/questions/41274007/anaconda-export-environment-file) with pyrpipe scripts to ensure reproducibility.**

### Install latest stable version

#### Through conda

```
conda install -c bioconda pyrpipe 
```
 
#### Through PIP

```
pip install pyrpipe --upgrade
```

If above command fails due to dependency issues, try: 
1. Download the [requirements.txt](https://github.com/urmi-21/pyrpipe/blob/master/requirements.txt)
2. `pip install -r requirements.txt`
3. `pip install pyrpipe`

To run tests:
1. Download the [test set](https://github.com/urmi-21/pyrpipe/tree/master/tests) ([direct link](https://minhaskamal.github.io/DownGit/#/home?url=https://github.com/urmi-21/pyrpipe/tree/master/tests))
2. `pip install pytest`
3. To build test_environment. Please [READ THIS](https://github.com/urmi-21/pyrpipe/blob/master/tests/README.md)
4. From pyrpipe root directory, run `py.test tests/test_*`


### Install current developer version
```
git clone https://github.com/urmi-21/pyrpipe.git
pip install -r pyrpipe/requirements.txt
pip install -e path_to/pyrpipe

#Running tests; From pyrpipe root perform
#To build test_environment (This will download tools): 
cd tests ; . ./build_test_env.sh
#in same terminal
py.test tests/test_*
```

## Setting NCBI SRA Toolkit
Use  ```vdb-config -i``` to configure SRA Toolkit. Make sure that:
* Under the **TOOLS** tab, prefetch downloads to is set to public user-repository
* Under the **CACHE** tab, location of public user-repository is not empty

## Funding

This work is funded in part by the National Science Foundation award IOS 1546858, "Orphan Genes: An Untapped Genetic Reservoir of Novel Traits".



