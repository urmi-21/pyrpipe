[![Build Status](https://travis-ci.org/urmi-21/pyrpipe.svg?branch=master)](https://travis-ci.org/urmi-21/pyrpipe)
[![Coverage Status](https://coveralls.io/repos/github/urmi-21/pyrpipe/badge.svg?branch=master)](https://coveralls.io/github/urmi-21/pyrpipe?branch=master)
[![Documentation Status](https://readthedocs.org/projects/pyrpipe/badge/?version=latest)](https://pyrpipe.readthedocs.io/en/latest/?badge=latest)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/pyrpipe)
![PyPI - Downloads](https://img.shields.io/pypi/dm/pyrpipe)
![PyPI - License](https://img.shields.io/pypi/l/pyrpipe)

# pyrpipe: python rna-seq pipeliner



## Introduction
pyrpipe (Pronounced as "pyre-pipe") is a python package to easily develop RNA-Seq analyses workflow by integrating popular RNA-Seq analysis programs in an object oriented manner.
pyrpipe can be used on local computers or on HPC environments to manage analysis of RNA-Seq data. Users can use the easy-to-use APIs to popular bioinformatic tools provided
with pyrpipe, or can use the methods provided in the pyrpipe_engine module to integrate any other third-party program in their pipeline.

Preprint is available [here](https://www.biorxiv.org/content/10.1101/2020.03.04.925818v1)

### Note: pyrpipe version 0.0.4 and above is not compatible with lower versions due to changes in the API design.


### What it does
Allows fast and easy development of bioinformatics pipelines in python by providing 
* a high level api to popular bioinformatics tools
* optimizes program parameters based on the data
* a general api to execute any linux command from python (using the subprocess module)
* comprehensive logging features to log all the commands, output and their return status
* report generating features for easy sharing, reproducing, benchmarking and debugging


## Prerequisites
* python 3.6 or higher
* OS: Linux, Mac


## pyrpipe provides API to:

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
| [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)                       | Transcript Assembly |
| [Samtools](https://github.com/samtools/samtools)                                     | Tools               |
| [Portcullis](https://github.com/maplesond/portcullis)                                | Tools               |
| [Mikado](https://github.com/EI-CoreBioinformatics/mikado)                            | Tools               |




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

