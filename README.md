[![Build Status](https://travis-ci.org/urmi-21/pyrpipe.svg?branch=master)](https://travis-ci.org/urmi-21/pyrpipe)
[![Coverage Status](https://coveralls.io/repos/github/urmi-21/pyrpipe/badge.svg?branch=master)](https://coveralls.io/github/urmi-21/pyrpipe?branch=master)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/pyrpipe)
![PyPI - Downloads](https://img.shields.io/pypi/dd/pyrpipe)
![PyPI - License](https://img.shields.io/pypi/l/pyrpipe)

# pyrpipe: python rna-seq pipeliner

### Version
Dev

## Introduction
pyrpipe (Pronounced as "pyre-pipe") is a python package to easily develop RNA-Seq analyses workflow by integrating popular RNA-Seq analysis programs.
pyrpipe can be used on local computers or on HPC environments to manage analysis of RNA-Seq data.

### What it does
Allows easy development of analyses pipelines in python by providing 
* a high level api to popular bioinformatics tools
* a general api to execute any linux command from python (using the subprocess module)
* logging features to log all the commands and output
* report generating features for easy sharing, reproducing, benchmarking and debugging


## Prerequisites
* python 3.5 or higher
* OS: Unix based

Required third-party programs

* [SRA Tools](https://github.com/ncbi/sra-tools) (v. 2.9.6 or higher)
* [Trimgalore](https://github.com/FelixKrueger/TrimGalore)
* [Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml)
* [Stringtie](https://github.com/gpertea/stringtie)
* [Samtools](https://github.com/samtools/samtools)





## Installation
Please follow these instruction 
```
git clone https://github.com/urmi-21/pyrpipe.git
pip install -r pyrpipe/requirements.txt
pip install -e path_to/pyrpipe
```

## Setting NCBI SRA Toolkit
Use  ```vdb-config -i``` to configure SRA Toolkit. Make sure that:
* Under the **TOOLS** tab, prefetch downloads to is set to public user-repository
* Under the **CACHE** tab, location of public user-repository is not empty

## Setting up the working directory
Please make sure that the working directory (where all data will be downloaded) does not contain any previously downloaded data to avoid problems with ```prefetch```.