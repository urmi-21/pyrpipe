[![Build Status](https://travis-ci.org/urmi-21/pyrpipe.svg?branch=master)](https://travis-ci.org/urmi-21/pyrpipe)
[![Coverage Status](https://coveralls.io/repos/github/urmi-21/pyrpipe/badge.svg?branch=master)](https://coveralls.io/github/urmi-21/pyrpipe?branch=master)
[![Documentation Status](https://readthedocs.org/projects/pyrpipe/badge/?version=latest)](https://pyrpipe.readthedocs.io/en/latest/?badge=latest)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/pyrpipe)
![PyPI - Downloads](https://img.shields.io/pypi/dd/pyrpipe)
![PyPI - License](https://img.shields.io/pypi/l/pyrpipe)

# pyrpipe: python rna-seq pipeliner


## Introduction
pyrpipe (Pronounced as "pyre-pipe") is a python package to easily develop RNA-Seq analyses workflow by integrating popular RNA-Seq analysis programs in an object oriented manner.
pyrpipe can be used on local computers or on HPC environments to manage analysis of RNA-Seq data. Users can use the easy-to-use APIs to popular bioinformatic tools provided
with pyrpipe, or can use the methods provided in the pyrpipe_engine module to integrate any other third-party program in their pipeline.

### What it does
Allows fast and easy development of bioinformatics pipelines in python by providing 
* a high level api to popular bioinformatics tools
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






## Installation
Please follow these instruction 

### Install latest stable version

```
pip install pyrpipe
```

### Install current developer version
```
git clone https://github.com/urmi-21/pyrpipe.git
pip install -r pyrpipe/requirements.txt
pip install -e path_to/pyrpipe

#Running tests; From pyrpipe root perform

py.test tests/
```

## Setting NCBI SRA Toolkit
Use  ```vdb-config -i``` to configure SRA Toolkit. Make sure that:
* Under the **TOOLS** tab, prefetch downloads to is set to public user-repository
* Under the **CACHE** tab, location of public user-repository is not empty

