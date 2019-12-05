[![Build Status](https://travis-ci.org/urmi-21/pyrpipe.svg?branch=master)](https://travis-ci.org/urmi-21/pyrpipe)
[![Coverage Status](https://coveralls.io/repos/github/urmi-21/pyrpipe/badge.svg?branch=master)](https://coveralls.io/github/urmi-21/pyrpipe?branch=master)

# pyrpipe: python rna-seq pipeliner

### Version
Dev

## Prerequisites
* python 3.6 or higher
* OS: Unix based

Required third-party programs

* [SRA Tools](https://github.com/ncbi/sra-tools) (v. 2.9.6 or higher)
* [Trimgalore](https://github.com/FelixKrueger/TrimGalore)
* [Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml)
* [Stringtie](https://github.com/gpertea/stringtie)
* [Samtools](https://github.com/samtools/samtools)


## Introduction
pyrpipe (Pronounced as "pyre-pipe") is a python package to easily develop RNA-Seq analyses workflow by integrating popular RNA-Seq analysis programs.
pyrpipe can be used on local computers or on HPC environments to manage analysis of RNA-Seq data.



## Installation
Please follow these instruction 
```
git clone https://github.com/urmi-21/pyrpipe.git
pip install -e path/to/dir/pyrpipe
```
