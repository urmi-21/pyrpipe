# Introduction

## *Arabidopsis thaliana* transcript assembly using pyrpipe checkpoints
We downloaded raw read RNA-Seq data for Arabidopsis from SRA, performed quality control using BBDuk, aligned reads to the genome using Hisat2 and
assembled transcripts using StringTie.

## To run the example:

1. Build the conda environment:

    `conda env create -f environment.yml`
    
2. Switch to the environment:

    `conda activate pyrpipe_examples`

3. Execute the code in `.ipynb` or `.py` file
