## Introduction

This basic example demonstrates importing pyrpipe libraries into Snakemake.
The example downloads data from NCBI-SRA, trims reads using Trim Galore, aligns reads to genome using STAR and assembles transcripts with stringtie

## Download reference data

Run the `prepare_data.sh` script to download the reference data.
Alternatively, edit the config.yaml to specify location of reference genome.

## Executing

Run the example using `snakemake` command
