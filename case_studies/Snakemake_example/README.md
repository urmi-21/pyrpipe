## Introduction

This basic example demonstrate importing pyrpipe libraries into snakemake.
This example downloads data from NCBI-SRA, trims read using Trim Galore, align reads to genome using STAR and assemble transcripts with stringtie

## Download reference data

Run the `prepare_data.sh` script to download the reference data.
Alternatively, edit the config.yaml to specify location of reference genome.

## Executing

Run the example using `snakemake` command
