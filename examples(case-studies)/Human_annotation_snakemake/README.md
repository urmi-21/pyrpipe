# Introduction

## Integrating pyrpipe scripts within a workflow management system (Snakemake)
We embedded pyrpipe into the Snakemake workflow management system, using it to download human RNA-Seq data with SRAtools, 
quality filter with BBDuk, align reads with Hisat2, assemble transcripts with StringTie and Cufflinks,
and merge the multiple assemblies with Mikado.



## Running this example

### First build the environment
* Please have conda installed
* Run `bash build_env.sh`
* Activate the environment: `conda activate pyrpipe_human`
* Execute Snakemake workflow: `snakemake -j <num threads>`
 
### To execute on a cluster with slurm

`snakemake -j 5 --cluster-config cluster.json --cluster "sbatch -t {cluster.time} -c 28 -C EGRESS`
