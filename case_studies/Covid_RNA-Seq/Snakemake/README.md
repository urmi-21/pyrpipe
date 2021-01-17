
## Integrating pyrpipe scripts within a workflow management system (Snakemake)
We embedded pyrpipe into the Snakemake workflow management system, using it to download COVID-19 infected human RNA-Seq data.
Data used here is published under accession [SRP287810](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP287810).

We use salmon to quantify the RNA-Seq samples and compile The transcripts TPMs for exploratory analysis with [MetaOmGraph](https://github.com/urmi-21/MetaOmGraph)
The final output file of this workflow, `results_TPM.tsv`, can be directly imported in MetaOmGraph together with the metadata file.



## Running this example

### Build the conda environment (optional)
* Please have conda installed
* Run `conda env create -f env.yaml`
* Activate the environment: `conda activate pyrpipe_human`

### Prepare reference data
To install the reference data used run the following command

`bash prepare_data.sh`

Alternatively, edit the config.yaml to change the reference files

### Execute using snakemake

`snakemake -j <cores>`
 
### To execute on a cluster with slurm

`snakemake -j 25 --cluster-config cluster.json --cluster "sbatch -t {cluster.time} -c 28"`
