
## Integrating pyrpipe scripts within a workflow management system (Nextflow)
This example shows usage of pyrpipe with Nextflow workflow management system. The example is identical to the [Snakemake example](https://github.com/urmi-21/pyrpipe/tree/master/case_studies/Covid_RNA-Seq).
This example downloads COVID-19 infected human RNA-Seq data from [Rother et. al.](https://www.medrxiv.org/content/10.1101/2020.06.08.20122143v1) published under accession [SRP287810](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP287810).

Uses salmon to quantify the RNA-Seq samples and compile The transcripts TPMs for exploratory analysis with [MetaOmGraph](https://github.com/urmi-21/MetaOmGraph)
The final output file of this workflow, `results_TPM.tsv`, can be directly imported in MetaOmGraph together with the metadata file.



## Running this example

### Build the conda environment (optional)
* Please have conda installed
* Run `conda env create -f env.yaml`
* Activate the environment: `conda activate pyrpipe_human`

### Prepare reference data
To install the reference data and build Salmon index, run [this file]() as 

`bash prepare_data.sh`

### Execute using Nextflow

`nextflow run run_quant.nf`