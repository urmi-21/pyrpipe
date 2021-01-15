## Integrating pyrpipe within Snakemake and Nextflow
We embedded pyrpipe into the Snakemake workflow management system, using it to download COVID-19 infected human RNA-Seq data.

Data used here is from [Rother et. al.](https://www.medrxiv.org/content/10.1101/2020.06.08.20122143v1) published under accession [SRP287810](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP287810).

We use salmon to quantify the RNA-Seq samples and compile The transcripts TPMs for exploratory analysis with [MetaOmGraph](https://github.com/urmi-21/MetaOmGraph)
The final output file of this workflow, `results_TPM.tsv`, can be directly imported in MetaOmGraph together with the metadata file.

