## Integrating pyrpipe within Snakemake and Nextflow
We embedded pyrpipe into the Snakemake workflow management system, using it to download COVID-19 infected human RNA-Seq data.

Data used here is from accession [SRP287810](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP287810).

We use salmon to quantify the RNA-Seq samples and compile The transcripts TPMs for exploratory analysis with [MetaOmGraph](https://github.com/urmi-21/MetaOmGraph)
The final output file of this workflow, `results_TPM.tsv`, can be directly imported in MetaOmGraph together with the metadata file.

To visualize and analyse data with MetOmGraph download *Covid_RNA-Seq/MOGproject-monocytes-60241genes-29samples-HCQtreat-2021-1-17.zip* file and import the .mog file into MetaOmGraph. [Download MetaOmGraph here](https://github.com/urmi-21/MetaOmGraph)
