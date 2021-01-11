# GTEx RNA-Seq processing pipeline
This pipeline implemented in pyrpipe processes GTEx RNA-Seq samples. 
The pipeline automatically downloads the GTEx raw data and submits multiple jobs to batch process the data.
This pipeline assumes users are working with slurm job scheduler. 
The function `submit_slurm_job()` in `driver_GTEx.py` can be modified to change how the jobs are submitted.

The pipeline can run in two modes
1. Alignment to genome followed by transcript assembly
2. Pseudo-alignment and quantification


## Requirements
* pyrpipe version 0.0.5 or higher
* Tools required: Biobambam2, Sambamba, STAR, Stringtie, Salmon
* Gen3-client to download GTEx data. **Users must have dbGap authorization to access the data**

## Steps
1. Download and setup [gen3-client](https://gen3.org/resources/user/gen3-client/)
2. Download the GTEx manifest file and save it as `manifest.json`. See [this for help](https://anvilproject.org/learn/reference/gtex-v8-free-egress-instructions).
3. Edit the `config.yaml` file as required.
4. Edit the yamls files under params directory to specify or modify tool parameters.
5. Run the pipeline as `python driver_GTEx.py <GTEx id file> <output directory> <slurm_job_name>`.
where `GTEx id file` is a file containing GTEx sample accessions. See sample_ids file for an example.