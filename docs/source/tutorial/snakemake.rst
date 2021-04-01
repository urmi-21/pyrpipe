Snakemake example
=====================
Since Snakemake direclty supports python, pyrpipe libraries can be driectly imorted into snakemake.
Advantage of using a workflow manager like Snakemake is that it can handle parallel job-scheduling and scale jobs on clusters.


Basic RNA-Seq example
^^^^^^^^^^^^^^^^^^^^^^

A basic example of directly using pyrpipe with Snakemake is provided here.
This example uses yeast RNA-Seq samples from the GEO accession GSE132425.

First run the following bash script to download the yeast reference genome from Ensembl.
The last command also generated a star index with the downloaded genome under the `\refdata\yeast_index` directory.

.. code-block:: bash
    :linenos:

    #!/bin/bash
    mkdir -p refdata
    wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-49/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz -O refdata/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
    wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-49/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.49.gff3.gz -O refdata/Saccharomyces_cerevisiae.R64-1-1.49.gff3.gz
    cd refdata
    gunzip -f *.gz
    #run star index
    mkdir -p yeast_index
    STAR --runThreadN 4 --runMode genomeGenerate --genomeDir yeast_index --genomeFastaFiles Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa --genomeSAindexNbases 10

Next, we create the following `snakefile`.

.. code-block:: python
    :linenos:

    import yaml
    import sys
    import os
    from pyrpipe import sra,qc,mapping,assembly

    ####Read config#####
    configfile: "config.yaml"
    DIR = config['DIR']
    THREADS=config['THREADS']
    ##check required files
    GENOME= config['genome']
    GTF=config['gtf']
    #####Read SRR ids######
    with open ("runids.txt") as f:
        SRR=f.read().splitlines()

    #Create pyrpipe objects
    #parameters defined in ./params will be automatically loaded, threads will be replaced with the supplied value
    tg=qc.Trimgalore(threads=THREADS)
    star=mapping.Star(threads=THREADS)
    st=assembly.Stringtie(threads=THREADS)

    rule all:
        input:
		expand("{wd}/{sample}/Aligned.sortedByCoord.out_star_stringtie.gtf",sample=SRR,wd=DIR),

    rule process:
        output:
		gtf="{wd}/{sample}/Aligned.sortedByCoord.out_star_stringtie.gtf"
	run:
		gtffile=str(output.gtf)
		srrid=gtffile.split("/")[1]
		sra.SRA(srrid,directory=DIR).trim(tg).align(star).assemble(st)

The above snakefile requires some additional files.
A config.yaml file contains path to data and threads information. 

.. code-block:: yaml

    DIR: "results"
    THREADS: 5
    genome: "refdata/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"
    gtf: "refdata/Saccharomyces_cerevisiae.R64-1-1.49.gff3"

Next the snakefile reads a file, `runids.txt` containing SRR accessions.

.. code-block:: yaml

    SRR9257163
    SRR9257164
    SRR9257165

Finally, we need to provide tool parameters for pyrpipe. Create a file `./params/star.yaml` and specify the index in it

.. code-block:: yaml

    --genomeDir: ./refdata/yeast_index/

Now the snakefile could be run using the snakemake command e.g `snakemake -j 8`

pyrpipe_conf.yaml
^^^^^^^^^^^^^^^^^

Users can create a yaml file, `pyrpipe_conf.yaml`, to specify pyrpipe parameters, instead of directly passing them as command line arguments.
When the `pyrpipe_conf.yaml` is found the pyrpipe specific arguments passed via the command-line are ignored.
An example of `pyrpipe_conf.yaml` is shown below with the pyrpipe default values


.. code-block:: yaml

    dry: False          # Only print pyrpipe's commands and not execute anything through pyrpipe_engine module
    threads: None       # Set the number of threads to use
    force: False        # Force execution of commands if their target files already exist
    params_dir: ./params # Directory containing parameters
    logs: True          # Enable or disable pyrpipe logs
    logs_dir: ./pyrpipe_logs    # Directory to save logs
    verbose: False      # Display pyrpipe messages
    memory: None        # Set memory to use (in GB)
    safe: False         # Disable file deletion via pyrpipe commands
    multiqc: False      # Automatically run multiqc after analysis



















    
