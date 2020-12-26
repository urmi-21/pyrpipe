
Basic RNA-Seq processing
-------------------------

After setting up the environment, one can import pyrpipe modules in python and start using it.
In this example we will use the *A. thaliana* Paired-end RNA-Seq run *SRR976159*.

Required files
================

We need to download the reference genome and annotation for *A. thalinia*. This can be done inside python script too.
For simplicity we just download these using the `wget` command from the terminal.

.. code-block:: bash
    :linenos:

    wget ftp://ftp.ensemblgenomes.org/pub/release-46/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
    gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
    wget ftp://ftp.ensemblgenomes.org/pub/release-46/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.46.gtf.gz
    gunzip Arabidopsis_thaliana.TAIR10.46.gtf.gz



Simple pipeline
================

RNA-Seq processing is as easy as creating required objects and executing required functions.
The following python script provides a basic example of using pyrpipe on publicly available RNA-Seq data.


.. code-block:: python
    :linenos:
    
    from pyrpipe import sra,qc,mapping,assembly
    #define some vaiables
    run_id='SRR976159'
    working_dir='example_output'
    gen='Arabidopsis_thaliana.TAIR10.dna.toplevel.fa'
    ann='Arabidopsis_thaliana.TAIR10.46.gtf'
    star_index='star_index/athaliana'
    #initialize objects
    #creates a star object to use with threads
    star=mapping.Star(index=star_index,genome=gen,threads=4)
    #use trim_galore for trimming
    trim_galore=qc.Trimgalore()
    #Stringtie for assembly
    stringtie=assembly.Stringtie(guide=ann)
    #create SRA object; this will download fastq if doesnt exist
    srr_object=sra.SRA(run_id,directory=working_dir)
    #create a pipeline using the objects
    srr_object.trim(trim_galore).align(star).assemble(stringtie)
    
    #The assembled transcripts are in srr_object.gtf
    print('Final result',srr_object.gtf)
    
The above code defines a simple pipeline (in a single line: Line 18) that: 

- Downloads fastq files from NCBI-SRA 
- Uses Trim Galore for trimming
- Uses Star for alignemnt to ref genome
- Uses Stringtie for assembly


**A line by line explanation:**

1. Imports the required pyrpipe modules

3. Lines 3 to 7 defines the variables for reference files, output directory, and star index

10. Creates a Star object. It takes index and genome as parameters. It will automatically verify the index and if an index is not found, it will use the genome to build one and save it to the index path provided.

12. Creates a Trimgalore object

14. Creates a Stringtie object

16. Creates an SRA object. This represents the RNA-Seq data. If the raw data is not available on disk it auto-downloads it via fasterq-dump.

18. This is the pipeline which describes a series of operations. The SRA class implements the *trim()*, *align()* and *assemble()* methods.
  -  *trim()* takes a qc type object, performs trimming via *qc.perform_qc* method and the trimmed fastq are updated in the SRA object
  -  *align()* takes a mappping type object, performs alignemnt via *mapping.perform_alignemnt* method. The resulting bam file is stored in *SRA.bam_path*
  -  *assemble()* takes a assembly type object, performs assembly via *mapping.perform_assembly* method. The resulting gtf file is stored in *SRA.gtf*


Executing the pipeline
=======================
To execute the pipeline defined above, save the python code in a file `pipeline.py`. 
The code can be executed as any python script using the python command:

.. code-block:: bash
    :linenos:

    python pipeline.py

Or it can be executed using the pyrpipe command and specifying the input script with the --in option

.. code-block:: bash
     :linenos:
 
     pyrpipe --in pipeline.py


One can specify pyrpipe specific options too

.. code-block:: bash
     :linenos:

     python pipeline.py --threads 10 --dry-run
        #OR
     pyrpipe --in pipeline.py --threads 10 --dry-run

The above two commands are equivalent and specifies pyrpipe to use 10 threads. Thus 10 threads will be used for each of the tool except for STAR where we explicitly
specied to use 4 threads during object creation.

The other option provided here is the `--dry-run` option and this option *turns off* the pyrpipe_engine and any command passed to the pyrpipe_engine is not actually executed but just displayed on screen and logged. During dry run the Runnable class also verifies the file dependencies (if any). More details are provided in the later chapters of the tutorial.

We recommend using the dry-run option before actually starting a job to make sure all parameters/dependencies are correct.


Specifying tool parameters
===========================

pyrpipe supports auto-loading of tool parameters specified in .yaml files. The .yaml files must be stored is a directory and can be specified using the `--param-dir` option. The default value is `./params`. The files must be named as `<tool>.yaml`, for example star.yaml for STAR parameters.
These parameters are loaded during object creation and user can easily override these during execution.

Create a directory `params` in the current directory and make a file `star.yaml` inside `params`. Add the following to `star.yaml` and rerun `pipeline.py` using the dry run option.

.. code-block:: yaml
    :linenos:

    --outSAMtype: BAM Unsorted SortedByCoordinate
    --outSAMunmapped: Within 
    --genomeLoad: NoSharedMemory 
    --chimSegmentMin: 15 
    --outSAMattributes: NH HI AS nM NM ch
    --outSAMattrRGline: ID:rg1 SM:sm1

If you did everything correctly, you wil notice that now the STAR commands contain these specified parameters.


Updating parameters dynamically
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The parameters specified in the yaml will be replaced by any parameters provided during object creation.
For example, consider the star.yaml specifying --runThreadN as 20

.. code-block:: yaml
    :linenos:
    
    --outSAMtype: BAM Unsorted SortedByCoordinate
    --outSAMunmapped: Within 
    --genomeLoad: NoSharedMemory 
    --runThreadN: 20

Now, consider creating a star object in the following scenarios

.. code-block:: python
    :linenos:

    star1=Star(index='index') #will use 20 threads as mentioned in the yaml file
    star2=Star(index='index',threads=5) #will use 5 threads
    star3=Star(index='index',**{'--runThreadN':'10'}) #will use 10 threads
    star4=Star(index='index') #initialized with 20 threads
    star4.run(...,**{'--runThreadN':'10'}) #will use 10 threads for this particular run; '--runThreadN':'20' remains in the star4 object


















