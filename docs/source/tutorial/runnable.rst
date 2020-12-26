Importing tools into python
-----------------------------

pyrpipe's `Runnable` can be used to import any Unix command into python.
The `Runnable` class implements the `run()` method which checks required dependencies, monitors execution, execute commands via pyrpipe_engine, and verify target files.
We will first download the *E. coli* genome file to use as input to orfipy.



.. code-block:: bash 
    :linenos:

    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
    gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz

Basic example
==============
In this tutorial we will consider the tool `orfipy`, a tool for fast and flexible ORFs, and import it in python.



.. code-block:: python
    :linenos:

    from pyrpipe.runnable import Runnable
    infile='GCF_000005845.2_ASM584v2_genomic.fna'
    #create a Runnable object
    orfipy=Runnable(command='orfipy')
    #specify orfipy options; these can be specified into orfipy.yaml too
    param={'--outdir':'orfipy_out','--procs':'3','--dna':'orfs.fa'}
    orfipy.run(infile,**param)

Save the above script in `example.py` and execute it using `python example.py --dry-run`.
pyrpipe should generate the orfipy command and display on screen during dry-run. 
Running it without `--dry-run` flag will generate the `orfipy_out/orfs.fa` file.

Requirements and targets
========================
One can specify required dependencies and expected target files in the run() method
Replacing the call to `run()` with the following will verify the required files and the target files.
If command is interrupted, pyrpipe will scan for `Locked` taget files and resume from where the pipeline was interrupted.

.. code-block:: python
    :linenos:
    
    orfipy.run(infile,requires=infile,target='orfipy_out/orfs.fa',**param)


Building APIs
==============
Users can extend the Runnable class and specify classes dedicated to tools.
Extra functionalities can be added by defining custom behaviour of classes.
An example is the RNA-Seq API built using pyrpipe framework.


















