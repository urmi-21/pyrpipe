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
The RNA-Seq API built using pyrpipe the framework.

A small example is presented here to build a class for orfipy tool

.. code-block:: python
    :linenos:
    
    class Orfipy(Runnable):
    """
    Extends Runnable class
    Attributes
    ----------
    
    """ 
    def __init__(self,*args,threads=None,mem=None,**kwargs):
        """
        init an Orfipy object

        Parameters
        ----------
        *args : tuple
            Positional arguements to orfipy
        threads : int, optional
            Threads to use for orfipy. This will override the global --threads parameter supplied to pyrpipe. The default is None.
        mem : int, optional
            Maximum memory to use in MB. The default is None.
        **kwargs : dict
            options for orfipy
        Returns
        -------
        None.

        """
        super().__init__(*args,**kwargs)
        self._command='orfipy'
        self._deps=[self._command]
        #self._param_yaml='orfipy.yaml' this is initialized by default; can override the path to yaml if required
        self._valid_args=['--dna','--pep','--bed','--bed12','--procs','--chunk-size'] #valid arguments for orfipy
        
        #resolve threads to use
        """
        orfipy parameter for threads is --procs
        if threads is passed in __init__() it will be used
        else if --procs is found in orfipy.yaml that will be used
        else if --procs is found in the passed **kwargs in __init__() it will be used
        else the default value i.e. _threads will be used
        if default value is None nothing will be done
        _threads variable will be stored in the Orfipy object. Use self._threads in all functions to access the threads parameters
        """
        self.resolve_parameter("--procs",threads,_threads,'_threads')
        #resolve memory to use
        """
        default value is None--> if mem is not supplied don't make the self._mem variable
        """
        self.resolve_parameter("--chunk-size",mem,None,'_mem')
        
        ##now we write a custom function that can be used with an SRA object
        
        def find_orfs(self,sra_object):
            out_dir=sra_object.directory
            
            out_file=os.path.join(out_dir,sra_object.srr_accession+"_ORFs.bed")
            
            if not _force and pu.check_files_exist(out_file):
                pu.print_green('Target files {} already exist.'.format(out_file))
                return out_file
            
            #find layout and fq file paths
            internal_args=(sra_object.fastq_path,)
            internal_kwargs={"--bed":sra_object.srr_accession+"_ORFs.bed","--outdir":out_dir}
            
            #call run
            status=self.run(internal_args,objectid=sra_object.srr_accession,target=out_file,**internal_kwargs)
            
            if status:
                return out_file
            
            return ""
















