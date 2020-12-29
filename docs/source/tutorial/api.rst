Building APIs
==============
Users can extend the Runnable class and specify classes dedicated to tools.
Extra functionalities can be added by defining custom behaviour of classes.
The RNA-Seq API built using pyrpipe the framework.

A small example is presented here to build a class for orfipy tool

.. code-block:: python
    :linenos:

    from pyrpipe import Runnable
    
    from pyrpipe.runnable import Runnable
    from pyrpipe import pyrpipe_utils as pu
    from pyrpipe import sra
    from pyrpipe import _threads,_mem,_force
    import os

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
            super().__init__(*args,command='orfipy',**kwargs)
            self._deps=[self._command]
            self._param_yaml='orfipy.yaml'  
            #valid arguments for orfipy
            self._valid_args=['--min','--between-stops','--include-stop','--dna','--pep','--bed','--bed12','--procs','--chunk-size','--outdir'] #valid arguments for orfipy
    
            #resolve threads to use
            """
            orfipy parameter for threads is --procs
            if threads is passed in __init__() it will be used
            else if --procs is found in orfipy.yaml that will be used
            else if --procs is found in the passed **kwargs in __init__() it will be used
            else the default value i.e. _threads will be used
            if default value is None nothing will be done
            after the function, --procs and its value will be stored in self._kwargs, and _threads variable will be stored in the Orfipy object.
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
    
            #In this example use orfipy on only first fastq file
            internal_args=(sra_object.fastq_path,)
            internal_kwargs={"--bed":sra_object.srr_accession+"_ORFs.bed","--outdir":out_dir}
    
            #call run
            status=self.run(*internal_args,objectid=sra_object.srr_accession,target=out_file,**internal_kwargs)
    
            if status:
                return out_file
    
            return ""


The above class, Orfipy, we created can be used directly with SRA type objects via the find_orfs function. 
We can also still use the `run()` method and provide any input to orfipy.

.. code-block:: python
    :linenos:

    #create object
    orfipy=Orfipy()
    #use run()
    orfipy.run('test.fa',**{'--dna':'d.fa','--outdir':'of_out'},requires='test.fa',target='of_out/d.fa')

    #use the api function to work with SRA
    srr=sra.SRA('SRR9257212')
    orfipy.find_orfs(srr)


Now try passing orfipy parameters from orfipy.yaml file. Create `params/orfipy.yaml` and add the following options into it.


.. code-block:: yaml
    :linenos:
    
    --min: 36
    --between-stops: True
    --include-stop: True


Now re-run the python code and it will automatically read orfipy options from `./params/orfipy.yaml`.

