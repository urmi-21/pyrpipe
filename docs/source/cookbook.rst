======================
Cookbook
======================

.. contents::

Using SRA objects
----------------------



.. code-block:: python
    :linenos:

    from pyrpipe.sra import SRA #imports the SRA class
    
    #create an SRA object using a valid run accession
    """
    this checks if fastq files already exist in the directory,
    otherwise downloads the fastq files and stores the path in the object
    """
    myob=SRA('SRR1168424',directory='./results')
    
    #create an SRA object using fastq paths
    myob=SRA('SRR1168424',fastq='./fastq/SRR1168424_1.fastq',fastq2='./fastq/SRR1168424_2.fastq')
    
    #create an SRA object using sra path
    myob=SRA('SRR1168424',sra='./sra/SRR1168424.sra')
    
    #accessing fastq files
    print(myob.fastq,myob.fastq2)
    
    #check if fastq files are present
    print (myob.fastq_exists())
    
    #check sra file
    print (myob.sra_exists())
    
    #delete fastq
    myob.delete_fastq()
    
    #delete sra
    myob.delete_sra()
    
    #download fastq 
    myob.download_fastq()
    
    #trim fastq files
    myob.trim(qcobject)
    

Using RNASeqQC objects
-----------------------

`RNASeqQC` objects can be used for quality control and trimming. These are defined in the `qc` module.
Following example uses Trimgalore class but is applicable to any class extending the `RNASeqQC` class

.. code-block:: python
    :linenos:

    from pyrpipe.qc import Trimgalore
    
    #create trimgalore object
    tgalore=Trimgalore()
    #print category
    print(tgalore._category) #should print Aligner
    
    #use with SRA object
    """
    Following will trim fastq and update fastq paths in the sraobject
    """
    sraobject.trim(tgalore)
    #following trims and returns qc
    fq1,fq2=tgalore.perform_qc(sraobject)
    
    #run trimgalore using user arguments; provide any arguments that Runnable.run() can take
    tgalore.run(*args,**kwargs)


    
Using Aligner objects
----------------------

`Aligner` objects from the `mapping` module can be used to perform alignment tasks.

.. code-block:: python
    :linenos:

    from pyrpipe.mapping import Star
    
    #create a star object
    star=Star(index='path_to_star_index')
    
    #print category
    print(star._category) #should print Aligner
        
    #perform alignment using SRA object
    bam=star.perform_alignment(sraobject)
    #or
    sraobject.align(star)
    bam=sraobject.bam_path
    
    #execute STAR with any arguments and parameters
    kwargs={'--outFilterType' : 'BySJout',
            '--runThreadN': '6',
            '--outSAMtype': 'BAM SortedByCoordinate',
            '--readFilesIn': 'SRR3098744_1.fastq SRR3098744_2.fastq'
            }
    star.run(**kwargs)
    
    

Using Assembler objects
-----------------------
    
`Assembler` objects are defined the the assembly module and can be used for transcript assembly.


.. code-block:: python
    :linenos:

    from pyrpipe.assembly import Stringtie
    
    #create a stringtie object
    stringtie=Stringtie(guide='path_to_ref_gtf')
    
    #perform assembly using SRA object
    """
    Note: following first runs star to perform alignment. After alignment the Sorted
    BAM file is stored in the sraobject.bam_path attribute and returns the modified sraobject
    The assemble function requires a valid bam_path attribute to work.
    """
    sraobject.align(star).assemble(stringtie)
    #Or manually set bam_path
    sraobject.bam_path='/path/to/sorted.bam'
    sraobject.assemble(stringtie)
    
    #use perform_assembly function
    result_gtf=stringtie.perform_assembly('/path/to/sorted.bam')
    
    #run stringtie with user arguments
    stringtie.run(verbose=True, **kwargs)
    
    

Using Quantification objects
-----------------------------

`Quantification` type objects can perform quantification and are defined inside the `quant` module.



.. code-block:: python
    :linenos:

    from pyrpipe.assembly import Salmon
    
    #create salmon object
    """
    A valid salmon idex is required. If index is not found it is built using the provided transcriptome
    """
    salmon=Salmon(index='path/to/index',transcriptome='path/to/tr')
    
    #directly quantify using SRA object
    sraobject.quant(salmon)
    #or trim reads before quant
    sraobject.trim(tgalore).quant(salmon)
    print('Result file',sraobject.abundance)
    
    #use perform quant function
    abundance_file=salmon.perform_quant(sraobject)
    
    #use salmon with user defined arguments
    salmon.run(**kwargs)
    
    
Using RNASeqTools objects
-----------------------------

The `RNASeqTools` type is defined in `tools` module. This contains various tools used routinely for RNA-Seq data processing/analysis


.. code-block:: python
    :linenos:

    from pyrpipe.tools import Samtools
    
    #create samtools object
    samtools=Samtools(threads=6)
    
    #convert sam to sorted bam
    bam=samtools.sam_sorted_bam('sam_file')
    
    #merge bam files
    mergedbam=samtools.merge_bam(bamfiles_list)
        
    
    #run samtools with used defined arguments
    """
    NOTE: the Runnable.run() method accepts a subcommand argument that allows user to procide a subcommand like samtools index or samtools merge
    """
    samtools.run(*args,subcommand='index',**kwargs)
    
    
Using Runnable objects
-----------------------

The `Runnable` class, defined inside the `runnable` module, is the main parent class in `pyrpipe` i.e. all other classes borrows its functionality.
User can directly create `Runnable` objects to define their own tools/commands.

A full example to build APIs is here: :ref:`API Overview`

.. code-block:: python
    :linenos:

    from pyrpipe.runnable import Runnable
    
    #say you want to use the tool orfipy
    orfipy=Runnable(command='orfipy')
    
    #execute orfipy as
    orfipy.run(*args,**kwargs)
    
    #another example using Unix grep 
    grep=Runnable(command='grep')
    grep.run('query1','file1.txt',verbose=True)
    grep.run('query2','file2.txt',verbose=True)
    
    #extend Runnable to build more complex APIs that fit with each other
    """
    One can create classes extending the Runnable class.
    Full example is given in the tutorial
    """
    class Orfipy(Runnable):
        def __init__(self,*args,threads=None,mem=None,**kwargs):
            super().__init__(*args,command='orfipy',**kwargs)
            self._deps=[self._command]
            self._param_yaml='orfipy.yaml'
        
        #create special API functions that can work with other objects
        def find_orfs(self,sra_object):
            #define logic here and gather command options and parameters 
            
            #call the self.run() function and check values
            
            #return a useful value


Using pyrpipe_engine module
----------------------------

The `pyrpipe_engine` module contain functions that creates new processes and enable executing commands.
User can directly import `pyrpipe_engine` module and start using the functions.
This is very useful for quickly executing commands without having to create a Runnable object.
A table describing functions implements in `pyrpipe_engine` is provided in the tutorial :ref:`Engine Overview`


.. code-block:: python
    :linenos:
    
    
    import pyrpipe_engine as pe
    
    #execute_command: Runs a command, logs the status and returns the status (True or False)
    pe.execute_command(['ls', '-l'],logs=False,verbose=True)
    
    #get_shell_output Runs a command and returns a tuple (returncode, stdout and stderr)
    """
    NOTE: only this function supports shell=True
    """
    result=pe.get_shell_output(['head','sample_file'])
    #result contains return code, stdout, stderr
    print(result)
    
    #make a function dry-run compatible
    """
    when --dry-run flag is used this function will be skipped and the first parameter 'cmd' will be printed to screen
    """
    @pe.dryable
    def func(cmd,...)
        #function logic here
        
    """
    Other way to work with dry-run flag is to directly import _dryrun flag
    """
    
    from from pyrpipe import _threads,_force,_dryrun
    def myfunction(...):
        if _dryrun: 
            print('This is a dry run')
            return
        
        #real code here...
        
    
    

Using pyrpipe_utils module
----------------------------

The `pyrpipe_utils` module contains several helpful functions that one needs to frequently access for typical computational pipelines.
A table describing functions implements in `pyrpipe_utils` is provided in the tutorial :ref:`Utils Overview`


.. code-block:: python
    :linenos:
    
    
    import pyrpipe_utils as pu
    
    #check if files exist
    pu.check_files_exist('path/f1','path/f2','path/f3') #returns bool
    
    #get filename without extention
    pu.get_file_basename('path/to/file.ext')
    
    #get file directory
    pu.get_file_directory('path/to/file.ext')
    
    #create a directory
    pu.mkdir('path/to/dir')
    
    #Search .txt files in a directly
    pu.find_files('path/to/directory','.*\.txt$',recursive=False,verbose=False)















    
    
    
