======================================================
Basic usage
======================================================

***********
Philosophy
***********
pyrpipe's simple and flexible framework allows any UNIX command to be imported in python. This can be used to build computational pipelines in pure python.
These tools then can be used and re-used flexibly in python as python objects.
pyrpipe follows an object oriented approach in which the tools, its parameters and data are encapsulated as classes.
This is made possible via the single `Runnable` class and the pyrpipe_engine module. 
Calling `run()` method of a `Runnable` object executes the commands via pyrpipe_engine.
A number of routinely used helper functions are provided in pyrpipe_utils module.


********************
API to RNA-Seq tools
********************

High-level APIs to bioinformatic RNA-Seq processing tools comes with pyrpipe for easy development of RNA-Seq processing pipelines.
Objects of these classes can be the used to build workflows. A quick explanation is below:

- The central class for RNA-Seq analysis is the `SRA` class defined inside the sra module. This represents raw RNA-Seq fastq data and can automatically download data from NCBI-SRA. 
- The qc module contains classes to represent quality-control tools: Trim galore and bbduk.sh from BBMAP.
- The mapping module contains classes for alignment tools: STAR, Hisat2 and Bowtie2
- The asembly module contains classes for transcript assembly: Stringtie and Cufflinks
- The quant module contains classes for quantification via pseudo-alignment: Salmon and Kallisto
- We have implemented methods such the basic RNA-Seq processing could be performed by applying series of operations


Importing Unix commands/tools in python
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each tool is represented as a class that provides access to various functionality via its methods. The methods implented are:

* Constructor (init()): The init() method is a python reserved methods to construct objects.
In pyrpipe each constructor can take some general, optional arguments like *threads* and *max_memory*.
For example a *Samtools* object can be constructed using the constructer::
        obj=Samtools(threads=8)
These arguments are stored with the object and are used by default.
User has the option to override them later when calling other methods, for example::
        obj.sort_bam(...threads=4...)
If the optional arguments are omitted from the constructor, pyrpipe initializes them to best available values, for example *threads* will be initialized to *os.cpu_count()*.

* API functions: These functions provide specific utility of a given tool e.g. the *sort_bam()* sorts bam file using *samtools sort* command.
These functions can have multiple arguments such as the required input/output files and other optional arguments.
These arguments are then parsed by `pyrpipe` and appropriate shell command is constructed to invoke the tool.
These functions also take arbitary dict arguments (kwargs) that can be used to specify arguments other than provided by the function definition.
If an argument is provided using both the function definition and kwargs, the latter will override the former::
        obj=Samtools(threads=8) #init object
        #examples of passing arguments
        obj.sort_bam(...) #will use 8 threads as specified in the constructor
        obj.sort_bam(...threads=4...) #will use 4 threads, override constructor
        obj.sort_bam(...threads=4...**{"-@":"6"}) #will use 6 threads, override the threads argument from function definition


* Run function: Each class representing a tool has a `run()` method. This method is responsible for actual execution of the tool i.e. all API functions call this function to execute the required shell command.
The run function takes a number of argument which specify the execution related behaviour. Most importantly it takes an arbitary dict argument (kwargs) which actually contains the arguments to be passed on to the tool. User can also directly use the run function to execute a tool.
For example::

        obj=Samtools(threads=8) #init object
        obj.run_samtools("sort",**{"-o":"out.bam","-@":"4","--":("in.bam",)})

*****************************
Examples and Case-Studies
*****************************


Many usage examples are available at https://github.com/urmi-21/pyrpipe/tree/master/case_studies


Transcript assembly
========================
A pyrpipe worflow for A thaliana transcript assembly

https://github.com/urmi-21/pyrpipe/tree/master/case_studies/Athaliana_transcript_assembly

LncRNA prediction
========================
A pyrpipe workflow for Maize lncRNA prediction.

https://github.com/urmi-21/pyrpipe/tree/master/case_studies/Maize_lncRNA_prediction

Integrating third-party tools
=============================
The module :py:mod:`pyrpipe_engine` contains helper methods to run any shell command from python.
See this examle for a detailed guide.

https://github.com/urmi-21/pyrpipe/blob/master/case_studies/Integrating%20third-party%20tools.ipynb
