pyrpipe: python RNA-Seq pipelines
==================================

:Author: Urminder Singh
:Date: |today|
:Version: |version|

*************
Introduction
*************
pyrpipe (pronounced as "pyre-pipe") is a python package to easily develop computational pipelines in pure python. pyrpipe provides an easy-to-use object oriented framework to import any UNIX command or tool in python. All commands executed via pyrpipe are automatically logged extensively. pyrpipe also provides automatic handling of tool options and parameters.
pyrpipe comes with specialized classes for easy implementation of RNA-Seq analysis pipelines. These classes provide high level APIs to many popular RNA-Seq tools for easier and faster development of pipelines.
pyrpipe can be used on local computers or on HPC environments to manage analysis of RNA-Seq data.


Key Features
#############

- Import any UNIX command in python
- Dry-run feature to check dependencies and commands before execution
- Flexible and robust handling of options and arguments
- Auto load command options from .yaml files
- Easily override threads and memory options using global values
- Extensive logging for all the commands
- Automatically verify Integrity of output targets
- Resume feature to restart pipelines/jobs from where interrupted
- Create reports, MultiQC reports for bioinformatic pipelines
- Easily integrated into workflow managers like Snakemake and NextFlow



To install the latest stable release, type::
	
	conda install -c bioconda pyrpipe 
	or
	pip install pyrpipe --upgrade

Install latest development version ::

	git clone https://github.com/urmi-21/pyrpipe.git
	pip install -r pyrpipe/requirements.txt
	pip install -e path_to/pyrpipe

See the :ref:`Installation notes <installation>` for details.


*****************************
Examples and Case-Studies
*****************************
Example usage and case-studies with real data is provided at `GitHub <https://github.com/urmi-21/pyrpipe/tree/master/examples(case-studies)>`_.

- Arabidopsis thaliana transcript assembly using pyrpipe checkpoints. `Link <https://github.com/urmi-21/pyrpipe/tree/master/examples(case-studies)/Athaliana_transcript_assembly>`_.
- Prediction of long non-coding RNAs (lncRNAs) in Zea mays by supplementing pyrpipe with a third-party tool. `Link <https://github.com/urmi-21/pyrpipe/tree/master/examples(case-studies)/Maize_lncRNA_prediction>`_.
- Integrating pyrpipe scripts within a workflow management system (Snakemake). `Link <https://github.com/urmi-21/pyrpipe/tree/master/examples(case-studies)/Human_annotation_snakemake>`_.
- Guide for integrating third-party tools into pyrpipe. `Link <https://github.com/urmi-21/pyrpipe/blob/master/examples(case-studies)/Integrating%20third-party%20tools.ipynb>`_.
 

***********
Philosophy
***********
pyrpipe provides dedicated high-level APIs to bioinformatic tools for easy development of workflows.
pyrpipe follows an object oriented approach in which the tools, its parameters and data are encapsulated as classes.
Objects of these classes can be the used to build workflows. Each class representing a tool has a `run()` method to invoke the tool by providing all necessary parameters.
Other methods, depending on the tool, are provided which allows integration and execution with other tools or objects. For example, `perform_alignment()` method of `Hisat2` class 
takes an `SRA` object as input and produces `SAM` files.

pyrpipe makes it easy to write RNA-Seq workflows in python. It also provides addtional features like logging, report generation, and benchmarking. 
pyrpipe is not a full workflow management system and does not directly manage parallel execution, memory managemet, and scaling to cluster. 
pyrpipe scripts can be easily integrated with systems like snakemake to acheive more scalable workflows.

Accessing tools with pyrpipe 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Contents
--------

.. toctree::
   :maxdepth: 2

   installation.rst
   usage.rst
   api.rst
..  faq.rst
..  developer.rst
..  release.rst
..  benchmarking.rst
..  glossary.rst

Indices and tables
------------------

Contents:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

References
----------

