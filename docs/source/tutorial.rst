Tutorial
=========

This tutorial provides introduction to pyrpipe. The first part describes using pyrpipe's RNA-Seq related modules to implement RNA-Seq processing pipelines.
The latter topics cover details on how to use pyrpipe framework to integrate any third-party UNIX command/tool in python.

.. contents::

Preparing environment
----------------------
Provide instruction XXXX
We recommend installing pyrpipe and all dependencies within a conda environment to ensure maximum reproducibility.
To create a new conda environment and switch

.. code-block:: bash
    
    conda create -n pyrpipe python=3.8
    conda activate pyrpipe


Basic RNA-Seq processing
-------------------------
After setting up the environment, one can import pyrpipe modules in python and start using it.
RNA-Seq processing can be as easy as creating required objects and executing required functions.


.. code-block:: python
    :linenos:
    
    from pyrpipe import sra,qc,mapping,assembly
    #define some vaiables
    working_dir='path_to_working_dir'
    genome='path_to_reference_genome_fasta'
    annotation='path_to_annotation'
    star_index='path_to_star_index'
    #initialize objects
    #creates a star object to use with 5 threads
    star=mapping.Star(index=star_index,threads=5)
    #use trim_galore for trimming
    trim_galore=qc.Trimgalore(threads=8)
    #Stringtie for assembly
    stringtie=assembly.Stringtie(guide=annotation,threads=10)
    #create SRA object and perform analysis
    sra.SRA('SRRXXXX',directory=working_dir).trim(trim_galore).align(star).assemble(stringtie)
    
Explanation of above code 
