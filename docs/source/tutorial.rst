Tutorial
=========

This tutorial provides introduction to pyrpipe. The first part describes using pyrpipe's RNA-Seq related modules to implement RNA-Seq processing pipelines.
The latter topics cover details on how to use pyrpipe framework to integrate any third-party UNIX command/tool in python.

.. contents::

Preparing environment
----------------------
We recommend installing pyrpipe and all dependencies within a conda environment to ensure maximum reproducibility.
To create a new conda environment and switch

.. code-block:: bash
    
    conda create -n pyrpipe python=3.8
    conda activate pyrpipe

Basic RNA-Seq processing
-------------------------
After setting up the environment, one can import pyrpipe modules in python and start using it.
RNA-Seq processing can be as easy as


.. code-block:: python
    :linenos:
    
    from pyrpipe import sra,qc,mapping,assembly
    #create objects to use
    star=mapping.Star(index='path_to_star_index',threads=5)
     
