pyrpipe: python RNA-Seq pipelines
==================================

:Author: Urminder Singh
:Date: |today|
:Version: |version|

*************
Introduction
*************
pyrpipe (Pronounced as "pyre-pipe") is a python package to easily develop RNA-Seq analyses workflow by integrating popular RNA-Seq analysis programs.
pyrpipe allows coding RNA-Seq workflows in pure python in an object oriented manner. pyrpipe provides a high level API to many popular RNA-Seq tools for easier and faster development.
pyrpipe can be used on local computers or on HPC environments to manage analysis of RNA-Seq data.

To install the latest release, type::

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

