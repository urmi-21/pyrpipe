pyrpipe: A python package for RNA-Seq workflows
================================================

:Author: Urminder Singh
:Date: |today|
:Version: |version|

*************
Introduction
*************
pyrpipe (pronounced as "pyre-pipe") is a python package to easily develop computational pipelines in pure python, in an object-oriented manner. 
pyrpipe provides an easy-to-use object oriented framework to import any UNIX command or tool in python. 
pyrpipe also provides flexible and easy handling of tool options and parameters and can automatically load the tool parameters from .yaml files.
This framework minimizes the commands user has to write, rather the tools are available as objects and are fully re-usable.
All commands executed via pyrpipe are automatically logged extensively and could be compiled into reports using  the pyrpipe_diagnostic tool.

To enable easy and easy processing of RNA-Seq data, we have implemented specialized classes build on top of the pyrpipe framework.
These classes provide high level APIs to many popular RNA-Seq tools for easier and faster development of RNA-Seq pipelines.
These pipelines can be fully customizable and users can easily add/replace the tools using the pyrpipe framework.
Finally, pyrpipe can be used on local computers or on HPC environments and pyrpipe scripts can be easily integrated into workflow management systems such as Snakemake and Nextflow.


Key Features
#############

- Import any UNIX command in python
- Dry-runnable pipelines to check dependencies and commands before execution
- Flexible and robust handling of tool arguments and parameters
- Auto load parameters from .yaml files
- Easily override threads and memory options using global values
- Extensive logging and reports with MultiQC reports for bioinformatic pipelines
- Specify GNU make like targets and verify the integrity of the targets
- Automatically resume pipelines/jobs from where interrupted
- Easily integrated into workflow managers like Snakemake and NextFlow


Installation
#############

To install the latest stable release via conda::
	conda install -c bioconda pyrpipe 

To install the latest stable release via pip::
	pip install pyrpipe --upgrade

Install latest development version ::
	git clone https://github.com/urmi-21/pyrpipe.git
	pip install -r pyrpipe/requirements.txt
	pip install -e path_to/pyrpipe

See the :ref:`Installation notes <installation>` for details.



*****************************
Examples and Case-Studies
*****************************
Example usage and case-studies with real data is provided at `GitHub <https://github.com/urmi-21/pyrpipe/tree/master/case_studies>`_.

 
Contents
--------

.. toctree::
   :maxdepth: 2

   installation.rst
   tutorial/index.rst
   cookbook.rst
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

