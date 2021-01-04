Setting up the environment
==========================

First step is to setup the environment in a way such that reproducibility is ensured or maximized.
In this tutorial we will use the conda environment manager to install python, pyrpipe and the required tools and dependencies into a single environment.
**Note:** Conda must be installed on the system. For help with setting up conda, please see `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_.



Create a new conda environment
-------------------------------
To create a new conda environment with python 3.8 execute the following commands.
We recommend sharing conda environment files with pipeline scripts to allow for reproducible analysis.

.. code-block:: bash
    
    conda create -n pyrpipe python=3.8

Activate the newly created conda environment and install required tools

.. code-block:: bash

    conda activate pyrpipe
    conda install -c bioconda pyrpipe star=2.7.7a sra-tools=2.10.9 stringtie=2.1.4 trim-galore=0.6.6 orfipy=0.0.3 salmon=1.4.0
    
To create a yaml file containing information about the conda environment, run the following command

.. code-block:: bash

    conda env export | grep -v "^prefix: " > environment.yml

To recreate the conda environment in the `environment.yml`, use

.. code-block:: bash

    conda env create -f environment.yml
	
