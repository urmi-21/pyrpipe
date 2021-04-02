
Setting up the environment
==========================

First step is to setup the environment in a way such that reproducibility is ensured or maximized.
In this tutorial we will use the conda environment manager to install python, pyrpipe and the required tools and dependencies into a single environment.
**Note:** Conda must be installed on the system. For help with setting up conda, please see `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_.


.. _conda_setup:

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


If the above command fails, please try adding conda channels (see commands below) in the right order and then try again.

.. code-block:: bash
    
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge



Using conda environment in yaml files
#####################################

We have provided a yaml file containing the conda packages required to reproduce pyrpipe environment. Users can also use this file to create a conda environment and run pyrpipe.
To create a conda environment, use the `pyrpipe_environment.yaml <https://github.com/urmi-21/pyrpipe/blob/master/pyrpipe_environment.yaml>`_:

.. code-block:: bash

    conda env create -f pyrpipe_environment.yml


Users can easily export and share their own conda environment yaml files containing information about the conda environment.
To export any conda environment as yaml, run the following command

.. code-block:: bash

    conda env export | grep -v "^prefix: " > environment.yml

To recreate the conda environment in the `environment.yml`, use

.. code-block:: bash

    conda env create -f environment.yml

Automated installation of required tools
########################################

We have also provided a utility to install required RNA-Seq tools via a single command:

.. code-block:: bash
    
    pyrpipe_diagnostic build-tools
    
**Note:** Users must verify the versions of the tools installed in the conda environment.


Setting up NCBI SRA-Tools
------------------------------

After installing sra-tools, please configure prefetch to save the downloads the the **public user-repository**.
This will ensure that the prefetch command will download the data to the user defined directory.
To do this

- Type `vdb-config -i` command in terminal to open the NCBI SRA-Tools configuration editor.
- Under the TOOLS tab, set prefetch downloads option to **public user-repository**

Users can easily test if SRA-Tools has been setup properly by invoking the following command

.. code-block:: bash
    
    pyrpipe_diagnostic test



