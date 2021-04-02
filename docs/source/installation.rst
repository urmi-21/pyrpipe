======================
Installation
======================
.. _installation:


**Note:** See :ref:`conda_setup` in tutorial, to learn how to install pyrpipe and required tools in an conda environment.

pyrpipe is available through bioconda could be installed via::

	conda install -c bioconda pyrpipe

pyrpipe is available on PyPI and could be installed via pip::

	pip install pyrpipe

To install from source, clone the git repo::

	git clone https://github.com/urmi-21/pyrpipe

Then install using pip::

	pip install -r pyrpipe/requirements.txt
	pip install -e pyrpipe

To run tests using pytest, from the pyrpipe root perform::

	py.test
	#or
	py.test tests/<specific test file>

