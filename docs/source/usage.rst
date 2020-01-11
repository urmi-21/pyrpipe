======================================================
Basic usage
======================================================

Specifying RNA-Seq data
=======================

Use the :py:mod:`sra` module to create an :class:`SRA` object::

	from pyrpipe import sra
	sra_obj = sra.SRA(srr_accession="SRR976159")

After the object is created, download the raw data to disk using the 
:meth:`SRA.download_sra` method::
	sra_obj.download_sra()
This will download the raw data in .sra format.
To convert .sra file to fastq, use :meth:`SRA.run_fasterqdump` method::

	sra_obj.run_fasterqdump()

The sra_obj will keep track of all the downloaded data. The location of downloaded data could be accesed by::

	sra_obj.location

To get the paths to sra of fastq files, use::

	sra_obj.localSRAFilePath
	sra_obj.localfastqPath
	#for paired
	sra_obj.localfastq1Path
	sra_obj.localfastq1Path


Performing read alignment
=========================
The :py:mod:`mapping` module contains several classes to access read alignment tools.
The method :meth:`perform_alignment` can be used with the SRA object.
An example using Hisat2::
	hs=mapping.Hisat2(hisat2_index="",**hsOpts)
	hs.perform_alignment(ob,**{"-p":"10"})
