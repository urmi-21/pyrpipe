======================================================
pyrpipe: python RNA-Seq pipelines
======================================================

Introduction
============

Pysam is a python module that makes it easy to read and manipulate
mapped short read sequence data stored in SAM/BAM files.  It is a
lightweight wrapper of the htslib_ C-API.

This page provides a quick introduction in using pysam followed by the
API. See :ref:`usage` for more detailed usage instructions.

To use the module to read a file in BAM format, create a
:class:`~pysam.AlignmentFile` object::

   import pysam
   samfile = pysam.AlignmentFile("ex1.bam", "rb")

Once a file is opened you can iterate over all of the read mapping to
a specified region using :meth:`~pysam.AlignmentFile.fetch`.  Each
iteration returns a :class:`~pysam.AlignedSegment` object which
represents a single read along with its fields and optional tags::

   for read in samfile.fetch('chr1', 100, 120):
	print read

   samfile.close()
