RNA-Seq API
------------

pyrpipe implements specialized classes for RNA-Seq processing. These classes are defined into different modules, each designed to capture steps integral to 
analysis of RNA-Seq data --from downloading raw data to trimming, alignment and assembly or quantification.
Each of these module implements classes coressponding to a RNA-Seq tool. These classes extend the `Runnable` class.
Specialized functions are implemented such that analysis of RNA-Seq data is intuitive and easy to code.
The following table provide details about the pyrpipe's RNA-Seq related modules.


=============================		=============================		====================================================================
	Module 					Class  					Purpose
-----------------------------		-----------------------------		--------------------------------------------------------------------
assembly    				Assembly			         Abstract class to represent Assebler type
assembly    				Stringtie			         API to Stringtie
assembly   				 Cufflinks       			 API to Cufflinks
mapping    				 Aligner            			 Abstract class for Aligner type
mapping     				Star                  			 API to Star
mapping     				Bowtie2               			  API to Bowtie2
mapping     				Hisat2                			  API to Hisat2
qc          				RNASeqQC              			  Abstract class for RNASeqQC type (quality control and trimming)
qc          				Trimgalore            			  API to Trim Galore
qc          				BBmap                 			  API to bbduk.sh
quant       				Quant                 			  Abstract class for Quantification type
quant       				Salmon                			  API to Salmon
quant       				Kallisto              			API to Kallisto
sra         				SRA                   			Class to represent RNA-Seq data and API to NCBI SRA-Tools
tools       				Samtools              			API to Samtools and other commonly used tools
=============================		=============================		====================================================================




The SRA class
^^^^^^^^^^^^^^^^^^
The SRA class contained in the sra module captures RNA-Seq data.
It can automatically download RNA-Seq data from the NCBI SRA servers via the prefetch command.
The SRA constructor can take the SRR accession, path to fastq or sra file as arguments.

The main attributes and functions are defined the following table.

=================       ===================================================================
Attribute        	   Description
-----------------       -------------------------------------------------------------------
fastq_path       	   Path to the fastq file. If single end this is the only fastq file.
fastq2_path      	   Path to the second fastq file for paired end data.
sra_path         	   Path to the sra file
srr_accession    	   The SRR accession for RNA-Seq run
layout           	   RNA-Seq layout, auto determined by SRA class.
bam_path         	   Path to bam file after running the align() function
gtf              	   Path to the gtf file after running assemble() function
=================       ===================================================================



================        ====================
Function                Description
----------------        --------------------
__init__()          	This is the constructor. It can take SRR accession, path to fastq files, or sra file as input. If accession if provided as input the files are downloaded via prefetch if they aren't preset on disk. It will automatically handle single-end and paired-end data.
download_sra()      	This function downloads the sra file via prefetch.
download_fastq()    	This function runs fasterq-dump on the sra file downloaded via prefetch.
sra_exists()        	Check if fastq sra files are present
fastq_exists()      	Check if fastq file exists
delete_sra()        	Delete the sra file
delete_fastq()      	Delete the fastq files
trim()              	This function takes a RNASeqQC type object and performs trimming. The trimmed fastq files are then stored in fastq_path and fastq2_path.
align()             	This function takes an Aligner type object and performs read alignemnt. The BAM file returned is stored in bam_path attribute.
assemble()          	This function takes an Assembly type object and performs transcript assembly. The result is soted on the SRA object as gtf attributes
quant()             	This function takes a Quantification type object and performs quant.
================        ====================



The RNASeqQC class
^^^^^^^^^^^^^^^^^^^^

The RNASeqQC is an abstract class defined in the qc module. RNASeqQC class extends the Runnable class and thus has all the attributes as in the Runnable class.
Classes Trimgalore and BBmap extends RNASeqQC class and share following attributes and functions.

============        ====================
Attribute            Description
------------        --------------------
_category           Represents the type: "RNASeqQC"
============        ====================

============        ====================
Function            Description
------------        --------------------
__init__()          The constructor function
perform_qc()        Takes a SRA object, performs qc and returns path to resultant fastq files
============        ====================



The Aligner class
^^^^^^^^^^^^^^^^^^^^^
The Aligner is an abstract class defined in the mapping module. Aligner class extends the Runnable class and thus has all the attributes as in the Runnable class.
Classes Star, Hisat2 and Bowtie2 extends the Aligner class and share following attributes and functions.

============        ====================
Attribute            Description
------------        --------------------
_category           Represents the type: "Aligner"
index               Index used by the aligner tool
genome              Reference genome used by the tool
============        ====================

===================         ===========================================
Function           		  Description
-------------------         -------------------------------------------
__init__()         		  The constructor function
build_index()      		  Build an index for the aligner tool using the `genome`.
check_index()      		  Checks if the index is valid
perform_alignment()		  Takes a sra object, performs alignemnt and returns path to the bam file
===================         ===========================================

The Assembly class
^^^^^^^^^^^^^^^^^^^^^
The Assembly is an abstract class defined in the assembly module. Assembly class extends the Runnable class and thus has all the attributes as in the Runnable class.
Classes Stringtie and Cufflinks extends the Assembly class and share following attributes and functions.

============        ====================
Attribute            Description
------------        --------------------
_category           Represents the type: "Assembler"
============        ====================


====================        ====================
Function            		Description
--------------------        --------------------
__init__()          		The constructor function
perform_assembly()  		Takes a SRA object, performs transcript assembly and returns path to resultant gtf/gff files
====================        ====================




The Quant class
^^^^^^^^^^^^^^^^^^
The Quant is an abstract class defined in the quant module. Quant class extends the Runnable class and thus has all the attributes as in the Runnable class.
Classes Salmon and Kallisto extends the Quant class and share following attributes and functions.


================        ====================
Attribute	            Description
----------------        --------------------
_category       	    Represents the type: "Quantification"
index           	    Index used by the aligner tool
transcriptome   	    Reference transcriptome used by the tool
================        ====================


================        ====================
Function	            Description
----------------        --------------------
__init__() 	         The constructor function
build_index()   	    Build an index for the quantification tool using the `transcriptome`.
check_index()   	    Checks if the index is valid
perform_quant() 	    Takes a sra object, performs quantification and returns path to the quantification results file
================        ====================





