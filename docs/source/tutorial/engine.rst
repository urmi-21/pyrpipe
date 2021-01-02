The pyrpipe_engine
---------------------
The pyrpipe_engine module is the key module responsible for handling execution of all the Unix commands.
The Runnable class calls the execute_comand() function in the pyrpipe_engine module.
All the commands executed via the pyrpipe_engine module are automatically logged.
All the functions responsible for executing Unix command are "decorated" with the dryable method, which allows using the --dry-run flag on any function using the pyrpipe_engine module.

The pyrpipe_engine can be directly used to execute Unix command or to import output of a Unix command in python. The functions defined in the pyrpipe_engine module are described below.


The pyrpipe_utils
---------------------
The pyrpipe_utils module define multiple helpful functions. Functions defined in the pyrpipe_utils modules are extensively used throughout pyrpipe modules.
Users can directly utilize these fuctions in their code and expedite development.
A description of these functions is provided below.


============        ====================
Function            Description
------------        -----------------
__init__()          This is the constructor. It can take SRR accession, path to fastq files, or sra file as input. If accession if provided as input the files are downloaded via prefetch if they aren't preset on disk. It will automatically handle single-end and paired-end data.
download_sra()      This function downloads the sra file via prefetch.