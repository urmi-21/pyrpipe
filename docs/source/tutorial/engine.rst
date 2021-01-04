pyrpipe_engine module
------------------------
The pyrpipe_engine module is the key module responsible for handling execution of all the Unix commands.
The Runnable class calls the execute_comand() function in the pyrpipe_engine module.
All the commands executed via the pyrpipe_engine module are automatically logged.
All the functions responsible for executing Unix command are "decorated" with the dryable method, which allows using the --dry-run flag on any function using the pyrpipe_engine module.

The pyrpipe_engine can be directly used to execute Unix command or to import output of a Unix command in python. The functions defined in the pyrpipe_engine module are described below.

=========================       ====================
Function           		 Description
-------------------------       --------------------
dryable()           			A decorater to make Unix functions dry when --dry-run is specified
parse_cmd()         			Parse a Unix command and return command as string
get_shell_output()  			Execute a command and return return code and output. These commands are not logged
get_return_status()			 Execute command and return status of the command
execute_commandRealtime()  		 Execute a command and print output in realtime
execute_command()   			Execute a command and log the status. This is the function used by the Runnable class.
is_paired()         			Check is sra file paired or single end. Uses fastq-dump
get_program_path() 			 Return path to a Unix command
check_dependencies()   			 Check a list of dependencies are present in Unix path
delete_files()       			Delete files, rm command
move_file()         			Moves files, mv command
=========================       ====================




pyrpipe_utils module
---------------------
The pyrpipe_utils module define multiple helpful functions. Functions defined in the pyrpipe_utils modules are extensively used throughout pyrpipe modules.
Users can directly utilize these fuctions in their code and expedite development.
A description of these functions is provided below.


======================        ====================
Function            		Description
----------------------        --------------------
pyrpipe_print()     		Prints in color
get_timestamp()     		Return current timestamp
check_paths_exist() 		Return true is paths are valid
check_files_exist() 		Return True if files are valid
check_hisatindex()  		Verify valid Hisat2 index
check_kallistoindex()   	Verify kallisto index
check_salmonindex() 		Verify salmon index
check_starindex()   		Verify STAR index
check_bowtie2index()   		Verify Bowtie2 index
get_file_size()     		Return file size in human readable format
parse_java_args()   		Parse tool options in JAVA style format
parse_unix_args()   		Parse tool options in Unix style format
get_file_directory()    	Return a file's directory
get_filename()      		Return filename with extension
get_fileext()       		Return file extension
get_file_basename()     	Return filename, without extension
mkdir()             		Create a directory
get_union()         		Return union of lists
find_files()        		Search files using regex patterns
get_mdf()           		Compute and return MD5 checksum of a file
======================        ====================





