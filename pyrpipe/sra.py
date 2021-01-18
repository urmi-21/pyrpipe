#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 15:45:26 2019

@author: usingh


"""

from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe
import os
import statistics 
from pyrpipe import valid_args
from pyrpipe import param_loader as pl
from pyrpipe import _dryrun
from pyrpipe import _threads
from pyrpipe import _params_dir

class SRA:
    """This class represents an SRA object
    
        Parameters
        ----------
        
        srr_accession: string
            A valid SRR accession
        directory: string
            Path where all data related to this object (e.g. .sra files, metadata, fastq files) will be stored.
            Default value of the path will be "./<SRR_accession>". <SRR_accession> is added at the end of the path
            so that final directory is directory/<SRR_accession>.
            For consistency, directory and SRR Accession id are not allowed to be modified.
            
        scan_path: string
            If RNA-Seq data already exists locally, provide the scan path to scan a directory and create an SRA object.
        
        Attributes
        -----------
        
        """
    def __init__(self,srr_accession=None, directory=None, fastq=None,fastq2=None,sra=None):
    
        #these attributes can be associated with an sra object
        self.fastq_path=fastq
        self.fastq2_path=fastq2
        self.sra_path=sra
        self.srr_accession=srr_accession
        self.layout=None
        
        
        if not self.init_object(srr_accession, directory, fastq,fastq2,sra):
            pu.print_boldred("ERROR: Creating SRA object")
            raise ValueError("Please check fastq files {} {}".format(fastq,fastq2))
            
    def __setattr__(self, name, value):
        """
        Make srr_accession and directory immutable
        """
        immutableFields=['srr_accession','directory']
        
        if (name in immutableFields and name in self.__dict__):
            raise TypeError("Can not modify "+name)
        else:
            self.__dict__[name] = value
        
        
    def init_object(self,srr_accession,directory,fastq,fastq2,sra):
        
        #if fastq are provided
        if fastq and fastq2:
            self.layout="PAIRED"
            self.directory=pu.get_file_directory(fastq)
            if not pu.check_files_exist(fastq,fastq2) and not _dryrun:
                pu.print_boldred("ERROR: File not found")
                raise ValueError("Please check fastq files {} {}".format(fastq,fastq2))
            return True
                
        #if only one fastq (single)
        if fastq:
            self.layout="SINGLE"
            self.directory=pu.get_file_directory(fastq)
            if not pu.check_files_exist(fastq) and not _dryrun:
                pu.print_boldred("ERROR: File not found")
                raise ValueError("Please check fastq files {}".format(fastq))
            return True
    
        #init from srr_accession and directory
        return self.init_from_accession(srr_accession,directory)
    
    def init_from_accession(self,srr_accession,directory):
        """
        Create SRA object using provided srr accession and directory, where data is downloaded/saved
        This functions inits srrid, and paths to srr/fastq if they already exist thus will not be downloaded again
        """
        #check if programs exist
        self.dep_list=['prefetch',"fasterq-dump"]
        if not pe.check_dependencies(self.dep_list):
            raise OSError("ERROR: Please install missing programs.")
        
        if not srr_accession:
            raise ValueError("Please provide a valid accession")
            
        if not directory:
            directory=os.getcwd()
        
        #create a dir named <srr_accession> and use as directory
        self.directory=os.path.join(directory,self.srr_accession)
        
        #sra file be stored here
        #self.sra_path=os.path.join(self.directory,self.srr_accession+".sra")
        
        #check if fastq files exist
        if not self.search_fastq(self.directory):
            #download sra and run fqdump
            if not self.download_sra():
                pu.print_boldred('prefetch failed!!! Trying fasterq-dump...')
            #run fasterqdump either on downloaded SRA file or direclty
            return self.download_fastq()
        
        return True
    
    
        
    def search_fastq(self,path):
        """Search .fastq file under a dir and create SRA object
        Return True if found otherwise False
        """
        
        #check files with names <SRR>_1.fastq and <SRR>_2.fastq
        fq=os.path.join(path,self.srr_accession+'_1.fastq')
        fq2=os.path.join(path,self.srr_accession+'_2.fastq')
        if pu.check_files_exist(fq,fq2):
            self.fastq_path=fq
            self.fastq2_path=fq2
            pu.print_green("Found .fastq "+self.fastq_path+" "+self.fastq2_path)
            self.layout="PAIRED"
            return True
        
        #check single end file
        fq=os.path.join(path,self.srr_accession+'.fastq')
        if pu.check_files_exist(fq):
            self.fastq_path=fq
            pu.print_green("Found .fastq "+self.fastq_path)
            self.layout="SINGLE"
            return True       
        
        #search files under the path
        #fq_files=pe.find_files(path,"*.fastq")
        fq_files=pu.find_files(path,".fastq$")
        
        if len(fq_files)<1:
            return False
        
        if len(fq_files)>2:
            return False
        
        fq_files.sort()
        #case with single fastq
        if len(fq_files)==1:
            self.fastq_path=fq_files[0]
            pu.print_green("Found .fastq "+self.fastq_path)
            self.layout="SINGLE"
        
        #case with paired fastq
        if len(fq_files)==2:
            self.fastq_path=fq_files[0]
            self.fastq2_path=fq_files[1]
            pu.print_green("Found .fastq "+self.fastq_path+" "+self.fastq2_path)
            self.layout="PAIRED"
        
        return True 
    
                
    
    def download_fastq(self,*args,**kwargs):
        """Function to download fastq files
        """
        
        #check if fastq files exists already
        if self.fastq_exists():
            pu.print_green("Fastq files exist already")
            return True        
        #internal_args are created by pyrpipe and will always replace external passed args
        #add the positional args
        if self.sra_exists():
            internal_args=(self.sra_path,)
        else:
            #fstrqd_Cmd.append(self.srr_accession)
            internal_args=(self.srr_accession,)
        #keyword args; boolean flags have empty values
        internal_kwargs={'-O':self.directory,
                         '-o':self.srr_accession+".fastq",
                         '-e':_threads,
                         '-f':""
                         }
        
        
        #merge args, kwargs, internal_args, internal_kwargs
        #If args and kwargs are present
        
        if args or kwargs:
            internal_kwargs={**kwargs,**internal_kwargs}
            internal_args=tuple(set(args+internal_args))
            #append the args to the kwargs using special key '--'
            internal_kwargs['--']=internal_args
        else:
            #check for yaml parameters        
            filepath=os.path.join(_params_dir,'fasterq-dump.yaml')
            yaml_params=pl.YAML_loader(filepath)
            yaml_kwargs=yaml_params.get_kwargs()
            #yaml_args=yaml_params.get_args()
            internal_kwargs={**yaml_kwargs,**internal_kwargs}
            #internal_args=tuple(set(yaml_args+internal_args))
            internal_kwargs['--']=internal_args
        
        
        
        params_list=pu.parse_unix_args(valid_args._args_FASTERQDUMP,internal_kwargs)
        
        fstrqd_Cmd=['fasterq-dump']
        
        #add command and params
        fstrqd_Cmd.extend(params_list)
        
        #execute command
        cmdStatus=pe.execute_command(fstrqd_Cmd,objectid=self.srr_accession)        
        
        if not cmdStatus:
            pu.print_boldred("fasterqdump failed for:"+self.srr_accession)
            return False        
        
        #self.search_fastq(self.directory)
        #determine layout
        self.layout='PAIRED'
        #check files with names <SRR>_1.fastq and <SRR>_2.fastq
        fq=os.path.join(self.directory,self.srr_accession+'_1.fastq')
        fq2=os.path.join(self.directory,self.srr_accession+'_2.fastq')
        self.fastq_path=fq
        self.fastq2_path=fq2
        
        #if dry run
        if _dryrun: return True
        
        if pu.check_files_exist(fq,fq2):
            self.fastq_path=fq
            self.fastq2_path=fq2
            self.layout="PAIRED"
            #remove SRA
            self.delete_sra()
            return True
        
        #check single end file
        fq=os.path.join(self.directory,self.srr_accession+'.fastq')
        if pu.check_files_exist(fq):
            self.fastq_path=fq
            self.layout="SINGLE"
            #remove SRA
            self.delete_sra()
            return True
        
        return False
        
        
    def sra_exists(self):
        """Function to check if sra file is present on disk
        """
        if hasattr(self,'sra_path') and self.sra_path:
            
            return(os.path.isfile(self.sra_path))
        else:
            return False
        
    def fastq_exists(self):
        """Function to check if fastq file is present on disk
        """
        
        if not hasattr(self,'layout'):
            return False
        
        if self.layout=='PAIRED':
            if hasattr(self,'fastq_path'):
                if not self.fastq_path:
                    #if None
                    return False
                if not pu.check_files_exist(self.fastq_path):
                    return False
            else:
                return False
            
            if hasattr(self,'fastq2_path'):
                if not self.fastq2_path:
                    #if None
                    return False
                if not pu.check_files_exist(self.fastq2_path):
                    return False
            else:
                return False
            
            return True
            
        else:
            if hasattr(self,'fastq_path'):
                if not self.fastq_path:
                    #if None
                    return False
                return pu.check_files_exist(self.fastq_path)
            else:            
                return False
    
    def align(self,mapping_object,**kwargs):
        """
        
        """
        #check a valid mapping_object
        if not (hasattr(mapping_object,'_category') and mapping_object._category=='Aligner'):
            raise TypeError("Error: No valid mapping object provided for "+self.srr_accession)
            
        status=mapping_object.perform_alignment(self,objectid=self.srr_accession,**kwargs)
        if not status:
            raise OSError("perform_mapping failed for: "+ self.srr_accession)
        
        #save the bam file
        self.bam_path=status
        return self
    
    def assemble(self,assembly_object,**kwargs):
        """
        
        """
        #check a valid mapping_object
        if not (hasattr(assembly_object,'_category') and assembly_object._category=='Assembler'):
            raise TypeError("Error: No valid assembly object provided for "+self.srr_accession)
            
        #must have a bam file
        if not self.bam_path:
            raise ValueError("No BAM file associated with "+ self.srr_accession)
        
        if not _dryrun and not pu.check_files_exist(self.bam_path):
            raise ValueError("BAM file {} not fount ".format(self.bam_path))
            
        status=assembly_object.perform_assembly(self.bam_path,objectid=self.srr_accession,**kwargs)
        
        if not status:
            raise OSError("perform_mapping failed for: "+ self.srr_accession)
        
        self.gtf=status    
        return self
    
    
    def quant(self,quant_object):
        """
        
        """
        #check a valid mapping_object
        if not (hasattr(quant_object,'_category') and quant_object._category=='Quantification'):
            raise ValueError("Error: No valid assembly object provided for "+self.srr_accession)
            
            
                   
        status=quant_object.perform_quant(self,objectid=self.srr_accession)
        if not status:
            raise OSError("perform_mapping failed for: "+ self.srr_accession)
            
        self.abundance=status    
        return self
        
    
    def trim(self,qc_object,delete_original=False,**kwargs):
        """Function to perform quality control with specified qc object.
        A qc object refers to one of the RNA-Seq qc program like trim_galore oe bbduk.
        The qc_object should be initialized with all the parameters.
        By default the trimmed/qc fastq files will be generated in the same directory as the original fastq files.
        After QC, this SRA object will update the fastq_path or fastq_path and fastq2_path variables to store the new fastq files.
        New variables localRawfastqPath or rawfastq_path and rawfastq2_path will be created to store the paths of original fastq files.
      
        Parameters
        ----------
        
        qc_object: RNASeqQC object
            qc_object specifying the program to be used. The object contains the necessary parametrs to execute the parameters
            
        deleteRawFastq: bool
            Delete the raw fastq files after QC
        
        kwargs: dict
            Arguments to pass on to perform_qc function

        :return: Return status of the QC. True if successful download and False if failed.
        :rtype: bool

        Examples
        --------
        >>> object.perform_qc(qc.BBmap())
        True
        """
        #check a valid qc_object
        if not (hasattr(qc_object,'_category') and qc_object._category=='RNASeqQC'):
            raise ValueError("Error: No valid QC object provided for "+self.srr_accession)
        
        #save thq qc object for later references
        #self.qc_object=qc_object
        #print("Performing QC using "+qc_object.programName)
        #each qc_object has a function run() to execute their method
        qcStatus=qc_object.perform_qc(self,objectid=self.srr_accession,**kwargs)
        
        
        #if job failed
        if not qcStatus[0]:
            #print ("Error performing QC for "+self.srr_accession)
            #return False
            raise OSError("perform_qc failed for: "+ self.srr_accession)
            
        
        if self.layout=='PAIRED':
            
            #delete old fastq files if specified
            if delete_original:
                self.delete_fastq()
            
            else:
                #create new fields to refer to older fastq files
                self.rawfastq_path=self.fastq_path
                self.rawfastq2_path=self.fastq2_path
            
            ##update local fastq path
            self.fastq_path=qcStatus[0]    
            self.fastq2_path=qcStatus[1]
        else:
            if delete_original:
                self.delete_fastq()
            else:
                self.localRawfastqPath=self.fastq_path
                
            self.fastq_path=qcStatus[0]
        
        #return True
        return self
    
    def delete_fastq(self):
        """Delte the fastq files from the disk.
        The files are referenced by self.fastq_path or self.fastq_path and self.fastq2_path
        """
        
        if self.layout=='PAIRED':
            if pe.delete_files(self.fastq_path,self.fastq2_path):
                del self.fastq_path
                del self.fastq2_path
                return True
        else:
            
            if pe.delete_file(self.fastq_path):
                del self.fastq_path
                return True
        
        return False
        
        
    def delete_sra(self):
        """Delete the downloaded SRA files.
        """
        if not hasattr(self,'sra_path'):
            return True
        
        if not self.sra_path:
            return True
        if(pe.delete_file(self.sra_path)):
            del self.sra_path
            return True
        return False
    
    def get_read_length(self,lines_to_examine=None):
        """Examine first lines_to_examine lines and return the mode of read lengths
        returns int
        """
        if _dryrun:
            return 100
        
        if not lines_to_examine:
            lines_to_examine=4*10000
        else:
            lines_to_examine=4*lines_to_examine
            
        if not self.fastq_exists():
            pu.print_boldred("Fastq files don't exist")
            return 0
            
        #check the read len
        if self.layout=='SINGLE':
            
            with open(self.fastq_path) as f:
                data=[next(f) for x in range(lines_to_examine)]
            reads_lens=[]
            
            for i in range(1,len(data),4):
                reads_lens.append(len(data[i].strip()))
                
            return statistics.mode(reads_lens)
        
        if self.layout=='PAIRED':
            reads_lens=[]
            with open(self.fastq_path) as f:
                data=[next(f) for x in range(lines_to_examine)]
            #get lens from first file
            for i in range(1,len(data),4):
                reads_lens.append(len(data[i].strip()))
            mode1=statistics.mode(reads_lens)
            #print("Mode:",mode1)
            
            return mode1
            """
            Not sure how to deal with unqeual read lengths
            
            #read 2nd file
            reads_lens=[]
            with open(self.fastq2_path) as f:
                data=[next(f) for x in range(lines_to_examine)]
            #get lens from second file
            for i in range(1,len(data),4):
                reads_lens.append(len(data[i].strip()))
            
            mode2=statistics.mode(reads_lens)
            print("Mode:",mode2)
            
            return (mode1+mode2)/2
            """
        return 0
    
    
    def download_sra(self,**kwargs):
        """This function downloads .sra file from NCBI SRA servers using the prefetch command.

        NCBI sra-toolkit 2.9 or higher must be installed on the system in order to use prefetch. 
        prefetch will create a folder with name same as <srr_accession> under the directory (path) specified.
        The path of downloaded file is saved in the object as localSRAPath. This localSRAPath is then used
        by other functions to access the downloaded data. 
        The **kwargs is for passing arguments to the prefetch command.
        
        Parameters
        ----------
        
        kwargs: dict
            dict containing additional prefetch arguments

        :return: Return status of the prefetch command. True if successful download and False if failed.
        :rtype: bool

        Examples
        --------
        >>> object.download_sra()
        True
        """     
        
        #store path to the downloaded sra file
        self.sra_path=os.path.join(self.directory,self.srr_accession+".sra")
        #check if already exists
        if pu.check_files_exist(self.sra_path):
            #pu.print_green("File already exists:"+self.sra_path)
            #save file .sra file size
            self.sraFileSize=pu.get_file_size(self.sra_path)
            #test if file is paired or single end
            if pe.is_paired(self.sra_path):
                self.layout="PAIRED"
            else:
                self.layout="SINGLE"
            return True
                
        #scan for prefetch arguments
        prefetchArgsList=['-f','-t','-l','-n','-s','-R','-N','-X','-o','-a','--ascp-options','-p','--eliminate-quals','-c','-o','-O','-h','-V','-L','-v','-q']
        
        
        #ignore directory and file name arguments if given
        if '-O' in kwargs:
            print("Ignoring -O flag."+" directory is: "+self.directory)
            #delete -O parameter
            del kwargs['-O']
        if '-o' in kwargs:
            print("Ignoring -o flag."+" File name is: "+self.srr_accession)
            #delete -o parameter
            del kwargs['-o']
            

        prefetch_Cmd=['prefetch']
        prefetch_Cmd.extend(pu.parse_unix_args(prefetchArgsList,kwargs))
        prefetch_Cmd.extend(['-O',self.directory])
        prefetch_Cmd.append(self.srr_accession)
                
        cmdStatus=pe.execute_command(prefetch_Cmd,objectid=self.srr_accession)
        
        #return if dryrun
        if _dryrun: return True
        
        if not cmdStatus:
            pu.print_boldred("prefetch failed for:"+self.srr_accession)
            return False
        
        
        #move file if downloaded inside the directory
        if not pu.check_files_exist(self.sra_path):
            #check outdir/SRR/SRR/SRR.sra
            if pu.check_files_exist(os.path.join(self.directory,self.srr_accession,self.srr_accession+".sra")):
                pe.move_file(os.path.join(self.directory,self.srr_accession,self.srr_accession+".sra"),self.sra_path)
                
           
        #validate path exists
        if not pu.check_files_exist(self.sra_path):            
            pu.print_boldred("Error downloading file. File "+self.sra_path+" does not exist!!!\n Please check you SRA-Tools config")
            return False
        
        #print ("Downloaded file: "+self.sra_path+" {0} ".format(pu.get_file_size(self.sra_path)))
        #save file .sra file size
        self.sraFileSize=pu.get_file_size(self.sra_path)
        #test if file is paired or single end
        if pe.is_paired(self.sra_path):
            self.layout="PAIRED"
        else:
            self.layout="SINGLE"
            
        
        return True
    
    
    #def destroy(self):
        """
        Delete everything for this object from memory and disk
        """
      #  pass
    
        