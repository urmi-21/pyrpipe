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


class SRA:
    """This class represents an SRA object
    
        Parameters
        ----------
        
        srr_accession: string
            A valid SRR accession
        location: string
            Path where all data related to this object (e.g. .sra files, metadata, fastq files) will be stored.
            Default value of the path will be "./<SRR_accession>". <SRR_accession> is added at the end of the path
            so that final location is location/<SRR_accession>.
            For consistency, location and SRR Accession id are not allowed to be modified.
            
        scan_path: string
            If RNA-Seq data already exists locally, provide the scan path to scan a directory and create an SRA object.
        
        Attributes
        -----------
        
        """
    def __init__(self,srr_accession=None, location=None, threads=None):
        
        if location is None and srr_accession is None:
            #pu.print_boldred("Please provide a valid srr accession or location, or both")
            raise Exception("Please provide a valid srr accession or location, or both")
            #use the scan_path to create sra object
            #self.init_from_path(scan_path)
        #else:
            #use the provided srr_accession
        self.init_from_accession(srr_accession,location)
        
        #default threads to use
        if not threads:
            threads=os.cpu_count()
        self.threads=str(threads)
    
    def init_from_path(self,path):
        if not pu.check_paths_exist(path):
            raise Exception("Please provide a valid path to scan for RNA-Seq data")
        
        #scan path for fastq
        self.search_fastq(path)
        #scan path for sra
        self.search_sra(path)
        if not (self.fastqFilesExistsLocally() or self.sraFileExistsLocally()):
                raise Exception("No files found at:"+ path+ "Please provide a valid path to scan for RNA-Seq data")
    
    def search_sra(self,path):
        """Search .sra file under a dir
        Return True if found otherwise False
        """
        #search files under the path
        
        sra_files=pe.find_files(path,"*.sra")
        
        if len(sra_files)<1:
            return False
        
        if len(sra_files)>1:
            pu.print_boldred("Found multiple .sra files. Using the first entry...")
        sra_path=sra_files[0]
        #self.location=path
        self.srr_accession=pu.get_file_basename(sra_path)
        self.localSRAFilePath=sra_path
        self.sraFileSize=pu.get_file_size(self.localSRAFilePath)
        #test if file is paired or single end
        if pe.is_paired(self.localSRAFilePath):
            self.layout="PAIRED"
        else:
            self.layout="SINGLE"
        
        pu.print_green("Found .sra "+self.localSRAFilePath)
        return True
        
    def search_fastq(self,path):
        """Search .fastq file under a dir and create SRA object
        Return True if found otherwise False
        """
        #search files under the path
        fq_files=pe.find_files(path,"*.fastq")
        
        if len(fq_files)<1:
            return False
        
        if len(fq_files)>2:
            pu.print_boldred("Can not determine .fastq. Downloading...")
            return False
        
        fq_files.sort()
        #case with single fastq
        if len(fq_files)==1:
            self.localfastqPath=fq_files[0]
            pu.print_green("Found .fastq "+self.localfastqPath)
            self.layout="SINGLE"
        
        #case with paired fastq
        if len(fq_files)==2:
            self.localfastq1Path=fq_files[0]
            self.localfastq2Path=fq_files[1]
            pu.print_green("Found .fastq "+self.localfastq1Path+" "+self.localfastq2Path)
            self.layout="PAIRED"
        
        #self.location=path
        #self.srr_accession=pu.get_file_basename(fq_files[0])
        return True
        
    
    def init_from_accession(self,srr_accession,location):
        """Create SRA object using provided srr accession and location to save the data
        """
        self.dep_list=['prefetch',"fasterq-dump"]
        if not pe.check_dependencies(self.dep_list):
            raise Exception("ERROR: Please install missing programs.")
        
        if srr_accession is None:
            raise Exception("Please provide a valid accession")
            
        if location is None:
            location=os.getcwd()
        #pu.print_info("Creating SRA: "+srr_accession)
        self.srr_accession=srr_accession
        #create a dir named <srr_accession> and use as location
        self.location=os.path.join(location,self.srr_accession)
    
        #search for existing files in location
        #self.search_fastq(self.location)
        #scan path for sra
        #self.search_sra(self.location)
        
        
        #check SRA file
        if pu.check_files_exist(os.path.join(self.location,self.srr_accession+".sra")):
            pu.print_green(self.srr_accession+".sra exists.")
            self.localSRAFilePath=os.path.join(self.location,self.srr_accession+".sra")
            self.sraFileSize=pu.get_file_size(self.localSRAFilePath)
            #test if file is paired or single end
            if pe.is_paired(self.localSRAFilePath):
                self.layout="PAIRED"
            else:
                self.layout="SINGLE"
        
        #check fastq file
        self.search_fastq(self.location)
        
        
        
        
    def __setattr__(self, name, value):
        """Make srr accession immutable
        """
        immutableFields=['srr_accession','location']
        
        if (name in immutableFields and name in self.__dict__):
            raise Exception("Can not modify "+name)
        else:
            self.__dict__[name] = value
                
    
    def download_fastq(self,verbose=False,quiet=False,logs=True,threads=None,**kwargs):
        """Function to download fastq files
        """
        
        #check if fastq files exists already
        if self.fastqFilesExistsLocally():
            pu.print_green("Fastq files exist already")
            return True
        
        #overwrite threads
        if not threads:
            threads=self.threads
        else:
            threads=str(threads)
        
        
        fasterqdumpArgsList=['-f','-t','-s','-N','-X','-a','-p','-c','-o','-O','-h','-V',
                             '-L','-v','-q','-b','-m','-x','-S','-3','-P','-M',
                             '-B','--option-file','--strict','--table','--include-technical',
                             '--skip-technical','--concatenate-reads']
        
        fstrqd_Cmd=['fasterq-dump']
        fstrqd_Cmd.extend(pu.parse_unix_args(fasterqdumpArgsList,kwargs))
        #add location
        fstrqd_Cmd.extend(['-O',self.location])
        #add output filename. output will be <srr_accession>.fastq or <srr_accession>_1.fastq and <srr_accession>_2.fastq
        fstrqd_Cmd.extend(['-o',self.srr_accession+".fastq"])
        fstrqd_Cmd.extend(['-e',str(threads)])
        fstrqd_Cmd.extend(['-f'])
        if self.sraFileExistsLocally():
            fstrqd_Cmd.append(self.localSRAFilePath)
        else:
            fstrqd_Cmd.append(self.srr_accession)
        
        #execute command
        cmdStatus=pe.execute_command(fstrqd_Cmd,objectid=self.srr_accession)
        if not cmdStatus:
            print("fasterqdump failed for:"+self.srr_accession)
            return False        
        
        
        if not hasattr(self,'layout'):
            fq_files=pe.find_files(self.location,self.srr_accession+"*.fastq")
            if len(fq_files)==1:
                self.layout='SINGLE'
            else:
                self.layout='PAIRED'
        
        #check if fastq files are downloaded        
        if(self.layout=="SINGLE"):
            self.localfastqPath=os.path.join(self.location,self.srr_accession+".fastq")
            if not pu.check_files_exist(self.localfastqPath):
                pu.print_boldred("Error running fasterq-dump file. File "+self.localfastqPath+" does not exist!!!")
                return False
        else:
            self.localfastq1Path=os.path.join(self.location,self.srr_accession+"_1.fastq")
            self.localfastq2Path=os.path.join(self.location,self.srr_accession+"_2.fastq")
            if not pu.check_files_exist(self.localfastq1Path,self.localfastq2Path):
                pu.print_boldred("Error running fasterq-dump file. File "+self.localfastq1Path+" does not exist!!!")
                return False
        
            
            
            
        return True
        
        
  
    def download_sra(self,verbose=False,quiet=False,logs=True,**kwargs):
        """This function downloads .sra file from NCBI SRA servers using the prefetch command.

        NCBI sra-toolkit 2.9 or higher must be installed on the system in order to use prefetch. 
        prefetch will create a folder with name same as <srr_accession> under the location (path) specified.
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
        self.localSRAFilePath=os.path.join(self.location,self.srr_accession+".sra")
        #check if already exists
        if pu.check_files_exist(self.localSRAFilePath):
            pu.print_green("File already exists:"+self.localSRAFilePath)
            #save file .sra file size
            self.sraFileSize=pu.get_file_size(self.localSRAFilePath)
            #test if file is paired or single end
            if pe.is_paired(self.localSRAFilePath):
                self.layout="PAIRED"
            else:
                self.layout="SINGLE"
            return True
            
        
        pu.print_info("Downloading "+self.srr_accession+" ...")
        
        #scan for prefetch arguments
        prefetchArgsList=['-f','-t','-l','-n','-s','-R','-N','-X','-o','-a','--ascp-options','-p','--eliminate-quals','-c','-o','-O','-h','-V','-L','-v','-q']
        
        
        #ignore location and file name arguments if given
        if '-O' in kwargs:
            print("Ignoring -O flag."+" location is: "+self.location)
            #delete -O parameter
            del kwargs['-O']
        if '-o' in kwargs:
            print("Ignoring -o flag."+" File name is: "+self.srr_accession)
            #delete -o parameter
            del kwargs['-o']
            

        prefetch_Cmd=['prefetch']
        prefetch_Cmd.extend(pu.parse_unix_args(prefetchArgsList,kwargs))
        prefetch_Cmd.extend(['-O',self.location])
        prefetch_Cmd.append(self.srr_accession)
                
        cmdStatus=pe.execute_command(prefetch_Cmd,objectid=self.srr_accession)
        if not cmdStatus:
            pu.print_boldred("prefetch failed for:"+self.srr_accession)
            return False
        
        
        #validate path exists
        if not pu.check_files_exist(self.localSRAFilePath):
            pu.print_boldred("Error downloading file. File "+self.localSRAFilePath+" does not exist!!!")
            return False
        
        print ("Downloaded file: "+self.localSRAFilePath+" {0} ".format(pu.get_file_size(self.localSRAFilePath)))
        #save file .sra file size
        self.sraFileSize=pu.get_file_size(self.localSRAFilePath)
        #test if file is paired or single end
        if pe.is_paired(self.localSRAFilePath):
            self.layout="PAIRED"
        else:
            self.layout="SINGLE"
            
        
        return True
    
    
           
    
    def sraFileExistsLocally(self):
        """Function to check if sra file is present on disk
        """
        if hasattr(self,'localSRAFilePath'):
            return(os.path.isfile(self.localSRAFilePath))
        else:
            return False
        
    def fastqFilesExistsLocally(self):
        """Function to check if fastq file is present on disk
        """
        
        if not hasattr(self,'layout'):
            return False
        
        if self.layout=='PAIRED':
            if hasattr(self,'localfastq1Path'):
                if not pu.check_files_exist(self.localfastq1Path):
                    return False
            else:
                return False
            
            if hasattr(self,'localfastq2Path'):
                if not pu.check_files_exist(self.localfastq2Path):
                    return False
            else:
                return False
            
            return True
            
        else:
            if hasattr(self,'localfastqPath'):
                return pu.check_files_exist(self.localfastqPath)
            else:            
                return False
    
    def run_fasterqdump(self,delete_sra=False,verbose=False,quiet=False,logs=True,**kwargs):
        """Execute fasterq-dump to convert .sra file to fastq files.
        The fastq files will be stored in the same directory as the sra file. All fastq files should be consistently named
        using the extension .fastq
        
        Parameters
        ----------
        
        delete_sra: bool
            delete sra file after completion
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        kwargs: dict
            A dict containing fasterq-dump arguments
        
        :return: Return status of the fasterq-dump command. True if successful download and False if failed.
        :rtype: bool

        Examples
        --------
        >>> object.run_fasterqdump()
        True
        """
        #check if fastq files exists already
        if self.fastqFilesExistsLocally():
            pu.print_green("Fastq files exist already")
            return True
        
        #first check is sra exists
        if not self.sraFileExistsLocally():
            pu.print_boldred("Error executing fasterq-dump: .sra file not found. Please run download_sra().")
            return False
        #else directly run fasterq-dump on accession ?
        
        fasterqdumpArgsList=['-f','-t','-s','-N','-X','-a','-p','-c','-o','-O','-h','-V',
                             '-L','-v','-q','-b','-m','-e','-x','-S','-3','-P','-M',
                             '-B','--option-file','--strict','--table','--include-technical',
                             '--skip-technical','--concatenate-reads']
        
        
        
        #ignore location and file name arguments if given
        if '-O' in kwargs:
            print("Ignoring -O flag."+" location is: "+self.location)
            #delete -O parameter
            del kwargs['-O']
        if '-o' in kwargs:
            print("Ignoring -o flag."+" File name is: "+self.srr_accession)
            #delete -o parameter
            del kwargs['-o']
        
        
        #determine threads if not specified
        if '-e' not in kwargs:
            kwargs['-e']=str(self.threads)
        #execute command
        
        fstrqd_Cmd=['fasterq-dump']
        fstrqd_Cmd.extend(pu.parse_unix_args(fasterqdumpArgsList,kwargs))
        
        #add location
        fstrqd_Cmd.extend(['-O',self.location])
        #add output filename. output will be <srr_accession>.fastq or <srr_accession>_1.fastq and <srr_accession>_2.fastq
        fstrqd_Cmd.extend(['-o',self.srr_accession+".fastq"])
        fstrqd_Cmd.append(self.localSRAFilePath)
        
        #execute command
        cmdStatus=pe.execute_command(fstrqd_Cmd,objectid=self.srr_accession)
        if not cmdStatus:
            print("fasterqdump failed for:"+self.srr_accession)
            return False
        
        
        #check if fastq files are downloaded 
        if(self.layout=="SINGLE"):
            self.localfastqPath=os.path.join(self.location,self.srr_accession+".fastq")
            
            if not pu.check_files_exist(self.localfastqPath):
                pu.print_boldred("Error running fasterq-dump file. File "+self.localfastqPath+" does not exist!!!")
                return False
            
        else:
            self.localfastq1Path=os.path.join(self.location,self.srr_accession+"_1.fastq")
            self.localfastq2Path=os.path.join(self.location,self.srr_accession+"_2.fastq")
            
            if not pu.check_files_exist(self.localfastq1Path,self.localfastq2Path):
                pu.print_boldred("Error running fasterq-dump file. File "+self.localfastq1Path+" does not exist!!!")
                return False
            
        #delete sra file if specified
        if delete_sra:
            self.delete_sra()
            
        return True
    
    def perform_qc(self,qcObject,deleteRawFastq=False,**kwargs):
        """Function to perform quality control with specified qc object.
        A qc object refers to one of the RNA-Seq qc program like trim_galore oe bbduk.
        The qcObject should be initialized with all the parameters.
        By default the trimmed/qc fastq files will be generated in the same directory as the original fastq files.
        After QC, this SRA object will update the localfastqPath or localfastq1Path and localfastq2Path variables to store the new fastq files.
        New variables localRawfastqPath or localRawfastq1Path and localRawfastq2Path will be created to store the paths of original fastq files.
      
        Parameters
        ----------
        
        qcObject: RNASeqQC object
            qcObject specifying the program to be used. The object contains the necessary parametrs to execute the parameters
            
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
        #check a valid qcObject
        if not (hasattr(qcObject,'category') and qcObject.category=='RNASeqQC'):
            print ("Error: No valid QC object provided. Skipping QC for "+self.srr_accession)
            return False
        
        #save thq qc object for later references
        self.QCObject=qcObject
        print("Performing QC using "+qcObject.programName)
        #each qcObject has a function run() to execute their method
        qcStatus=qcObject.perform_qc(self,objectid=self.srr_accession,**kwargs)
        
        #if job failed
        if not qcStatus[0]:
            print ("Error performing QC for "+self.srr_accession)
            return False
        
        if self.layout=='PAIRED':
            
            #delete old fastq files if specified
            if deleteRawFastq:
                self.delete_fastq()
            
            else:
                #create new fields to refer to older fastq files
                self.localRawfastq1Path=self.localfastq1Path
                self.localRawfastq2Path=self.localfastq2Path
            
            ##update local fastq path
            self.localfastq1Path=qcStatus[0]    
            self.localfastq2Path=qcStatus[1]
        else:
            if deleteRawFastq:
                self.delete_fastq()
            else:
                self.localRawfastqPath=self.localfastqPath
                
            self.localfastqPath=qcStatus[0]
        
        
        
        return True
    
    def delete_fastq(self):
        """Delte the fastq files from the disk.
        The files are referenced by self.localfastqPath or self.localfastq1Path and self.localfastq2Path
        """
        if self.layout=='PAIRED':
            if pe.deleteMultipleFilesFromDisk(self.localfastq1Path,self.localfastq2Path):
                del self.localfastq1Path
                del self.localfastq2Path
                return True
        else:
            if pe.deleteFileFromDisk(self.localfastqPath):
                del self.localfastqPath
                return True
        return False
        
        
    def delete_sra(self):
        """Delete the downloaded SRA files.
        """
        if(pe.deleteFileFromDisk(self.localSRAFilePath)):
            del self.localSRAFilePath
            return True
        return False
    
    def get_read_length(self,lines_to_examine=None):
        """Examine first lines_to_examine lines and return the mode of read lengths
        
        """
        if not lines_to_examine:
            lines_to_examine=4*10000
        else:
            lines_to_examine=4*lines_to_examine
            
        if not self.fastqFilesExistsLocally():
            pu.print_boldred("Fastq files don't exist")
            return 0
            
        #check the read len
        if self.layout=='SINGLE':
            
            with open(self.localfastqPath) as f:
                data=[next(f) for x in range(lines_to_examine)]
            reads_lens=[]
            
            for i in range(1,len(data),4):
                reads_lens.append(len(data[i].strip()))
                
            return statistics.mode(reads_lens)
        
        if self.layout=='PAIRED':
            reads_lens=[]
            with open(self.localfastq1Path) as f:
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
            with open(self.localfastq2Path) as f:
                data=[next(f) for x in range(lines_to_examine)]
            #get lens from second file
            for i in range(1,len(data),4):
                reads_lens.append(len(data[i].strip()))
            
            mode2=statistics.mode(reads_lens)
            print("Mode:",mode2)
            
            return (mode1+mode2)/2
            """
            
        return 0
            
        
    
