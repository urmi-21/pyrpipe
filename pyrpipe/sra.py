#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 15:45:26 2019

@author: usingh


"""

from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe
import os


class SRA:
    """This class represents an SRA object
        Parameters
        ----------
        srr_accession(string): A valid SRR accession
        location (String): Path where all data related to this object (e.g. .sra files, metadata, fastq files) will be stored. 
                Default value of the path will be "./<SRR_accession>". <SRR_accession> is added at the end of the path
                so that final location is location/<SRR_accession>.
                For consistency, location and SRR Accession id are not allowed to be modified.
        
        Attributes
        -----------
        
        """
    def __init__(self,srr_accession,location=os.getcwd()):
        
        
        self.dep_list=['prefetch',"fasterq-dump"]
        if not pe.check_dependencies(self.dep_list):
            raise Exception("ERROR: Please install missing programs.")
        
        
        pu.print_info("Creating SRA: "+srr_accession)
        self.srr_accession=srr_accession
        #append the SRR accession to the location
        self.location=os.path.join(location,self.srr_accession)
        
        ##check if sra, fastq files already exists
    
        
    def __setattr__(self, name, value):
        """Make srr accession immutable
        """
        immutableFields=['srr_accession','location']
        
        if (name in immutableFields and name in self.__dict__):
            raise Exception("Can not modify "+name)
        else:
            self.__dict__[name] = value
                
    
  
    def download_sra(self,verbose=False,quiet=False,logs=True,**kwargs):
        """This function downloads .sra file from NCBI SRA servers using the prefetch command.

        NCBI sra-toolkit 2.9 or higher must be installed on the system in order to use prefetch. 
        prefetch will create a folder with name same as <srr_accession> under the location (path) specified.
        The path of downloaded file is saved in the object as localSRAPath. This localSRAPath is then used
        by other functions to access the downloaded data. 
        The **kwargs is for passing arguments to the prefetch command.
        
        Parameters
        ----------
        kwargs(dict): dict containing additional prefetch arguments
        
        Returns
        -------
        bool
            Return status of the prefetch command. True if successful download and False if failed.

        Examples
        --------
        >>> object.download_sra()
        True
        """
        
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
        
        #store path to the downloaded sra file
        self.localSRAFilePath=os.path.join(self.location,self.srr_accession+".sra")
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
        if hasattr(self,'localSRAFilePath'):
            return(os.path.isfile(self.localSRAFilePath))
        else:
            return False
        
    def fastqFilesExistsLocally(self):
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
        arg1: bool
            delete sra file after completion
            
        arg2: dict
            A dict containing fasterq-dump arguments
        
        Returns
        -------
        bool
            Return status of the fasterq-dump command. True if successful download and False if failed.

        Examples
        --------
        >>> object.run_fasterqdump()
        True
        """
        
        #first check is sra exists
        if not self.sraFileExistsLocally():
            print ("Error executing fasterq-dump: .sra file not found. Please run download_sra().")
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
            print("prefetch failed for:"+self.srr_accession)
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
    
    def perform_qc(self,qcObject,deleteRawFastq=False):
        """Function to perform quality control with specified qc object.
        A qc object refers to one of the RNA-Seq qc program like trim_galore oe bbduk.
        The qcObject should be initialized with all the parameters.
        By default the trimmed/qc fastq files will be generated in the same directory as the original fastq files.
        After QC, this SRA object will update the localfastqPath or localfastq1Path and localfastq2Path variables to store the new fastq files.
        New variables localRawfastqPath or localRawfastq1Path and localRawfastq2Path will be created to store the paths of original fastq files.
      
        Parameters
        ----------
        arg1: object
            qcObject specifying the program to be used. The object contains the necessary parametrs to execute the parameters
            
        arg2: bool
            Delete the raw fastq files after QC
        
        Returns
        -------
        bool
            Return status of the QC. True if successful download and False if failed.

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
        qcStatus=qcObject.perform_qc(self,objectid=self.srr_accession)
        
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
        
        

if __name__ == "__main__":
    #test
    print("main")
        
    
    
    