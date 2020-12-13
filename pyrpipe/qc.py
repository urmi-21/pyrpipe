#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 17:48:00 2019

@author: usingh
"""
import os
from pyrpipe.runnable import Runnable
from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe
from pyrpipe import valid_args
from pyrpipe import _dryrun
from pyrpipe import _threads


class RNASeqQC(Runnable):
    """This is an abstract parent class for fastq quality control programs.
    """
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self._category="RNASeqQC"
        self._command=None
        
    def perform_qc(self):
        """Perform qc on a SRA object
        """
        pass

class Trimgalore(RNASeqQC):
    """This class represents trimgalore
    
    Parameters
    ----------
        
        kwargs:
            trim_galore arguments.
    """
    def __init__(self,*args,**kwargs):
        #run super to inherit parent class properties
        super().__init__(*args,**kwargs)
        self._command='trim_galore'
        self._deps=[self._command,'cutadapt']
        self._param_yaml='trim_galore.yaml'
        self._valid_args=valid_args._args_TRIM_GALORE
        self.check_dependency()
        self.init_parameters(*args,**kwargs)
            
    def perform_qc(self,sra_object,out_dir="",out_suffix="_trimgalore",objectid="NA"):
        """Function to perform qc using trimgalore.
        The function perform_qc() is consistent for all QC classess.
        
        Parameters
        ----------
        
        sra_object: SRA
            An SRA object whose fastq files will be used
        out_dir: str
            Path to output directory
        out_suffix: string
            Suffix for the output sam file
        threads: int
            Num threads to use
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to trimgalore. This will override the existing options 

            :return: Returns the path of fastq files after QC. tuple has one item for single end files and two for paired.
            :rtype: tuple
        """
        if not out_dir:
            out_dir=sra_object.directory
        else:
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
        
        #get layout
        if sra_object.layout=='PAIRED':
            fq1=sra_object.fastq_path
            fq2=sra_object.fastq2_path
            out_file1=os.path.join(out_dir,pu.get_file_basename(fq1)+out_suffix+".fastq")
            out_file2=os.path.join(out_dir,pu.get_file_basename(fq2)+out_suffix+".fastq")
            internal_args=(fq1,fq2)
            internal_kwargs={"--paired":"","-o":out_dir,"--cores":_threads}

            
            #run trimgalore
            self.run(*internal_args,objectid=objectid,**internal_kwargs)
            """
            running trim galore will create two files named <input>_val_1.fq and <input>_val_2.fq
            move these files to the specified out files
            """
            oldFile1=os.path.join(out_dir,pu.get_file_basename(fq1)+"_val_1.fq")
            oldFile2=os.path.join(out_dir,pu.get_file_basename(fq2)+"_val_2.fq")
            
            pe.move_file(oldFile1,out_file1,verbose=False)
            pe.move_file(oldFile2,out_file2,verbose=False)
            
            if not _dryrun and not pu.check_files_exist(out_file1,out_file2):
                print("Trimgalore failed")
                return ("",)
            return out_file1,out_file2
            
        else:
            fq=sra_object.localfastqPath
            out_file=os.path.join(out_dir, pu.get_file_basename(fq)+out_suffix+".fastq")
            internal_args=(fq,)
            internal_kwargs={"-o":out_dir,"--cores":_threads}

            #run trimgalore
            self.run(*internal_args,objectid=objectid,**internal_kwargs)
            """
            running trim galore will create one file named <input>_trimmed.fq
            move these files to the specified out files
            """
            oldFile=os.path.join(out_dir,pu.get_file_basename(fq)+"_trimmed.fq")
            
            pe.move_file(oldFile,out_file)
            
            if not _dryrun and not pu.check_files_exist(out_file):
                print("Trimgalore failed")
                return ("",)
            
            return (out_file,)
            
            

class BBmap(RNASeqQC):
    """This class represents bbmap programs
    """
    def __init__(self,*args,**kwargs):
        """
        """
        #run super to inherit parent class properties
        super().__init__(*args,**kwargs)
        self._command='bbduk.sh'
        self._args_style='JAVA'
        self._deps=[self._command]
        self._param_yaml='bbmap.yaml'
        self._valid_args=valid_args._args_BBDUK
        self.check_dependency()
        self.init_parameters(*args,**kwargs)
            
            
    def perform_qc(self,sra_object,out_dir="",out_suffix="_bbduk",objectid="NA"):
        """Run bbduk on fastq files specified by the sra_object
       
        :return: Returns the path of fastq files after QC. tuple has one item for single end files and 2 for paired.
        :rtype: tuple
            
        """
        #make out_dir
        if not out_dir:
                out_dir=sra_object.directory
        else:
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
               
        if sra_object.layout=='PAIRED':
            fq1=sra_object.fastq_path
            fq2=sra_object.fastq2_path
            out_fileName1=pu.get_file_basename(fq1)+out_suffix+".fastq"
            out_fileName2=pu.get_file_basename(fq2)+out_suffix+".fastq"
            out_file1Path=os.path.join(out_dir,out_fileName1)
            out_file2Path=os.path.join(out_dir,out_fileName2)
            
            internal_args=()
            internal_kwargs={"in":fq1,"in2":fq2,"out":out_file1Path,"out2":out_file2Path,"threads":_threads}
                        
            #run bbduk
            status=self.run(*internal_args,objectid=objectid,**internal_kwargs)
            
            if status:
                if not pu.check_files_exist(out_file1Path,out_file2Path) and not _dryrun:
                        return("",)
                        
            return(out_file1Path,out_file2Path)
            
            
        else:
            fq=sra_object.localfastqPath
            out_fileName=pu.get_file_basename(fq)+out_suffix+".fastq"
            out_filePath=os.path.join(out_dir,out_fileName)
            internal_args=()
            internal_kwargs={"in":fq,"out":out_filePath,"threads":_threads}
            
            #run bbduk
            status=self.run(*internal_args,objectid=objectid,**internal_kwargs)
            if status:
                if not pu.check_files_exist(out_file1Path,out_file2Path) and not _dryrun:
                    return("",)
                
            return(out_filePath,) 
    
    def perform_cleaning(self,sra_object,bbsplit_index,out_dir="",out_suffix="_bbsplit",overwrite=True,objectid="NA",**kwargs):
        """
        Remove contaminated reads mapping to given reference using bbsplit
        
        Parameters
        ----------
        
        sra_object: SRA
            an SRA object
        bbsplit_index: string
            Path to bbsplit index or fasta file which will generate index
        out_dir: string
            Path to output dir. Default: sra_object.directory
        out_suffix: string
            Suffix for output file name
        overwrite: bool
            overwrite existing files
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        
        kwargs: dict
            options passed to bbsplit

            :return: Returns the path of fastq files after QC. tuple has one item for single end files and 2 for paired.
            :rtype: tuple
        """
        
        #check index
        indexPath=""
        if not pu.check_paths_exist(bbsplit_index):
            #index folder doesn't exist
            #check if input is path to fasta
            if not pu.check_files_exist(bbsplit_index):
                print("Error: Please check bbsplit index")
                return ("",)
            #check if index folder "ref" exists in this directory
            indexPath=os.path.join(pu.get_file_directory(bbsplit_index),"ref")
            if pu.check_paths_exist(indexPath):
                print("Using bbsplit index: "+indexPath)
            else:
                #create new index
                print("Creating new index"+indexPath)
                newOpts={"ref_x":bbsplit_index,"path": pu.get_file_directory(bbsplit_index)}
                mergedOpts={**kwargs,**newOpts}
                #run bbduk
                if not self.run_bbsplit(objectid=objectid,**mergedOpts):
                    print("Error creating bbsplit index.")
                    return ("",)
                if not pu.check_paths_exist(indexPath):
                    print("Error creating bbsplit index.")
                    return ("",)
        else:
            indexPath=bbsplit_index
                
        
        #indexPath point to the ref directory, go one directory higher
        indexPath=os.path.dirname(indexPath)
        
        
        #make out_dir
        if not out_dir:
                out_dir=sra_object.directory
        else:
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
        
        if sra_object.layout=='PAIRED':
            fq1=sra_object.fastq_path
            fq2=sra_object.fastq2_path
            #append input and output options
            
            out_fileName1=pu.get_file_basename(fq1)+out_suffix+".fastq"
            out_fileName2=pu.get_file_basename(fq2)+out_suffix+".fastq"
            out_file1Path=os.path.join(out_dir,out_fileName1)
            out_file2Path=os.path.join(out_dir,out_fileName2)
            
            newOpts={"in1":fq1,"in2":fq2,"outu1":out_file1Path,"outu2":out_file2Path,"path":indexPath}
            mergedOpts={**kwargs,**newOpts}
            
            #run bbsplit
            if self.run_bbsplit(objectid=objectid,**mergedOpts):
                if pu.check_files_exist(out_file1Path,out_file2Path):
                    return(out_file1Path,out_file2Path)
            return("",)
            
            
        else:
            fq=sra_object.localfastqPath
            #append input and output options
           
            out_fileName=pu.get_file_basename(fq)+out_suffix+".fastq"
            out_filePath=os.path.join(out_dir,out_fileName)
            newOpts={"in":fq,"outu":out_filePath,"path":indexPath}
            mergedOpts={**kwargs,**newOpts}
            
            #run bbsplit
            if self.run_bbsplit(objectid=objectid,**mergedOpts):
                if pu.check_files_exist(out_filePath):
                    return(out_filePath,)
            
            return("",)
    
    
    
    def run_bbsplit(self,objectid="NA",**kwargs):
        """wrapper to run bbsplit
        
        :return: Status of bbsplit command
        :rtype: bool
        """
        
        bbsplit_args=['ref','ref_x','build','path','in','in1','in2','outu','outu2','outu1','qin','interleaved',
                          'maxindel','minratio','minhits','ambiguous','ambiguous2',
                          'qtrim','untrim','out_','basename','bs','scafstats',
                          'refstats','nzo','-Xmx','-eoom','-da']
        
        
        #create command to run
        bbsp_cmd=["bbsplit.sh"]
        
        #bbduk.sh follows java style arguments
        bbsp_cmd.extend(pu.parse_java_args(bbsplit_args,kwargs))
        
        
        #start ececution
        status=pe.execute_command(bbsp_cmd,objectid=objectid)
        if not status:
            pu.print_boldred("bbsplit failed")
        #return status
        return status
    
    
    
    
    
    
        
