#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 14:54:22 2019

@author: usingh
"""
import os
from pyrpipe.runnable import Runnable
from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe



class RNASeqTools(Runnable):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.category="RNASeqTools"
        self._command=None
        
    

class Samtools(RNASeqTools):
    """init function allows to specify commonly tuned parameters like threads, out_dir, memory etc.
    More specific parameters can be passed when calling the specific functions. There is option to override the parameters specified here in init later.
    
    threads: int
        Number of threads samtools will use
    max_memory: int
        Max memory to use in GB
    
    """
    def __init__(self,threads=None,max_memory=None):
        self.programName="samtools"
        #check if hisat2 exists
        if not pe.check_dependencies([self.programName]):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        self.threads=threads
        #Default: if threads are None use 80% of threads to avaoid memory issues
        if not self.threads:
            self.threads=int(os.cpu_count()*0.8)
            
        self.max_memory=max_memory
        
        
    def sam_to_bam(self,sam_file,out_dir="",out_suffix="",threads=None, delete_sam=False,objectid="NA",**kwargs):
        """Convert sam file to a bam file. 
        Output bam file will have same name as input sam.
        sam_file: string
            Path to input Sam file
        out_suffix: string
            Suffix for the output sam file
        threads: int
            Number of threads. Default: Use self.threads initialized in init().
        delete_sam: bool
            delete the sam file after conversion
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to samtools. This will override the existing options 

        :return: Returns the path to the bam file. Returns empty string if operation failed.
        :rtype: string
        """        
        if not out_dir:            
            out_dir=pu.get_file_directory(sam_file)
        else:
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
        
        fname=pu.get_file_basename(sam_file)
        
        #output will be out_bam
        out_bam=os.path.join(out_dir,fname+out_suffix+'.bam')
        
        #handle threads
        if not threads:
            threads=self.threads
        
        newOpts={"--":(sam_file,),"-o":out_bam,"-@":str(threads),"-b":""}
        
        #add (and override) any arguments provided via kwargs
        mergedOpts={**newOpts,**kwargs}
        
        status=self.run_samtools("view",objectid=objectid,**mergedOpts)
                
        if not status:
            print("Sam to bam failed for:"+sam_file)
            return ""
        
        #check if bam file exists
        if not pu.check_files_exist(out_bam):
            return ""
        
        #delete_sam_file
        if delete_sam:
            if not pe.deleteFileFromDisk(sam_file):
                print("Error deleting sam file:"+sam_file)
                
        #return path to file
        return out_bam
        
        
        
        
    #sort bam file.output will be bam_file_sorted.bam
    def sort_bam(self,bam_file,out_dir="",out_suffix="",threads=None,delete_bam=False,objectid="NA",**kwargs):
        """Sorts an input bam file. Outpufile will end in _sorted.bam
        bam_file: str
            Path to the input bam file
        out_dir: str
            Path to output directory
        out_suffix: str
            Output file suffix
        threads: int
            Number of threads. Default: Use self.threads initialized in init().
        delete_bam: bool
            Delete input bam_file
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to samtools. This will override the existing options 

        :return: Returns path to the sorted bam file. Returns empty string if operation failed.
        :rtype: string
        
        """
        if not out_dir:
            out_dir=pu.get_file_directory(bam_file)
        else:
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
                
        fname=pu.get_file_basename(bam_file)
        #output will be out_bam
        outSortedbam_file=os.path.join(out_dir,fname+out_suffix+'_sorted.bam')
        
        #handle threads
        if not threads:
            threads=self.threads
        
        newOpts={"--":(bam_file,),"-o":outSortedbam_file,"-@":str(threads)}
        mergedOpts={**newOpts,**kwargs}
        
        status=self.run_samtools("sort",objectid=objectid,**mergedOpts)
        
        if not status:
            print("Bam sort failed for:"+bam_file)
            return ""
        
        #check if bam file exists
        if not pu.check_files_exist(outSortedbam_file):
            return ""

        if delete_bam:
            if not pe.deleteFileFromDisk(bam_file):
                print("Error deleting sam file:"+bam_file)
                
        #return path to file
        return outSortedbam_file
    
    def sam_sorted_bam(self,sam_file,out_dir="",out_suffix="",threads=None,delete_sam=False,delete_bam=False,objectid="NA",**kwargs):
        """Convert sam file to bam and sort the bam file.
        sam_file: str
            Path to the input sam file
        out_dir: str
            Path to output directory
        out_suffix: str
            Output file suffix
        threads: int
            Number of threads. Default: Use self.threads initialized in init().
        delete_sam: bool
            Delete input sam_file
        delete_bam: bool
            Delete the intermediate unsorted bam_file
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to samtools. This will override the existing options 

        :return: Returns path to the sorted bam file. Returns empty string if operation failed.
        :rtype: string
        """
        
        sam2bam_file=self.sam_to_bam(sam_file,threads=threads,delete_sam=delete_sam,objectid=objectid,**kwargs)
        
        if not sam2bam_file:
            return ""
            

        bamSorted=self.sort_bam(sam2bam_file,out_dir, out_suffix,threads=threads,delete_bam=delete_bam,objectid=objectid,**kwargs)
        
        if not bamSorted:
            return ""
        
        return bamSorted
    
    
    def merge_bam(self,*args,out_file="merged",out_dir="",threads=None,delete_bams=False,objectid="NA",**kwargs):
        """Merge multiple bam files into a single file
        
        Parameters
        ----------
        args: *args
            Input bam files to merge
        out_file: string
            Output file name to save the results. .bam will be added at the end.
        out_dir: string
            Path where to save the merged bam file. Default path is the same as the first bam_file's
        threads: int
            Number of threads. Default: Use self.threads initialized in init().
        delete_bams: bool
            Delete input bam files after merging.
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to samtools. This will override the existing options 

        :return: Returns the path to the merged bam file.
        :rtype: string
        """
               
        if len(args) < 2:
            print("Please supply at least 2 files to merge")
            return ""
        
        if not out_dir:
            out_dir=pu.get_file_directory(args[0])
        else:
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
        
        outMergedFile=os.path.join(out_dir,out_file+".bam")
        
        #handle threads
        if not threads:
            threads=self.threads
        
        newOpts={"-@":str(threads),"--":(outMergedFile,)+args}
        
        #override parameters by supplying as kwargs
        mergedOpts={**newOpts,**kwargs}
        
        status=self.run_samtools("merge",objectid=objectid,**mergedOpts)
        
        if not status:
            print("Bam merge failed for:"+outMergedFile)
            return ""
        
        #check if bam file exists
        if not pu.check_files_exist(outMergedFile):
            return ""
        

        if delete_bams:
            for bam_file in args:
                if not pe.deleteFileFromDisk(bam_file):
                    print("Error deleting sam file:"+bam_file)
                    
        return outMergedFile
        
        
        
    def run_samtools(self,sub_command,valid_args=None,objectid="NA",**kwargs):
        """A wrapper to run samtools.
        
        Parameters
        ----------
        
        sub_command: string
            sub_command to pass to samtools e.g. sort, merge etc
        valid_args: list
            A list containing valid parameters. Parameters in kwargs not in this list will be ignored. Default: None
        arg1: dict
            arguments to pass to samtools. 
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to samtools. This will override the existing options 

        :return: Returns the status of samtools. True is passed, False if failed.
        :rtype: bool
        """
        samtools_cmd=['samtools',sub_command]
        #add options
        samtools_cmd.extend(pu.parse_unix_args(valid_args,kwargs))
                
        #start ececution
        status=pe.execute_command(samtools_cmd,objectid=objectid)
        if not status:
            pu.print_boldred("samtools failed")
        
        #return status
        return status
    
    
    