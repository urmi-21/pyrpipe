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
from pyrpipe import valid_args
from pyrpipe import _threads


class RNASeqTools(Runnable):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.category="RNASeqTools"
        self._command=None
        
    

class Samtools(RNASeqTools):
    """
    Class to access samtools
    args: tuple
        arguments to samtools
    threads: int
        Number of threads samtools will use
    kwargs: dict
        Options to samtools
    
    """
    def __init__(self,*args,threads=None,**kwargs):
        super().__init__(*args,**kwargs)
        self._command='samtools'
        self._deps=[self._command]
        self._param_yaml='samtools.yaml'
        self._valid_args=valid_args._args_SAMTOOLS
        #resolve threads to use
        self.resolve_parameter("-@",threads,_threads,'_threads')

    def sam_to_bam(self,sam_file,out_dir=None,out_suffix=None, delete_sam=False,objectid=None):
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
                
        if not out_suffix:
            out_suffix=""
        
        fname=pu.get_file_basename(sam_file)
        
        
        #target: output will be out_bam
        out_bam=os.path.join(out_dir,fname+out_suffix+'.bam')
        internal_args=(sam_file,)
        internal_kwargs={"-o":out_bam,"-b":""}
        
        
        status=self.run(*internal_args,subcommand="view",target=out_bam,objectid=objectid,**internal_kwargs)
                
        if not status:
            return ""
        
        #delete_sam_file
        if delete_sam:
            if not pe.deleteFileFromDisk(sam_file):
                pu.print_boldred("Error deleting sam file:"+sam_file)
                
        #return path to file
        return out_bam
        
        
        
        
    #sort bam file.output will be bam_file_sorted.bam
    def sort_bam(self,bam_file,out_dir=None,out_suffix=None,delete_bam=False,objectid=None):
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
        if not out_suffix:
            out_suffix=""
                
        fname=pu.get_file_basename(bam_file)
        #output will be out_bam
        outSortedbam_file=os.path.join(out_dir,fname+out_suffix+'_sorted.bam')
        
        
        internal_args=(bam_file,)
        
        internal_kwargs={"-o":outSortedbam_file}
        
        
        status=self.run(*internal_args,subcommand="sort",target=outSortedbam_file,objectid=objectid,**internal_kwargs)
        
        if not status:
            print("Bam sort failed for:"+bam_file)
            return ""
        

        if delete_bam:
            if not pe.deleteFileFromDisk(bam_file):
                print("Error deleting sam file:"+bam_file)
                
        #return path to file
        return outSortedbam_file
    
    def sam_sorted_bam(self,sam_file,out_dir=None,out_suffix=None,delete_sam=False,delete_bam=False,objectid=None):
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
        
            
        sam2bam_file=self.sam_to_bam(sam_file,delete_sam=delete_sam,objectid=objectid)
        
        if not sam2bam_file:
            return ""
            

        bamSorted=self.sort_bam(sam2bam_file,out_dir, out_suffix,delete_bam=delete_bam,objectid=objectid)
        
        if not bamSorted:
            return ""
        
        return bamSorted
    
    
    def merge_bam(self,bam_list,out_file="merged",out_dir=None,delete_bams=False,objectid=None):
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
               
        if len(bam_list) < 2:
            pu.print_boldred("Error: Please supply at least 2 files to merge")
            return ""
        
        if not out_dir:
            out_dir=pu.get_file_directory(bam_list[0])
        else:
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
        
        outMergedFile=os.path.join(out_dir,out_file+".bam")
        
        internal_args=(outMergedFile,)+tuple(bam_list)        
        internal_kwargs={}
        
        
        status=self.run(*internal_args,subcommand="merge",target=outMergedFile,objectid=objectid,**internal_kwargs)
        
        if not status:
            print("Bam merge failed for:"+outMergedFile)
            return ""
        
        #check if bam file exists
        if not pu.check_files_exist(outMergedFile):
            return ""
        

        if delete_bams:
            for bam_file in bam_list:
                if not pe.deleteFileFromDisk(bam_file):
                    print("Error deleting sam file:"+bam_file)
                    
        return outMergedFile
        
    
    
    