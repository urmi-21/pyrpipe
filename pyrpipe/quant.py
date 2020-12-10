#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 18:20:41 2020

@author: usingh
"""

from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe
from pyrpipe import valid_args
from pyrpipe import param_loader as pl
from pyrpipe import _dryrun
from pyrpipe import _threads
from pyrpipe import _params_dir
import os


class Quant:
    """This is an abstract class for quantification programs.
    """
    def __init__(self,index=""):
        self.category="Quantification"
        self.passedArgumentDict={}
        self.index=index
        
    def build_index(self):
        """function to create an index used by the quantification program
        """
        pass
    
    def check_index(self):
        """Function to check if index of this object is valid and exists
        """
    
    def perform_quant(self,sra_object):
        """Function to perform quant taking and sraobject as input
        
        """
        pass
    
class Kallisto(Quant):
    """This class represents kallisto
    
    index: string
        path to kallisto index
    threads: int
        num threads to use
    """
    
    def __init__(self,index,*args,transcriptome=None,**kwargs):
        super().__init__() 
        self.programName="kallisto"
        self.index=index
        self.transcriptome=transcriptome
        self.dep_list=[self.programName]        
        if not pe.check_dependencies(self.dep_list):
            raise Exception("ERROR: "+ self.programName+" not found.")
  
        #init the parameters for the object
        if args:
            self._args=args
        else:
            self._args=()
        if kwargs:
            self._kwargs=kwargs
        else:
            self._kwargs={}
        #read yaml parameters
        yamlfile=os.path.join(_params_dir,'kallisto.yaml')
        #override yaml parameters by kwargs
        if pu.check_files_exist(yamlfile):
            yaml_params=pl.YAML_loader(yamlfile)
            yaml_kwargs=yaml_params.get_kwargs()
            self._kwargs={**yaml_kwargs,**self._kwargs}
               
        #check index
        #kallisto index is a single file
        if not pu.check_files_exist(index):
            if not (pu.check_files_exist(self.transcriptome)):
                raise Exception("Please provide a valid Kallisto index or fasta file to generate the index")
            else:
                #call build index to generate index
                self.build_index(index,self.transcriptome)
                
            
    def build_index(self,index_path,transcriptome,overwrite=False,verbose=False,quiet=False,logs=True,objectid="NA"):
        """Function to  build kallisto index
        
        index_path: str
            path to the output directory
        index_name: str
            index name
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to kallisto. This will override the existing options in self.passed_args_dict (only replace existing arguments and not replace all the arguments).
            
        :return: Status of kallisto index
        :rtype: bool
        """
        #if index already exists then exit
        if not overwrite:
            #check if files exists
            if pu.check_files_exist(index_path):
                pu.print_green("Kallisto index {} already exists.".format(index_path))
                self.index=index_path
                return True
            
        #check input
        if not pu.check_files_exist(transcriptome):
            pu.print_boldred("{} does not exist. Exiting".format(transcriptome))
            raise Exception("Please check input to kallisto index")
            
        
        #create out dir
        indexdir=pu.get_file_directory(index_path)
        #create the out dir
        if not pu.check_paths_exist(indexdir):
            if not pu.mkdir(indexdir):
                raise Exception("Error creating kallisto index. Failed to create index directory.")
        
        args=(transcriptome,)
        internal_kwargs={"-i":index_path}
        #read build parameters
        yamlfile=os.path.join(_params_dir,'kallisto_index.yaml')
        if pu.check_files_exist(yamlfile):
            yaml_params=pl.YAML_loader(yamlfile)
            yaml_kwargs=yaml_params.get_kwargs()
            internal_kwargs={**yaml_kwargs,**internal_kwargs}
        
        #add positional args
        internal_kwargs['--']=args
                
        kallisto_cmd=['kallisto','index']
        kallisto_cmd.extend(pu.parse_unix_args(None,internal_kwargs))
        
        #call kallisto
        status=pe.execute_command(kallisto_cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
                
        if status:
            if pu.check_files_exist(index_path) and not _dryrun:
                #update object's index
                self.index=index_path
                if self.check_index():
                    return True
        else:
            raise Exception("Error building kallisto index")
        
        return False
    
    
    def perform_quant(self,sra_object,out_suffix="",out_dir="",verbose=False,quiet=False,logs=True,objectid="NA"):
        """Run kallisto quant
        
        sra_object: SRA
            SRA object contatining paths to fastq files
        index_path: str
            path to the output directory
        index_name: str
            index name
        threads: int
            Number of threads
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to kallisto. This will override the existing options

        :return: Path to kallisto out directory
        :rtype: string
        """
        
        if not out_dir:
            out_dir=os.path.join(sra_object.directory,"kallisto_out")
        else:
            #create out_dir if not exists
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
        
        
        if sra_object.layout == 'PAIRED':
            args=(sra_object.fastq_path,sra_object.fastq2_path)
            internal_kwargs={"--threads":_threads,"-o":out_dir,"-i":self.index}
        else:
            args=(sra_object.fastq_path,)
            internal_kwargs={"--threads":_threads,"-o":out_dir,"--single":"","-i":self.index}
            
                      
        #call kallisto
        status=self.run(*args,verbose=verbose,quiet=quiet,logs=logs,objectid=sra_object.srr_accession,**internal_kwargs)
        
        if status:
            outfile=os.path.join(out_dir,"abundance.tsv")
            if not pu.check_files_exist(outfile) and not _dryrun:
                return ""
            
            #move
            newfile=os.path.join(out_dir,"abundance"+out_suffix+".tsv")
            pe.move_file(outfile,newfile)
            return newfile
        
        return ""
        
    
    
    def run(self,*args,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper for running kallisto.
        
        Parameters
        ----------
        
        subcommand: str
            subcommand for kallisto
        valid_args: list
            List of valid arguments, arguments in kwargs not in this list will be ignored
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to kallisto. This will override the existing options

        :return: Returns the status of kallisto. True is passed, False if failed.
        :rtype: bool
        """
        
        #check for a valid index
        if not self.check_index() and not _dryrun:
            raise Exception("ERROR: Invalid kallisto index. Please run build index to generate an index.")
        
        #override class kwargs by passed
        kwargs={**self._kwargs,**kwargs}
        #if no args provided use constructor's
        if not args[0]:
            args=self._args
        #if args exist
        if args:
            #add args
            kwargs['--']=args
        
        kallisto_Cmd=['kallisto','quant']
        kallisto_Cmd.extend(pu.parse_unix_args(valid_args._args_KALLISTO,kwargs))
        
        #start ececution
        status=pe.execute_command(kallisto_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,command_name=kallisto_Cmd[0])
        if not status:
            pu.print_boldred("kallisto failed: "+ " ".join(kallisto_Cmd))
            
        return status       
    
    def check_index(self):
        """Check valid kallisto index
        """
        if hasattr(self,'index'):
            return(pu.check_files_exist(self.index))
        return False
            


class Salmon(Quant):
    """This class represents salmon
    
    index: string
        Path to salmon index
    threads: int
        Number of threads
    """      
    def __init__(self,index,*args,transcriptome=None,**kwargs):    
        super().__init__() 
        self.programName="salmon"
        self.index=index
        self.transcriptome=transcriptome
        self.dep_list=[self.programName]        
        if not pe.check_dependencies(self.dep_list):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
       
        
        #init the parameters for the object
        if args:
            self._args=args
        else:
            self._args=()
        if kwargs:
            self._kwargs=kwargs
        else:
            self._kwargs={}
        #read yaml parameters
        yamlfile=os.path.join(_params_dir,'salmon.yaml')
        #override yaml parameters by kwargs
        if pu.check_files_exist(yamlfile):
            yaml_params=pl.YAML_loader(yamlfile)
            yaml_kwargs=yaml_params.get_kwargs()
            self._kwargs={**yaml_kwargs,**self._kwargs}
               
        #check index
        if not pu.check_salmonindex(index):
            if not (pu.check_files_exist(self.transcriptome)):
                raise Exception("Please provide a valid Salmon index or fasta file to generate the index")
            else:
                #call build index to generate index
                self.build_index(index,self.transcriptome)
            
            
    def build_index(self,index_path,transcriptome,overwrite=False,verbose=False,quiet=False,logs=True,objectid="NA"):
        """
        build salmon index and store the path to index in self
        
        index_path: str
            path to the output directory
        index_name: str
            index name
        fasta: str
            Path to fasta file
        threads: int
            Number of threads
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to salmon. This will override the existing options
            
        :return: status of salmon index
        :rtype: bool
        """
        
        #if index already exists then exit
        if not overwrite:
            #check if files exists
            if pu.check_salmonindex(index_path):
                pu.print_green("Salmon index {} already exists.".format(index_path))
                self.index=index_path
                return True
            
        #check input
        if not pu.check_files_exist(transcriptome):
            pu.print_boldred("{} does not exist. Exiting".format(transcriptome))
            return False
        
        #create out dir
        indexdir=pu.get_file_directory(index_path)
        #create the out dir
        if not pu.check_paths_exist(indexdir):
            if not pu.mkdir(indexdir):
                raise Exception("Error creating salmon index. Failed to create index directory.")
        
        
        validArgsIndex=['-v','--version','-h','--help','-t','--transcripts','-k','--kmerLen','-i',
                             '--index','--gencode','--keepDuplicates','-p','--threads','--perfectHash',
                             '--type','-s','--sasamp']
            
        internal_kwargs={"--threads":_threads,"-t":transcriptome,"-i":index_path}
        #read build parameters
        yamlfile=os.path.join(_params_dir,'index.yaml')
        if pu.check_files_exist(yamlfile):
            yaml_params=pl.YAML_loader(yamlfile)
            yaml_kwargs=yaml_params.get_kwargs()
            internal_kwargs={**yaml_kwargs,**internal_kwargs}
            
        salmon_cmd=['salmon','index']
        salmon_cmd.extend(pu.parse_unix_args(validArgsIndex,internal_kwargs))
        
        #call salmon
        status=pe.execute_command(salmon_cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        
        if status:
            if pu.check_salmonindex(index_path) and not _dryrun:
                #update object's index
                self.index=index_path
                if self.check_index():
                    return True
        else:
            raise Exception("Error building salmon index")
        
        return False
        
        
        
    
    def perform_quant(self,sra_object,out_suffix="",out_dir="",verbose=False,quiet=False,logs=True,objectid="NA"):
        """run salmon quant
        sra_object: SRA
            An SRA object with valid fastq files
        lib_type: str
            Library type. Default:A
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
            Options to pass to salmon. This will override the existing options

        :return: Path to salmon out directory
        :rtype: string
        """
            
        if not out_dir:
            out_dir=os.path.join(sra_object.directory,"salmon_out")
        else:
            #create out_dir if not exists
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
        
        
        if sra_object.layout == 'PAIRED':
            internal_kwargs={"--threads":_threads,"-o":out_dir,"-l":"A","-1":sra_object.fastq_path,"-2":sra_object.fastq2_path,"-i":self.index}
        else:
            internal_kwargs={"--threads":_threads,"-o":out_dir,"-l":"A","-r":sra_object.fastq_path,"-i":self.index}
        
        #call salmon
        status=self.run(None,verbose=verbose,quiet=quiet,logs=logs,objectid=sra_object.srr_accession,**internal_kwargs)
        
        if status:
            outfile=os.path.join(out_dir,"quant.sf")
            if not pu.check_files_exist(outfile) and not _dryrun:
                return ""
            #move
            newfile=os.path.join(out_dir,"quant"+out_suffix+".sf")
            pe.move_file(outfile,newfile)
            return newfile
        
        return ""
        
        
        
    def run(self,*args,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper for running salmon.
        
        Parameters
        ----------
        
        subcommand: str
            subcommand for salmon
        valid_args: list
            List of valid arguments
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to salmon. This will override the existing options

        :return: Returns the status of salmon. True is passed, False if failed.
        :rtype: bool
        """
        
        #check for a valid index
        if not self.check_index() and not _dryrun:
            raise Exception("ERROR: Invalid salmon index. Please run build index to generate an index.")
        
        #override class kwargs by passed
        kwargs={**self._kwargs,**kwargs}
            
        salmon_Cmd=['salmon','quant']
        salmon_Cmd.extend(pu.parse_unix_args(valid_args._args_SALMON,kwargs))
        
        #start ececution
        status=pe.execute_command(salmon_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,command_name=" ".join(salmon_Cmd[0:2]))
        if not status:
            pu.print_boldred("salmon failed: "+" ".join(salmon_Cmd))
        return status 

    def check_index(self):
        if hasattr(self,'index'):
            return pu.check_salmonindex(self.index)
        return False
    
    
