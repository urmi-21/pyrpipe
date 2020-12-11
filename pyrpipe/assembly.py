#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:21:01 2019

@author: usingh
"""

from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe
import os
from pyrpipe import valid_args
from pyrpipe import param_loader as pl
from pyrpipe import _dryrun
from pyrpipe import _threads
from pyrpipe import _params_dir

class Assembly:
    """This class represents an abstract parent class for all programs which can perfrom transcripts assembly.
    """
    def __init__(self):
        self.category="Assembler"
    def perform_assembly(bam_file):
        """Function to perform assembly using a bam file as input. Inherited by all children.
        
                
        :param bam_file: path to input BAM
        :type bam_file: string
        
        :return: path to output GTF or output directory depending on the specific assembly program.
        :rtype: string
        """
        pass

class Stringtie(Assembly):
    """This class represents Stringtie program for transcript assembly.
        
        Parameters
        ----------
        
        threads: int
            number of threads
        """
    def __init__(self,*args,**kwargs):
        
        super().__init__()
        self.program_name="stringtie"
        #check if stringtie exists
        if not pe.check_dependencies([self.program_name]):
            raise Exception("ERROR: "+ self.program_name+" not found.")
        
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
        yamlfile=os.path.join(_params_dir,'stringtie.yaml')
        #override yaml parameters by kwargs
        if pu.check_files_exist(yamlfile):
            yaml_params=pl.YAML_loader(yamlfile)
            yaml_kwargs=yaml_params.get_kwargs()
            self._kwargs={**yaml_kwargs,**self._kwargs}

                         
    def perform_assembly(self,bam_file,out_dir=None,out_suffix="_stringtie",verbose=False,quiet=False,logs=True,objectid="NA"):
        """Function to run stringtie using a bam file.
                
        Parameters
        ----------
        
        bam_file: string
            path to the bam file
        out_suffix: string
            Suffix for the output gtf file
        reference_gtf: str
            Path to the reference gtf used as guide
        threads: int
            Number of threads to use
        overwrite: bool
            Overwrite if output file already exists.
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to stringtie. 
        :return: Returns the path to output GTF file
        :rtype: string
        """
        
        #create path to output file
        fname=pu.get_file_basename(bam_file)
        
        if not out_dir:
            out_dir=pu.get_file_directory(bam_file)
        
        if not pu.check_paths_exist(out_dir):
            pu.mkdir(out_dir)
            
        out_gtf_file=os.path.join(out_dir,fname+out_suffix+".gtf")
        
        
        #Add output file name and input bam
        internal_args=(bam_file,)
        internal_kwargs={"-o":out_gtf_file,"-p":_threads}
        #add positional args
        internal_kwargs['--']=internal_args
        
        
        #call stringtie
        status=self.run(verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**internal_kwargs)
        
        if status:
            #check if sam file is present in the location directory of sraOb
            if not pu.check_files_exist(out_gtf_file) and not _dryrun:
                return ""
            return out_gtf_file
        
        return ""        
            
    
    def run(self,*args,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper for running stringtie. This can be used to run stringtie without using perform_assembly() function.
        
        Parameters
        ----------
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
            Options to pass to stringtie. 
            
        :return: Returns the status of stringtie command.
        :rtype: bool
        """
            
        #override class kwargs by passed
        kwargs={**self._kwargs,**kwargs}
        
       
        stie_cmd=['stringtie']
        #add options
        stie_cmd.extend(pu.parse_unix_args(valid_args._args_STRINGTIE,kwargs))        
        
                
        #start ececution
        status=pe.execute_command(stie_cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            pu.print_boldred("stringtie failed: "+ " ".join(stie_cmd))
        
        #return status
        return status
    
    
    
class Cufflinks(Assembly):
    """This class represents cufflinks
    
    threads: int
            Number of threads to use
    """
    def __init__(self,*args,**kwargs):
        
        super().__init__()
        self.program_name="cufflinks"
        #check if stringtie exists
        if not pe.check_dependencies([self.program_name]):
            raise Exception("ERROR: "+ self.program_name+" not found.")
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
        yamlfile=os.path.join(_params_dir,'cufflinks.yaml')
        #override yaml parameters by kwargs
        if pu.check_files_exist(yamlfile):
            yaml_params=pl.YAML_loader(yamlfile)
            yaml_kwargs=yaml_params.get_kwargs()
            self._kwargs={**yaml_kwargs,**self._kwargs}       
        
    
    
    def perform_assembly(self,bam_file,out_dir=None,out_suffix="_cufflinks",verbose=False,quiet=False,logs=True,objectid="NA"):
        """Function to run cufflinks with BAM file as input.
                
        Parameters
        ----------
        bam_file: string
            path to bam file
        out_dir: 
            output directory
        out_suffix: string
            Suffix for the output gtf file
        reference_gtf: str
            Path to reference gtf 
        threads: int
            Number of threads to use
        overwrite: bool
            Overwrite if output file already exists.
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession.
        kwargs: dict
            Options to pass to cufflinks. 
            
        :return: Returns the path to output GTF file
        :rtype: string       
        """
        
        #create path to output file
        fname=pu.get_file_basename(bam_file)
        if not out_dir:
            out_dir=pu.get_file_directory(bam_file)
        else:
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
        out_gtf_file=os.path.join(out_dir,fname+out_suffix+".gtf")
        
            
        #Add output file name and input bam
        internal_args=(bam_file,)
        internal_kwargs={"-o":out_dir,"-p":_threads}
        #add positional args
        internal_kwargs['--']=internal_args
        
        #call cufflinks
        status=self.run(verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**internal_kwargs)
        
        
        if status:
            outfile=os.path.join(out_dir,"transcripts.gtf")
            #if file exist
            if not pu.check_files_exist(outfile) and not _dryrun:
                return ""
            
            #move out_dir/transcripts.gtf to outfile
            pe.move_file(outfile,out_gtf_file)
            if pu.check_files_exist(out_gtf_file):
                return out_gtf_file
            else: return ""
        return ""
    
      
    
    def run(self,*args,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper for running cufflinks
        
        Parameters
        ----------
        valid_args: list
            Liast of valid arguments
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession.
        kwargs: dict
            Options passed to cufflinks
        
        :return: Returns the status of cufflinks command.
        :rtype: bool
        """
            
        #override class kwargs by passed
        kwargs={**self._kwargs,**kwargs}
        
        cufflinks_cmd=['cufflinks']
        #add options
        cufflinks_cmd.extend(pu.parse_unix_args(valid_args._args_CUFFLINKS,kwargs))        
        
        
        #start ececution
        status=pe.execute_command(cufflinks_cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            pu.print_boldred("cufflinks failed: "+" ".join(cufflinks_cmd))
        #return status
        return status
    
    
    
class Trinity(Assembly):
    """This class represents Trinity RNA-Seq assembler
    
    threads: int
            Number of threads to use
    max_memory: int
            Max memory in GB to use
    """
    def __init__(self,*args,**kwargs):
        
        super().__init__()
        self.program_name="Trinity"
        self.dep_list=[self.program_name,'jellyfish','bowtie2']
        #check if trinity exists
        if not pe.check_dependencies(self.dep_list):
            raise Exception("ERROR: "+ self.program_name+" not found.")
        
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
        yamlfile=os.path.join(_params_dir,'trinity.yaml')
        #override yaml parameters by kwargs
        if pu.check_files_exist(yamlfile):
            yaml_params=pl.YAML_loader(yamlfile)
            yaml_kwargs=yaml_params.get_kwargs()
            self._kwargs={**yaml_kwargs,**self._kwargs}
        
    
    
    def perform_assembly(self,bam_file,out_dir=None,out_suffix="_trinity",verbose=False,quiet=False,logs=True,objectid="NA"):
        """Function to run trinity with sra object or BAM file as input.
                
        Parameters
        ----------
        
        sra_object: SRA
            object of SRA class
        bam_file: string
            path to bam file
        out_dir: string
            path to out directory
        max_memory: string
            Max memory argument e.g. "2G"
        max_intron: int
            specify the "--genome_guided_max_intron" argument
        threads: int
            Number of threads to use
        overwrite: bool
            Overwrite if output file already exists
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession.
        
        kwargs: dict
            Options to pass to Trinity. 
            
        :return: Return the path to output GTF file
        :rtype: string
        """
        
        #add trinity to outdir
        if "trinity" not in out_dir:
            out_dir+="_trinity"
            
        new_opts={}
        if sra_object is not None:
            parent_dir=sra_object.directory
            out_dir=os.path.join(parent_dir,out_dir)
            if sra_object.layout == 'PAIRED':
                new_opts={"--seqType":"fq","--left":sra_object.fastq_path,"--right":sra_object.fastq2_path,"--output":out_dir,"--CPU":_threads}
            else:
                new_opts={"--seqType":"fq","--single":sra_object.fastq_path,"--output":out_dir,"--CPU":_threads}
        elif bam_file is not None:
            if not pu.check_files_exist(bam_file):
                pu.print_boldred("Input to trinity does not exist:"+bam_file)
                return ""
            parent_dir=pu.get_file_directory(bam_file)
            out_dir=os.path.join(parent_dir,out_dir)
            new_opts={"--genome_guided_bam":bam_file,"--output":out_dir,"--max_memory":str(max_memory)+"G","--CPU":str(threads),"--genome_guided_max_intron":max_intron}
        else:
            pu.print_boldred("Please provide valid input to run trinity")
            return ""
        
        merged_opts={**new_opts,**kwargs}
        
        #call trinity
        status=self.run_trinity(valid_args_list=None,verbose=False,quiet=False,logs=True,objectid="NA",**merged_opts)
        
        if status:
            #check out dir
            if pu.check_paths_exist(out_dir):
                return out_dir
        else:
            return ""

    
    
    def run(self,valid_args_list=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper for running trinity
        
        Parameters
        ----------
        valid_args: list
            list of valid arguments
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession.
        kwargs: dict
            Options passed to trinity

        :return: Return the status of trinity command.
        :rtype: bool
        """
            
               
        trinity_cmd=['Trinity']
        #add options
        trinity_cmd.extend(pu.parse_unix_args(valid_args_list,kwargs))        
        
        
        #start ececution
        status=pe.execute_command(trinity_cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            pu.print_boldred("trinity failed")
        #return status
        return status
    
    
    
    
